#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <optional>
#include <ostream>
#include <unistd.h>

#include "QASMTransPrimitives.hpp"
#include "circuit_passes/transpiler.hpp"
#include "dump_pulses.hpp"
#include "dump_qasm.hpp"
#include "IR/chip.hpp"
#include "parser/qasm_parser.hpp"
#include <pybind11/embed.h>

namespace py = pybind11;
namespace fs = std::filesystem;

// Module-local globals for shared headers.
namespace QASMTrans
{
    std::unordered_set<std::string> g_device_basis_gates;
    std::unordered_map<std::string, std::string> g_merged_gate_aliases;
}

struct TranspileOptions
{
    std::string mode = "ibmq";
    std::string backend_config;       // Path to backend JSON (required)
    std::string pulse_template_path;  // Optional pulse template JSON
    std::string output_path;          // Optional QASM output path (file or directory)
    std::string pulse_output_path;    // Optional pulse JSON output path
    bool emit_run = false;            // Placeholder for future QICK wiring
    bool limited_qubits = false;
    bool disable_mapomatic = false;
    bool full_fidelity = false;
    bool allow_parameterized_merge = true;
    bool optimize_1q = false;
    bool optimize_2q_cancel = false;
    bool optimize_commute_2q = false;
    bool optimize_2q_synth = false;
    bool sabre_layout = false;
    bool fast_quality_routing = false;
    std::string routing_mode;    // "sabre" or "fast" (empty uses fast_quality_routing)
    bool routing_decay = false;
    double routing_decay_increment = 0.001;
    long long routing_decay_reset = 5;
    long long routing_trials = 1;
    std::size_t mapomatic_limit = 1000;
    std::size_t fast_quality_max_embeddings = 50;
    long long seed = -1; // Deterministic routing seed; <0 uses random device
    int verbose = 0;
};

struct TranspileResult
{
    std::string output_qasm;      // Transpiled OpenQASM text
    std::string output_qasm_path; // Where the QASM was written
    std::string pulse_schedule;   // JSON pulse schedule (when available)
    py::object pulse_doc;         // Parsed pulse JSON (None if not produced)
    std::string pulse_path;       // Where the pulse JSON was written (if any)
    std::string log;              // Aggregated log/output
    std::vector<IdxType> logical_to_physical;  // Logical -> physical mapping after routing
    std::vector<IdxType> measurement_mapping;  // Mapping used for final measurement ordering

    TranspileResult()
        : pulse_doc(py::none())
    {
    }
};

namespace
{
    using QASMTrans::g_device_basis_gates;
    using QASMTrans::g_merged_gate_aliases;

    struct StreamCapture
    {
        std::ostringstream out;
        std::ostringstream err;
        std::streambuf *old_out = nullptr;
        std::streambuf *old_err = nullptr;

        void start()
        {
            old_out = std::cout.rdbuf(out.rdbuf());
            old_err = std::cerr.rdbuf(err.rdbuf());
        }

        void stop()
        {
            if (old_out)
            {
                std::cout.rdbuf(old_out);
            }
            if (old_err)
            {
                std::cerr.rdbuf(old_err);
            }
        }
    };

    IdxType mode_from_string(std::string mode_name)
    {
        auto to_lower = [](std::string value)
        {
            std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
                           { return static_cast<char>(std::tolower(c)); });
            return value;
        };
        std::string lowered = to_lower(std::move(mode_name));
        if (lowered == "ibmq")
            return 0;
        if (lowered == "ionq")
            return 1;
        if (lowered == "quantinuum")
            return 2;
        if (lowered == "rigetti")
            return 3;
        if (lowered == "quafu")
            return 4;
        if (lowered == "iqm")
            return 5;
        throw std::invalid_argument("Unknown mode '" + lowered + "'");
    }

    void ingest_backend_metadata(const std::string &backend_path)
    {
        g_device_basis_gates.clear();
        g_merged_gate_aliases.clear();
        try
        {
            std::ifstream backend_stream(backend_path);
            if (!backend_stream.is_open())
            {
                throw std::runtime_error("Unable to open backend config at " + backend_path);
            }
            nlohmann::json backend_config = nlohmann::json::parse(backend_stream, nullptr, true, true);
            auto to_lower = [](std::string value)
            {
                std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
                               { return static_cast<char>(std::tolower(c)); });
                return value;
            };
            auto ingest_aliases = [&](const nlohmann::json &aliases)
            {
                if (!aliases.is_object())
                {
                    return;
                }
                for (const auto &item : aliases.items())
                {
                    std::string alias_name = to_lower(item.key());
                    if (alias_name.empty())
                    {
                        continue;
                    }
                    g_device_basis_gates.insert(alias_name);
                    const nlohmann::json &info = item.value();
                    if (!info.is_object())
                    {
                        continue;
                    }
                    std::string logical = to_lower(info.value("logical_gate", std::string{}));
                    if (logical.empty())
                    {
                        continue;
                    }
                    if (g_merged_gate_aliases.find(logical) == g_merged_gate_aliases.end())
                    {
                        g_merged_gate_aliases[logical] = alias_name;
                    }
                }
            };
            auto ingest_basis = [&](const nlohmann::json &arr)
            {
                if (!arr.is_array())
                {
                    return;
                }
                for (const auto &entry : arr)
                {
                    if (entry.is_string())
                    {
                        std::string gate = to_lower(entry.get<std::string>());
                        if (!gate.empty())
                        {
                            g_device_basis_gates.insert(gate);
                        }
                    }
                }
            };
            auto ingest_gate_lens = [&](const nlohmann::json &obj)
            {
                if (!obj.is_object())
                {
                    return;
                }
                for (auto it = obj.begin(); it != obj.end(); ++it)
                {
                    const std::string key = it.key();
                    size_t first_digit = key.find_first_of("0123456789");
                    if (first_digit == std::string::npos)
                    {
                        continue;
                    }
                    std::string gate = to_lower(key.substr(0, first_digit));
                    if (!gate.empty())
                    {
                        g_device_basis_gates.insert(gate);
                    }
                }
            };
            ingest_basis(backend_config.value("basis_gates", nlohmann::json::array()));
            if (backend_config.contains("metadata") && backend_config["metadata"].is_object())
            {
                ingest_basis(backend_config["metadata"].value("basis_gates", nlohmann::json::array()));
                ingest_aliases(backend_config["metadata"].value("merged_gate_aliases", nlohmann::json::object()));
            }
            ingest_aliases(backend_config.value("merged_gate_aliases", nlohmann::json::object()));
            if (g_device_basis_gates.empty())
            {
                ingest_gate_lens(backend_config.value("gate_lens", nlohmann::json::object()));
            }
        }
        catch (const std::exception &ex)
        {
            throw std::runtime_error(std::string("Failed to parse backend config: ") + ex.what());
        }
    }

    std::vector<IdxType> build_measurement_mapping(const std::map<std::string, creg> &cregs,
                                                   const std::vector<IdxType> &logical_to_physical)
    {
        std::vector<IdxType> mapping;
        for (const auto &entry : cregs)
        {
            const auto &qubit_indices = entry.second.qubit_indices;
            for (std::size_t pos = 0; pos < qubit_indices.size(); ++pos)
            {
                IdxType logical_index = qubit_indices[pos];
                IdxType physical = 0;
                if (logical_index == UN_DEF || logical_index < 0)
                {
                    mapping.push_back(UN_DEF);
                    continue;
                }
                if (logical_index >= static_cast<IdxType>(logical_to_physical.size()))
                {
                    throw std::runtime_error("Invalid measurement mapping for creg '" +
                                             entry.first + "' (logical index " +
                                             std::to_string(logical_index) + ")");
                }
                if (logical_index >= 0)
                {
                    physical = logical_to_physical[static_cast<std::size_t>(logical_index)];
                }
                mapping.push_back(physical);
            }
        }
        if (mapping.empty())
        {
            return logical_to_physical;
        }
        return mapping;
    }

    std::string derive_pulse_output_path(const std::string &qasm_output_path,
                                         const std::string &supplied)
    {
        if (!supplied.empty())
        {
            return supplied;
        }
        fs::path qasm_path(qasm_output_path);
        fs::path directory = qasm_path.parent_path();
        std::string stem = qasm_path.stem().string();
        if (stem.empty())
        {
            stem = qasm_path.filename().string();
        }
        fs::path candidate = directory / (stem + "_pulses.json");
        return candidate.string();
    }

    TranspileResult transpile_qasm(const std::string &qasm_source, const TranspileOptions &options)
    {
        if (options.backend_config.empty())
        {
            throw std::invalid_argument("backend_config is required");
        }

        IdxType mode = mode_from_string(options.mode);
        ingest_backend_metadata(options.backend_config);

        fs::path input_path;
        bool created_temp = false;
        const bool looks_like_path = qasm_source.find('\n') == std::string::npos &&
                                     qasm_source.find('\r') == std::string::npos &&
                                     fs::exists(qasm_source);
        if (looks_like_path)
        {
            input_path = fs::absolute(qasm_source);
        }
        else
        {
            fs::path temp_dir = fs::temp_directory_path();
            fs::path temp_file;
            // Create a pseudo-unique filename without relying on unique_path (C++17 portability).
            {
                auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
                auto pid = static_cast<long long>(::getpid());
                std::ostringstream oss;
                oss << "qasmtrans_input_" << pid << "_" << now << ".qasm";
                temp_file = temp_dir / fs::path(oss.str());
            }
            std::ofstream temp_out(temp_file);
            if (!temp_out.is_open())
            {
                throw std::runtime_error("Failed to create temporary QASM file at " + temp_file.string());
            }
            temp_out << qasm_source;
            temp_out.close();
            input_path = temp_file;
            created_temp = true;
        }

        TranspileResult result;
        std::ostringstream log;
        StreamCapture capture;
        capture.start();
        try
        {
            qasm_parser parser(input_path.c_str());
            IdxType n_qubits = parser.num_qubits();
            auto circuit = std::make_shared<Circuit>(n_qubits);
            parser.loadin_circuit(circuit);
            std::map<std::string, creg> cregs = parser.get_list_cregs();

            auto chip = constructChip(n_qubits, options.backend_config, options.limited_qubits, options.verbose);
            if (!chip)
            {
                throw std::runtime_error("Failed to construct device from backend config");
            }

            CriticalPathHeuristicMode cp_mode = CriticalPathHeuristicMode::LogProduct;
            std::optional<uint64_t> routing_seed;
            if (options.seed >= 0)
            {
                routing_seed = static_cast<uint64_t>(options.seed);
            }
            RoutingMode routing_mode = RoutingMode::Sabre;
            if (!options.routing_mode.empty())
            {
                std::string mode_value = options.routing_mode;
                std::transform(mode_value.begin(), mode_value.end(), mode_value.begin(),
                               [](unsigned char ch)
                               { return static_cast<char>(std::tolower(ch)); });
                if (mode_value == "sabre")
                {
                    routing_mode = RoutingMode::Sabre;
                }
                else if (mode_value == "fast")
                {
                    routing_mode = RoutingMode::FastQuality;
                }
                else
                {
                    throw std::invalid_argument("Unknown routing_mode '" + options.routing_mode + "'. Use 'sabre' or 'fast'.");
                }
            }
            else if (options.fast_quality_routing)
            {
                routing_mode = RoutingMode::FastQuality;
            }
            transpiler(circuit,
                       chip,
                       cregs,
                       options.verbose,
                       mode,
                       options.full_fidelity,
                       cp_mode,
                       options.disable_mapomatic,
                       options.mapomatic_limit,
                       options.optimize_1q,
                       options.optimize_2q_cancel,
                       options.optimize_commute_2q,
                       options.optimize_2q_synth,
                       options.routing_decay,
                       options.routing_decay_increment,
                       static_cast<IdxType>(options.routing_decay_reset),
                       static_cast<IdxType>(options.routing_trials),
                       options.sabre_layout,
                       routing_mode,
                       options.fast_quality_max_embeddings,
                       g_device_basis_gates,
                       routing_seed);

            result.logical_to_physical = circuit->get_mapping();

            // Build measurement mapping if needed (single-circuit case).
            std::vector<IdxType> measurement = build_measurement_mapping(cregs, circuit->get_mapping());
            if (!measurement.empty())
            {
                circuit->set_mapping(measurement);
                result.measurement_mapping = measurement;
            }

            std::vector<QASMTrans::Gate> expanded_gates;
            bool rigetti_mode = false;
            if (!options.pulse_template_path.empty())
            {
                expanded_gates = QASMTrans::pulses::expandGatesForPulseDump(circuit,
                                                                             options.backend_config,
                                                                             options.pulse_template_path,
                                                                             options.allow_parameterized_merge,
                                                                             &rigetti_mode);
            }
            std::string qasm_path;
            if (!expanded_gates.empty() && rigetti_mode)
            {
                qasm_path = dumpQASMFromGates(circuit,
                                              expanded_gates,
                                              input_path.filename().c_str(),
                                              options.output_path,
                                              options.verbose,
                                              mode);
            }
            else
            {
                qasm_path = dumpQASM(circuit,
                                     input_path.filename().c_str(),
                                     options.output_path,
                                     options.verbose,
                                     mode);
            }
            result.output_qasm_path = qasm_path;
            {
                std::ifstream qasm_in(qasm_path);
                std::stringstream buffer;
                buffer << qasm_in.rdbuf();
                result.output_qasm = buffer.str();
            }

            if (!options.pulse_template_path.empty())
            {
                std::string pulse_path = derive_pulse_output_path(qasm_path, options.pulse_output_path);
                dumpPulses(circuit,
                           input_path.filename().c_str(),
                           options.backend_config,
                           options.pulse_template_path,
                           pulse_path,
                           options.verbose,
                           options.allow_parameterized_merge);
                result.pulse_path = pulse_path;
                std::ifstream pulse_in(pulse_path);
                std::stringstream buffer;
                buffer << pulse_in.rdbuf();
                result.pulse_schedule = buffer.str();
                try
                {
                    result.pulse_doc = py::module::import("json").attr("loads")(result.pulse_schedule);
                }
                catch (const std::exception &)
                {
                    result.pulse_doc = py::none();
                }
            }

            log << "Transpiled " << input_path.filename().string() << " with mode=" << options.mode
                << ", qubits=" << n_qubits;
        }
        catch (...)
        {
            capture.stop();
            if (created_temp)
            {
                std::error_code ec;
                fs::remove(input_path, ec);
            }
            throw;
        }
        capture.stop();

        if (created_temp)
        {
            std::error_code ec;
            fs::remove(input_path, ec);
        }
        result.log = log.str();
        std::string captured = capture.out.str();
        if (!captured.empty())
        {
            result.log += "\n" + captured;
        }
        return result;
    }
} // namespace

PYBIND11_MODULE(qasmtrans_core, m)
{
    m.doc() = "Prototype Python bindings for QASMTrans.";

    py::class_<TranspileOptions>(m, "TranspileOptions", "Configuration options for the transpiler.")
        .def(py::init<>())
        .def_readwrite("mode", &TranspileOptions::mode)
        .def_readwrite("backend_config", &TranspileOptions::backend_config)
        .def_readwrite("pulse_template_path", &TranspileOptions::pulse_template_path)
        .def_readwrite("output_path", &TranspileOptions::output_path)
        .def_readwrite("pulse_output_path", &TranspileOptions::pulse_output_path)
        .def_readwrite("emit_run", &TranspileOptions::emit_run)
        .def_readwrite("limited_qubits", &TranspileOptions::limited_qubits)
        .def_readwrite("disable_mapomatic", &TranspileOptions::disable_mapomatic)
        .def_readwrite("full_fidelity", &TranspileOptions::full_fidelity)
        .def_readwrite("allow_parameterized_merge", &TranspileOptions::allow_parameterized_merge)
        .def_readwrite("optimize_1q", &TranspileOptions::optimize_1q)
        .def_readwrite("optimize_2q_cancel", &TranspileOptions::optimize_2q_cancel)
        .def_readwrite("optimize_commute_2q", &TranspileOptions::optimize_commute_2q)
        .def_readwrite("optimize_2q_synth", &TranspileOptions::optimize_2q_synth)
        .def_readwrite("sabre_layout", &TranspileOptions::sabre_layout)
        .def_readwrite("fast_quality_routing", &TranspileOptions::fast_quality_routing)
        .def_readwrite("routing_mode", &TranspileOptions::routing_mode)
        .def_readwrite("routing_decay", &TranspileOptions::routing_decay)
        .def_readwrite("routing_decay_increment", &TranspileOptions::routing_decay_increment)
        .def_readwrite("routing_decay_reset", &TranspileOptions::routing_decay_reset)
        .def_readwrite("routing_trials", &TranspileOptions::routing_trials)
        .def_readwrite("mapomatic_limit", &TranspileOptions::mapomatic_limit)
        .def_readwrite("fast_quality_max_embeddings", &TranspileOptions::fast_quality_max_embeddings)
        .def_readwrite("seed", &TranspileOptions::seed)
        .def_readwrite("verbose", &TranspileOptions::verbose);

    py::class_<TranspileResult>(m, "TranspileResult", "Outputs from a transpilation run.")
        .def(py::init<>())
        .def_readonly("output_qasm", &TranspileResult::output_qasm)
        .def_readonly("output_qasm_path", &TranspileResult::output_qasm_path)
        .def_readonly("pulse_schedule", &TranspileResult::pulse_schedule)
        .def_readonly("pulse_doc", &TranspileResult::pulse_doc)
        .def_readonly("pulse_path", &TranspileResult::pulse_path)
        .def_readonly("log", &TranspileResult::log)
        .def_readonly("logical_to_physical", &TranspileResult::logical_to_physical)
        .def_readonly("measurement_mapping", &TranspileResult::measurement_mapping);

    m.def(
        "transpile_qasm",
        &transpile_qasm,
        py::arg("qasm_source"),
        py::arg("options") = TranspileOptions{},
        R"pbdoc(
            Transpile an OpenQASM circuit (file path or raw string).

            Parameters
            ----------
            qasm_source:
                Either a path to a QASM file or the QASM program text.
            options:
                TranspileOptions controlling backend mode, qubit limits, and optional pulse emission paths.
        )pbdoc");
    // QICK utilities
    m.def(
        "load_qick_config",
        [](const std::string &path)
        {
            py::module json = py::module::import("json");
            py::object loader = json.attr("load");
            std::ifstream in(path);
            if (!in.is_open())
            {
                throw std::runtime_error("Unable to open QICK config at " + path);
            }
            py::object py_file = py::cast(in, py::return_value_policy::reference);
            return loader(py_file);
        },
        py::arg("path"),
        "Load a QICK config JSON into a Python dict.");

    m.def(
        "emit_qick",
        [](py::object pulse_doc, py::object qick_config, bool run_enabled, bool summary_only)
        {
            py::object fn;
            try
            {
                py::module runner = py::module::import("qasmtrans.qick_pulse_runner");
                fn = runner.attr("run_program_obj");
            }
            catch (const py::error_already_set &)
            {
                py::module runner = py::module::import("qick_pulse_runner");
                fn = runner.attr("run_program_obj");
            }
            return fn(pulse_doc, qick_config, run_enabled, summary_only);
        },
        py::arg("pulse_doc"),
        py::arg("qick_config"),
        py::arg("run_enabled") = false,
        py::arg("summary_only") = false,
        R"pbdoc(
            Emit a pulse schedule to QICK using in-memory documents.

            pulse_doc: dict parsed from a QASMTrans pulse JSON (e.g., result.pulse_doc).
            qick_config: dict with qubit_gen_map and generator/readout settings.
            run_enabled: if True, execute on hardware; otherwise only prepare summary.
            summary_only: if True, skip hardware execution even if run_enabled is True.
        )pbdoc");
}
