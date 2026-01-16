#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../include/QASMTransPrimitives.hpp"
#include "../include/cli_support.hpp"
#include "../include/IR/chip.hpp"
#include "../include/circuit_passes/transpiler.hpp"
#include "../include/dump_pulses.hpp"
#include "../include/parser/parser_util.hpp"
#include "../include/parser/qasm_parser.hpp"
#include "../include/qick_emitter.hpp"
#include "../include/util/chip_partition.hpp"

using namespace QASMTrans;
namespace fs = std::filesystem;

namespace QASMTrans
{
    std::unordered_set<std::string> g_device_basis_gates;
    std::unordered_map<std::string, std::string> g_merged_gate_aliases;
}

namespace
{
    struct SubchipArtifact
    {
        // Snapshot of per-circuit data needed to export standalone subchip artifacts later.
        std::size_t index = 0;
        std::string circuit_prefix;
        std::string source_filename;
        std::shared_ptr<Chip> subchip;
        std::vector<IdxType> local_to_global;
        std::vector<QASMTrans::Gate> local_gates;
        std::vector<IdxType> logical_mapping;
        std::vector<IdxType> measurement_mapping;
        std::map<std::string, creg> cregs;
    };

    std::vector<IdxType> build_measurement_mapping(const std::map<std::string, creg> &cregs,
                                                   const std::vector<IdxType> &logical_to_physical)
    {
        // Classical bits are stored by creg; map them to the physical qubit order we ended up using.
        std::vector<IdxType> mapping;
        for (const auto &entry : cregs)
        {
            const auto &qubit_indices = entry.second.qubit_indices;
            for (std::size_t pos = 0; pos < qubit_indices.size(); ++pos)
            {
                IdxType logical_index = qubit_indices[pos];
                IdxType physical = 0;
                if (logical_index != UN_DEF && logical_index >= 0 &&
                    logical_index < static_cast<IdxType>(logical_to_physical.size()))
                {
                    physical = logical_to_physical[static_cast<std::size_t>(logical_index)];
                }
                else if (!logical_to_physical.empty())
                {
                    physical = logical_to_physical[pos % logical_to_physical.size()];
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

    void print_help()
    {
        // Minimal usage text; aligns with legacy CLI expectations.
        std::cout << "Usage: ./qasmtrans [options]" << std::endl;
        std::cout << "Option            Description" << std::endl;
        std::cout << "-i <path>         Input QASM file (repeat -i for multiple circuits)" << std::endl;
        std::cout << "-c <backend>      Path to backend configuration json file" << std::endl;
        std::cout << "-limited          Run the transpiler with limited physical qubits usage" << std::endl;
        std::cout << "-limited          Limit qubit usage to circuit than device. "
                  << "It reduces qubit usage but may introduce extra routing cost or unable to route" << std::endl;
        std::cout << "-backend_list     Print the available device backends" << std::endl;
        std::cout << "-m <name>         Set the transpiler targeted device (ibmq, ionq, quantinuum, rigetti, quafu, iqm), default is ibmq" << std::endl;
        std::cout << "-v <0/1/2>        Set the output level, default is 0" << std::endl;
        std::cout << "-full_fidelity    Score Mapomatic candidates on the entire circuit instead of its critical path" << std::endl;
        std::cout << "-cp_mode <product|hybrid>  Choose scoring strategy (default product)" << std::endl;
        std::cout << "-mapomatic_limit <N>       Limit the number of candidate embeddings Mapomatic evaluates (default 1000)" << std::endl;
        std::cout << "--disable_mapomatic       Skip the calibration-aware Mapomatic pass" << std::endl;
        std::cout << "-o <path>         Set the output file, "
                  << "default is data/output/transpiled_modename_filename.qasm" << std::endl;
        std::cout << "-p <path>         Pulse template json (optional; enables pulse dumping)" << std::endl;
        std::cout << "-e <config>       Emit pulses via QICK using the provided QICK config" << std::endl;
        std::cout << "--emit-run        When paired with -e, stream pulses to hardware (otherwise summary only)" << std::endl;
        std::cout << "--merge-allow-params     Include parameterised logical gates as merge candidates (default)" << std::endl;
        std::cout << "--merge-disallow-params  Exclude parameterised logical gates from merge candidate analysis" << std::endl;
        std::cout << "-h                print the help function" << std::endl;
    }
} // namespace

using namespace QASMTrans::cli;

struct CliConfig
{
    bool run_with_limit = false;
    IdxType mode = 0;
    std::string mode_name = "ibmq";
    IdxType debug_level = 0;
    std::string output_path = "../data/output/";
    std::string pulse_template_path;
    bool emit_requested = false;
    bool emit_run = false;
    std::string qick_config_path;
    bool allow_parameterized_merge_candidates = true;
    std::string backendpath;
    std::vector<std::string> input_files;
    bool disable_mapomatic = false;
    bool use_full_fidelity = false;
    CriticalPathHeuristicMode cp_mode = CriticalPathHeuristicMode::LogProduct;
    std::size_t mapomatic_limit = 1000;
    std::string program_name = "qasmtrans";
};

struct CircuitBatch
{
    std::vector<std::shared_ptr<Circuit>> circuits;
    std::vector<std::map<std::string, creg>> circuit_cregs;
    std::vector<IdxType> circuit_sizes;
};

struct CombinedArtifacts
{
    std::vector<Gate> combined_gates;
    std::vector<IdxType> combined_mapping;
    std::map<std::string, creg> combined_cregs;
    std::vector<SubchipArtifact> subchip_artifacts;
};

struct OutputArtifacts
{
    std::string qasm_path;
    std::string pulses_path;
};

bool parse_cli(int argc, char **argv, CliConfig &config, int &exit_code)
{
    if (argc > 0 && argv && argv[0])
    {
        config.program_name = argv[0];
    }

    if (argc == 1)
    {
        print_help();
        exit_code = 0;
        return false;
    }
    if (cmdOptionExists(argv, argv + argc, "-h"))
    {
        print_help();
        exit_code = 0;
        return false;
    }

    if (cmdOptionExists(argv, argv + argc, "-limited"))
    {
        config.run_with_limit = true;
    }
    if (cmdOptionExists(argv, argv + argc, "-cp_mode"))
    {
        std::string cp_mode_value = std::string(getCmdOption(argv, argv + argc, "-cp_mode"));
        std::string lowered = cp_mode_value;
        std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                       [](unsigned char ch)
                       { return static_cast<char>(std::tolower(ch)); });
        if (lowered == "product" || lowered == "log" || lowered == "multiplicative")
        {
            config.cp_mode = CriticalPathHeuristicMode::LogProduct;
        }
        else if (lowered == "hybrid")
        {
            config.cp_mode = CriticalPathHeuristicMode::Hybrid;
        }
        else
        {
            std::cerr << "Error: unknown -cp_mode value '" << cp_mode_value
                      << "'. Expected 'product' or 'hybrid'." << std::endl;
            exit_code = 1;
            return false;
        }
    }
    if (cmdOptionExists(argv, argv + argc, "-mapomatic_limit"))
    {
        const char *opt = getCmdOption(argv, argv + argc, "-mapomatic_limit");
        if (!opt)
        {
            std::cerr << "Error: -mapomatic_limit requires a positive integer argument." << std::endl;
            exit_code = 1;
            return false;
        }
        try
        {
            long long parsed = std::stoll(opt);
            if (parsed <= 0)
            {
                std::cerr << "Error: -mapomatic_limit must be greater than zero (got " << parsed << ")." << std::endl;
                exit_code = 1;
                return false;
            }
            config.mapomatic_limit = static_cast<std::size_t>(parsed);
        }
        catch (const std::exception &)
        {
            std::cerr << "Error: failed to parse -mapomatic_limit argument '" << opt << "'." << std::endl;
            exit_code = 1;
            return false;
        }
    }
    if (cmdOptionExists(argv, argv + argc, "-v"))
    {
        config.debug_level = IdxType(std::stoi(getCmdOption(argv, argv + argc, "-v")));
    }
    if (cmdOptionExists(argv, argv + argc, "-full_fidelity"))
    {
        config.use_full_fidelity = true;
    }
    if (cmdOptionExists(argv, argv + argc, "-o"))
    {
        config.output_path = std::string(getCmdOption(argv, argv + argc, "-o"));
    }
    if (cmdOptionExists(argv, argv + argc, "-p"))
    {
        config.pulse_template_path = std::string(getCmdOption(argv, argv + argc, "-p"));
    }
    if (cmdOptionExists(argv, argv + argc, "-e"))
    {
        const char *emit_config = getCmdOption(argv, argv + argc, "-e");
        if (emit_config == nullptr)
        {
            std::cerr << "Error: missing QICK config path after -e" << std::endl;
            exit_code = 1;
            return false;
        }
        config.qick_config_path = std::string(emit_config);
        config.emit_requested = true;
    }
    if (cmdOptionExists(argv, argv + argc, "--emit-run"))
    {
        config.emit_run = true;
    }
    if (cmdOptionExists(argv, argv + argc, "--merge-disallow-params"))
    {
        config.allow_parameterized_merge_candidates = false;
    }
    if (cmdOptionExists(argv, argv + argc, "--merge-allow-params"))
    {
        config.allow_parameterized_merge_candidates = true;
    }

    // Manual scan to pick up repeatable flags like multiple -i inputs.
    int argi = 1;
    while (argi < argc)
    {
        std::string current = argv[argi];
        if (current == "-i")
        {
            if (argi + 1 >= argc)
            {
                std::cerr << "Error: -i requires a following QASM file path." << std::endl;
                exit_code = 1;
                return false;
            }
            config.input_files.emplace_back(argv[argi + 1]);
            argi += 2;
            continue;
        }
        if (current == "--disable_mapomatic")
        {
            config.disable_mapomatic = true;
            ++argi;
            continue;
        }
        if (current == "--merge-disallow-params")
        {
            config.allow_parameterized_merge_candidates = false;
            ++argi;
            continue;
        }
        if (current == "--merge-allow-params")
        {
            config.allow_parameterized_merge_candidates = true;
            ++argi;
            continue;
        }
        if (current == "--emit-run")
        {
            config.emit_run = true;
            ++argi;
            continue;
        }
        ++argi;
    }

    if (cmdOptionExists(argv, argv + argc, "-c"))
    {
        config.backendpath = std::string(getCmdOption(argv, argv + argc, "-c"));
    }
    if (cmdOptionExists(argv, argv + argc, "-backend_list"))
    {
        std::cout << "The available backends are:" << std::endl;
        std::cout << "ibmq_toronto (27 qubits)" << std::endl;
        std::cout << "ibmq_jakarta (7 qubits)" << std::endl;
        std::cout << "ibmq_guadalupe (16 qubits)" << std::endl;
        std::cout << "ibm_seattle (433 qubits)" << std::endl;
        std::cout << "ibm_cairo (27 qubits)" << std::endl;
        std::cout << "ibm_brisbane (127 qubits)" << std::endl;
        std::cout << "aspen_m3 (80 qubits)" << std::endl;
        std::cout << "h1_2 (12 qubits)" << std::endl;
        std::cout << "h1_1 (20 qubits)" << std::endl;
        std::cout << "dummy_ibmq12 (12 qubits)" << std::endl;
        std::cout << "dummy_ibmq14 (14 qubits)" << std::endl;
        std::cout << "dummy_ibmq15 (15 qubits)" << std::endl;
        std::cout << "dummy_ibmq16 (16 qubits)" << std::endl;
        std::cout << "dummy_ibmq30 (30 qubits)" << std::endl;
        std::cout << "You can manually add new machine in json file at data/device" << std::endl;
        exit_code = 0;
        return false;
    }
    if (config.input_files.empty())
    {
        std::cerr << "Error: missing input QASM file(s) (-i)." << std::endl;
        exit_code = 1;
        return false;
    }
    if (config.backendpath.empty())
    {
        std::cerr << "Error: missing machine backend file via -c" << std::endl;
        exit_code = 1;
        return false;
    }
    if (config.emit_requested && config.pulse_template_path.empty())
    {
        std::cerr << "Error: -e requires a pulse template via -p to generate pulses." << std::endl;
        exit_code = 1;
        return false;
    }
    if (cmdOptionExists(argv, argv + argc, "-m"))
    {
        config.mode_name = std::string(getCmdOption(argv, argv + argc, "-m"));
        if (config.mode_name == "ibmq" || config.mode_name == "IBMQ")
        {
            config.mode = 0;
        }
        else if (config.mode_name == "ionq" || config.mode_name == "IonQ")
        {
            config.mode = 1;
        }
        else if (config.mode_name == "Quantinuum" || config.mode_name == "quantinuum")
        {
            config.mode = 2;
        }
        else if (config.mode_name == "Rigetti" || config.mode_name == "rigetti")
        {
            config.mode = 3;
        }
        else if (config.mode_name == "Quafu" || config.mode_name == "quafu")
        {
            config.mode = 4;
        }
        else if (config.mode_name == "IQM" || config.mode_name == "iqm")
        {
            config.mode = 5;
        }
        else
        {
            std::cout << "Invalid mode name, please check" << std::endl;
            exit_code = 0;
            return false;
        }
    }
    exit_code = 0;
    return true;
}

CircuitBatch load_circuits(const std::vector<std::string> &input_files)
{
    CircuitBatch batch;
    batch.circuits.reserve(input_files.size());
    batch.circuit_cregs.reserve(input_files.size());
    batch.circuit_sizes.reserve(input_files.size());

    // Parse each input QASM into its own circuit and collect metadata.
    for (const auto &file : input_files)
    {
        qasm_parser parser(file.c_str());
        IdxType n_qubits = parser.num_qubits();
        auto circuit = make_shared<Circuit>(n_qubits);
        parser.loadin_circuit(circuit);
        map<string, creg> cregs = parser.get_list_cregs();
        batch.circuits.push_back(circuit);
        batch.circuit_cregs.push_back(cregs);
        batch.circuit_sizes.push_back(n_qubits);
    }
    return batch;
}

std::vector<IdxType> compute_partition_sizes(const std::vector<IdxType> &circuit_sizes,
                                             IdxType device_capacity,
                                             IdxType total_requested_qubits)
{
    std::vector<IdxType> partition_sizes = circuit_sizes;
    if (partition_sizes.empty())
    {
        return partition_sizes;
    }

    // Allocate device capacity proportionally (with rounding) across requested circuits.
    const double total_logical = static_cast<double>(total_requested_qubits);
    std::vector<double> remainders(partition_sizes.size(), 0.0);
    IdxType allocated = 0;
    for (std::size_t idx = 0; idx < partition_sizes.size(); ++idx)
    {
        double exact = total_logical > 0.0
                           ? (static_cast<double>(partition_sizes[idx]) / total_logical) * static_cast<double>(device_capacity)
                           : static_cast<double>(partition_sizes[idx]);
        double base = std::floor(exact);
        IdxType candidate = static_cast<IdxType>(base);
        if (candidate < circuit_sizes[idx])
        {
            candidate = circuit_sizes[idx];
        }
        partition_sizes[idx] = candidate;
        remainders[idx] = exact - base;
        allocated += candidate;
    }

    if (allocated < device_capacity)
    {
        IdxType leftover = device_capacity - allocated;
        if (leftover > 0)
        {
            std::vector<std::size_t> order(partition_sizes.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(),
                      [&](std::size_t lhs, std::size_t rhs)
                      {
                          if (remainders[lhs] != remainders[rhs])
                          {
                              return remainders[lhs] > remainders[rhs];
                          }
                          return partition_sizes[lhs] < partition_sizes[rhs];
                      });
            std::size_t cursor = 0;
            while (leftover > 0 && !order.empty())
            {
                std::size_t idx = order[cursor % order.size()];
                partition_sizes[idx] += 1;
                ++cursor;
                --leftover;
            }
        }
    }
    else if (allocated > device_capacity)
    {
        IdxType overshoot = allocated - device_capacity;
        std::vector<std::size_t> order(partition_sizes.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](std::size_t lhs, std::size_t rhs)
                  {
                      IdxType slack_lhs = partition_sizes[lhs] - circuit_sizes[lhs];
                      IdxType slack_rhs = partition_sizes[rhs] - circuit_sizes[rhs];
                      if (slack_lhs != slack_rhs)
                      {
                          return slack_lhs > slack_rhs;
                      }
                      return remainders[lhs] < remainders[rhs];
                  });
        for (std::size_t pos = 0; overshoot > 0 && pos < order.size(); ++pos)
        {
            std::size_t idx = order[pos];
            while (overshoot > 0 && partition_sizes[idx] > circuit_sizes[idx])
            {
                partition_sizes[idx] -= 1;
                --overshoot;
            }
        }
        if (overshoot > 0)
        {
            throw std::runtime_error("Unable to fit requested circuit sizes within device capacity.");
        }
    }
    return partition_sizes;
}

std::vector<std::vector<IdxType>> allocate_partitions(const std::shared_ptr<Chip> &chip,
                                                      const std::vector<IdxType> &partition_sizes,
                                                      IdxType debug_level)
{
    std::vector<std::vector<IdxType>> partitions;
    auto partition_start = std::chrono::steady_clock::now();
    try
    {
        // Try to find well-distributed subchips; fall back to contiguous allocation.
        partitions = partition_chip(chip, partition_sizes);
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Warning: advanced partitioning failed (" << ex.what()
                  << "), falling back to contiguous allocation." << std::endl;
        partitions.clear();
        partitions.reserve(partition_sizes.size());
        IdxType next_qubit = 0;
        bool fallback_ok = true;
        for (IdxType requested : partition_sizes)
        {
            std::vector<IdxType> part;
            part.reserve(static_cast<std::size_t>(requested));
            for (IdxType offset = 0; offset < requested; ++offset)
            {
                if (next_qubit >= chip->chip_qubit_num)
                {
                    fallback_ok = false;
                    break;
                }
                part.push_back(next_qubit++);
            }
            if (static_cast<IdxType>(part.size()) != requested)
            {
                fallback_ok = false;
                break;
            }
            partitions.push_back(std::move(part));
        }
        if (!fallback_ok || partitions.size() != partition_sizes.size())
        {
            throw std::runtime_error("Error: contiguous fallback partitioning failed.");
        }
    }
    auto partition_end = std::chrono::steady_clock::now();
    auto partition_ms = std::chrono::duration_cast<std::chrono::milliseconds>(partition_end - partition_start).count();
    std::cout << "[timing] partition_ms=" << partition_ms << std::endl;
    if (debug_level > 0 && !partitions.empty())
    {
        std::cout << "Partitioned device into " << partitions.size() << " subchip(s)." << std::endl;
    }
    return partitions;
}

CombinedArtifacts transpile_and_merge(const std::vector<std::shared_ptr<Circuit>> &circuits,
                                      const std::vector<std::map<std::string, creg>> &circuit_cregs,
                                      const std::vector<std::vector<IdxType>> &partitions,
                                      const std::vector<std::string> &input_files,
                                      const std::shared_ptr<Chip> &chip,
                                      const CliConfig &config,
                                      IdxType total_requested_qubits)
{
    CombinedArtifacts artifacts;
    artifacts.combined_gates.reserve(1024);
    artifacts.combined_mapping.reserve(static_cast<std::size_t>(total_requested_qubits));
    artifacts.subchip_artifacts.reserve(circuits.size());

    // Keep creg prefixes aligned (e.g., circuit00_, circuit01_) for readability.
    const std::size_t prefix_width = circuits.empty() ? 1 : std::to_string(circuits.size() - 1).size();

    for (std::size_t ci = 0; ci < circuits.size(); ++ci)
    {
        // Bind each circuit to its subchip and run the transpiler.
        auto circuit = circuits[ci];
        const auto &cregs = circuit_cregs[ci];

        std::vector<IdxType> local_to_global;
        auto subchip = make_subchip(chip, partitions[ci], local_to_global);

        if (config.debug_level > 0)
        {
            std::cout << "Circuit " << ci << ": allocated " << partitions[ci].size()
                      << " physical qubits." << std::endl;
        }
        if (config.debug_level > 1)
        {
            std::cout << "  Classical registers:" << std::endl;
            for (const auto &entry : cregs)
            {
                std::cout << "    " << entry.first << " (width=" << entry.second.width << ") indices:";
                for (IdxType idx_val : entry.second.qubit_indices)
                {
                    std::cout << " " << idx_val;
                }
                std::cout << std::endl;
            }
        }

        transpiler(circuit,
                   subchip,
                   cregs,
                   config.debug_level,
                   config.mode,
                   config.use_full_fidelity,
                   config.cp_mode,
                   config.disable_mapomatic,
                   config.mapomatic_limit,
                   std::nullopt);

        // Capture per-circuit artifacts before we rewrite everything to global space.
        std::vector<IdxType> local_mapping = circuit->get_mapping();
        std::vector<IdxType> local_measurement = build_measurement_mapping(cregs, local_mapping);
        std::vector<QASMTrans::Gate> local_gates = circuit->get_gates();

        SubchipArtifact artifact;
        artifact.index = ci;
        artifact.source_filename = input_files[ci];
        artifact.subchip = subchip;
        artifact.local_to_global = local_to_global;
        artifact.local_gates = local_gates;
        artifact.logical_mapping = local_mapping;
        artifact.measurement_mapping = local_measurement;
        artifact.cregs = cregs;
        artifacts.subchip_artifacts.push_back(std::move(artifact));

        std::size_t classical_bits = 0;
        for (const auto &entry : cregs)
        {
            classical_bits += entry.second.width;
        }
        std::vector<IdxType> global_mapping(local_mapping.size(), -1);
        if (config.debug_level > 0)
        {
            std::cout << "  Logical qubits: " << local_mapping.size()
                      << ", classical bits: " << classical_bits << std::endl;
        }
        // Translate each circuit's local physical qubits into device-global indices.
        for (std::size_t idx = 0; idx < local_mapping.size(); ++idx)
        {
            IdxType local_phys = local_mapping[idx];
            if (local_phys < 0)
            {
                continue;
            }
            if (local_phys >= static_cast<IdxType>(local_to_global.size()))
            {
                global_mapping[idx] = local_phys;
            }
            else
            {
                global_mapping[idx] = local_to_global[static_cast<std::size_t>(local_phys)];
            }
        }

        circuit->set_mapping(global_mapping);
        std::vector<IdxType> global_measurement = build_measurement_mapping(cregs, global_mapping);

        if (config.debug_level > 1)
        {
            std::cout << "  Global mapping:";
            for (IdxType value : global_mapping)
            {
                std::cout << " " << value;
            }
            std::cout << std::endl;
        }

        std::vector<Gate> gates = circuit->get_gates();
        for (auto &gate : gates)
        {
            if (gate.qubit >= 0)
            {
                gate.qubit = local_to_global[static_cast<std::size_t>(gate.qubit)];
            }
            if (gate.ctrl >= 0)
            {
                gate.ctrl = local_to_global[static_cast<std::size_t>(gate.ctrl)];
            }
            if (gate.extra >= 0)
            {
                gate.extra = local_to_global[static_cast<std::size_t>(gate.extra)];
            }
        }
        circuit->set_gates(gates);

        artifacts.combined_gates.insert(artifacts.combined_gates.end(), gates.begin(), gates.end());

        if (global_measurement.empty())
        {
            global_measurement = global_mapping;
        }
        artifacts.combined_mapping.insert(artifacts.combined_mapping.end(), global_measurement.begin(), global_measurement.end());

        // Prefix creg names so multiple circuits can be merged into one QASM.
        std::ostringstream prefix_builder;
        prefix_builder << "circuit" << std::setw(static_cast<int>(std::max<std::size_t>(2, prefix_width))) << std::setfill('0') << ci << "_";
        const std::string prefix = prefix_builder.str();
        if (!artifacts.subchip_artifacts.empty())
        {
            artifacts.subchip_artifacts.back().circuit_prefix = prefix;
        }

        for (const auto &entry : cregs)
        {
            creg renamed = entry.second;
            renamed.name = prefix + entry.first;
            artifacts.combined_cregs.emplace(renamed.name, renamed);
        }
    }
    return artifacts;
}

OutputArtifacts emit_outputs(const CliConfig &config,
                             const std::shared_ptr<Circuit> &combined_circuit,
                             const std::vector<std::string> &input_files,
                             const std::vector<Gate> &combined_gates,
                             IdxType total_requested_qubits,
                             const std::shared_ptr<Chip> &chip)
{
    OutputArtifacts outputs;
    std::string requested_output_path = config.output_path;

    // Emit combined QASM and optional metrics.
    if (config.debug_level > 0)
    {
        std::cout << "======== QASMTrans ========" << std::endl;
        std::cout << "Processed " << input_files.size() << " circuit(s) on backend: " << config.backendpath
                  << " (" << chip->chip_qubit_num << " physical qubits)" << std::endl;
        std::cout << "Combined logical qubits: " << total_requested_qubits << std::endl;
        std::cout << "Basis gate mode: " << config.mode_name << std::endl;
        std::cout << "Limit mode: " << (config.run_with_limit ? "True" : "False") << std::endl;
        if (!config.pulse_template_path.empty())
        {
            std::cout << "Pulse template: " << config.pulse_template_path << std::endl;
        }
        else
        {
            std::cout << "Pulse template: (not provided; skipping pulse dump)" << std::endl;
        }
        if (config.emit_requested)
        {
            std::cout << "QICK emission: enabled (" << (config.emit_run ? "run" : "summary") << " mode)" << std::endl;
        }
    }

    if (config.debug_level > 0)
    {
        const auto summary = compute_gate_summary(combined_gates, combined_circuit->num_qubits());
        std::cout << "[metrics] one_qubit_gates=" << summary.single_qubit
                  << " two_qubit_gates=" << summary.two_qubit
                  << " depth=" << summary.depth << std::endl;
    }

    std::vector<QASMTrans::Gate> expanded_gates;
    bool rigetti_mode = false;
    if (!config.pulse_template_path.empty())
    {
        expanded_gates = QASMTrans::pulses::expandGatesForPulseDump(combined_circuit,
                                                                     config.backendpath,
                                                                     config.pulse_template_path,
                                                                     config.allow_parameterized_merge_candidates,
                                                                     &rigetti_mode);
    }
    if (!expanded_gates.empty() && rigetti_mode)
    {
        outputs.qasm_path = dumpQASMFromGates(combined_circuit,
                                              expanded_gates,
                                              input_files.front().c_str(),
                                              requested_output_path,
                                              config.debug_level,
                                              config.mode);
    }
    else
    {
        outputs.qasm_path = dumpQASM(combined_circuit, input_files.front().c_str(), requested_output_path, config.debug_level, config.mode);
    }
    std::cout << "Saving output qasm to: " << outputs.qasm_path << std::endl;

    // Optional pulse dump.
    if (!config.pulse_template_path.empty())
    {
        outputs.pulses_path = derive_pulse_output_path(outputs.qasm_path);
        dumpPulses(combined_circuit,
                   input_files.front().c_str(),
                   config.backendpath,
                   config.pulse_template_path,
                   outputs.pulses_path,
                   config.debug_level,
                   config.allow_parameterized_merge_candidates);
        std::cout << "Saving output pulses to: " << outputs.pulses_path << std::endl;
    }
    else if (config.debug_level > 0)
    {
        std::cout << "Pulse template not provided; skipping pulse dump." << std::endl;
    }

    // Optional QICK emission against generated pulse schedules.
    if (config.emit_requested)
    {
        if (outputs.pulses_path.empty())
        {
            throw std::runtime_error("Error: unable to emit pulses because the pulse schedule was not generated.");
        }
        QASMTrans::pulses::emit_with_qick(outputs.pulses_path,
                                          config.qick_config_path,
                                          config.emit_run,
                                          !config.emit_run,
                                          config.program_name.c_str(),
                                          config.debug_level > 0);
    }

    return outputs;
}

void export_subchips(const std::vector<SubchipArtifact> &subchip_artifacts,
                     const OutputArtifacts &outputs,
                     const CliConfig &config)
{
    // Ensure we have a spot to stash per-circuit subchip artifacts.
    fs::path final_output_path(outputs.qasm_path);
    fs::path base_output_dir = final_output_path.has_parent_path() ? final_output_path.parent_path() : fs::current_path();
    std::string subchip_dir_name = final_output_path.stem().string() + "_subchips";
    fs::path subchip_dir = base_output_dir / subchip_dir_name;
    std::error_code subchip_ec;
    fs::create_directories(subchip_dir, subchip_ec);
    if (subchip_ec)
    {
        std::cerr << "Warning: failed to create subchip directory '" << subchip_dir << "' (" << subchip_ec.message() << ")" << std::endl;
    }

    // Export each subchip (JSON + standalone QASM) with the original circuit cregs.
    for (const auto &artifact : subchip_artifacts)
    {
        std::string cleaned_prefix = artifact.circuit_prefix;
        if (!cleaned_prefix.empty() && cleaned_prefix.back() == '_')
        {
            cleaned_prefix.pop_back();
        }
        std::string subchip_name = cleaned_prefix.empty() ? "subchip" : cleaned_prefix + "_subchip";
        fs::path json_path = subchip_dir / (subchip_name + ".json");
        fs::path qasm_path = subchip_dir / (subchip_name + ".qasm");

        auto pruned = prune_subchip_artifact(artifact.subchip,
                                             artifact.local_to_global,
                                             artifact.local_gates,
                                             artifact.logical_mapping,
                                             artifact.measurement_mapping,
                                             subchip_name,
                                             config.mode);
        write_json_file(pruned.device_json, json_path);

        auto subchip_circuit = std::make_shared<Circuit>(static_cast<IdxType>(pruned.new_to_old.size()));
        subchip_circuit->set_gates(pruned.gates);
        if (!pruned.measurement_mapping.empty())
        {
            subchip_circuit->set_mapping(pruned.measurement_mapping);
        }
        else
        {
            subchip_circuit->set_mapping(pruned.logical_mapping);
        }
        subchip_circuit->set_creg(artifact.cregs);
        dumpQASM(subchip_circuit, artifact.source_filename.c_str(), qasm_path.string(), config.debug_level > 1 ? config.debug_level : 0, config.mode);
    }
}

int main(int argc, char **argv)
{
    CliConfig config;
    int early_exit = 0;
    if (!parse_cli(argc, argv, config, early_exit))
    {
        return early_exit;
    }

    try
    {
        ingest_backend_metadata(config.backendpath, g_device_basis_gates, g_merged_gate_aliases);

        auto batch = load_circuits(config.input_files);

        IdxType total_requested_qubits = 0;
        for (auto size : batch.circuit_sizes)
        {
            total_requested_qubits += size;
        }
        // Guard against empty/invalid circuits before constructing device model.
        if (total_requested_qubits == 0)
        {
            std::cerr << "Error: no qubits found in provided circuit(s)." << std::endl;
            return 1;
        }

        auto chip = constructChip(total_requested_qubits, config.backendpath,
                                  config.run_with_limit, config.debug_level);
        if (!chip)
        {
            std::cerr << "Error: failed to construct chip from backend." << std::endl;
            return 1;
        }

        if (total_requested_qubits > chip->chip_qubit_num)
        {
            std::cerr << "Error: total logical qubits (" << total_requested_qubits
                      << ") exceed device capacity (" << chip->chip_qubit_num << ")." << std::endl;
            return 1;
        }

        auto overall_start = std::chrono::steady_clock::now();

        auto partition_sizes = compute_partition_sizes(batch.circuit_sizes, chip->chip_qubit_num, total_requested_qubits);
        auto partitions = allocate_partitions(chip, partition_sizes, config.debug_level);
        auto artifacts = transpile_and_merge(batch.circuits,
                                             batch.circuit_cregs,
                                             partitions,
                                             config.input_files,
                                             chip,
                                             config,
                                             total_requested_qubits);

        auto combined_circuit = std::make_shared<Circuit>(chip->chip_qubit_num);
        combined_circuit->set_gates(artifacts.combined_gates);
        combined_circuit->set_creg(artifacts.combined_cregs);
        combined_circuit->set_mapping(artifacts.combined_mapping);

        auto outputs = emit_outputs(config,
                                    combined_circuit,
                                    config.input_files,
                                    artifacts.combined_gates,
                                    total_requested_qubits,
                                    chip);

        export_subchips(artifacts.subchip_artifacts, outputs, config);

        auto overall_end = std::chrono::steady_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(overall_end - overall_start).count();
        std::cout << "[timing] total_ms=" << total_ms << std::endl;
        return 0;
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 1;
    }
}
