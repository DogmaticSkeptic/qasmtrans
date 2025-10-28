#include <memory>
#include <string>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <system_error>
#include <chrono>
#include <optional>
#include <array>

#include "../include/QASMTransPrimitives.hpp"
#include "../include/IR/chip.hpp"
#include "../include/parser/parser_util.hpp"
#include "../include/parser/qasm_parser.hpp"
#include "../include/circuit_passes/transpiler.hpp"
#include "../include/util/chip_partition.hpp"
#include "../include/nlomann/json.hpp"

using namespace QASMTrans;
namespace fs = std::filesystem;

namespace
{

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
        mapping = logical_to_physical;
    }
    return mapping;
}

std::vector<std::string> basis_gates_for_mode(IdxType mode)
{
    switch (mode)
    {
    case 0:
        return {"rz", "sx", "x", "cx"};
    case 1:
        return {"rx", "ry", "rz", "rxx"};
    case 2:
        return {"rx", "rz", "zz"};
    case 3:
        return {"rx", "ry", "cz"};
    case 4:
        return {"cz", "rx", "ry", "rz", "h"};
    default:
        return {};
    }
}

struct GateSummary
{
    std::size_t single_qubit = 0;
    std::size_t two_qubit = 0;
    std::size_t depth = 0;
};

GateSummary compute_gate_summary(const std::vector<Gate> &gates, IdxType initial_capacity)
{
    GateSummary summary;
    if (initial_capacity < 0)
    {
        initial_capacity = 0;
    }
    std::vector<std::size_t> qubit_depth(static_cast<std::size_t>(initial_capacity), 0);

    for (const auto &gate : gates)
    {
        switch (gate.op_name)
        {
        case OP::M:
        case OP::MA:
        case OP::RESET:
            continue;
        default:
            break;
        }

        std::array<IdxType, 3> raw_indices = {gate.qubit, gate.ctrl, gate.extra};
        std::vector<std::size_t> involved;
        involved.reserve(raw_indices.size());

        for (IdxType raw_index : raw_indices)
        {
            if (raw_index < 0)
            {
                continue;
            }
            std::size_t index = static_cast<std::size_t>(raw_index);
            if (index >= qubit_depth.size())
            {
                qubit_depth.resize(index + 1, 0);
            }
            if (std::find(involved.begin(), involved.end(), index) == involved.end())
            {
                involved.push_back(index);
            }
        }

        if (involved.empty())
        {
            continue;
        }

        if (involved.size() == 1)
        {
            summary.single_qubit += 1;
        }
        else if (involved.size() == 2)
        {
            summary.two_qubit += 1;
        }

        std::size_t gate_depth = 0;
        for (std::size_t idx : involved)
        {
            gate_depth = std::max(gate_depth, qubit_depth[idx]);
        }
        gate_depth += 1;
        for (std::size_t idx : involved)
        {
            qubit_depth[idx] = gate_depth;
        }
        summary.depth = std::max(summary.depth, gate_depth);
    }

    return summary;
}

nlohmann::json build_subchip_json(const std::shared_ptr<Chip> &subchip,
                                  const std::vector<IdxType> &new_to_old,
                                  const std::vector<IdxType> &old_to_new,
                                  const std::vector<IdxType> &local_to_global,
                                  const std::vector<IdxType> &logical_mapping,
                                  const std::string &name,
                                  IdxType mode)
{
    using nlohmann::json;
    json result;
    result["name"] = name;
    result["version"] = "generated";
    result["num_qubits"] = static_cast<IdxType>(new_to_old.size());
    result["basis_gates"] = basis_gates_for_mode(mode);

    std::vector<std::string> coupling;
    if (subchip)
    {
        for (std::size_t new_src = 0; new_src < new_to_old.size(); ++new_src)
        {
            IdxType old_src = new_to_old[new_src];
            if (old_src < 0 || old_src >= static_cast<IdxType>(subchip->edge_list.size()))
            {
                continue;
            }
            for (IdxType old_dst : subchip->edge_list[old_src])
            {
                if (old_dst < 0 || old_dst >= static_cast<IdxType>(old_to_new.size()))
                {
                    continue;
                }
                IdxType new_dst = old_to_new[old_dst];
                if (new_dst >= 0)
                {
                    coupling.emplace_back(std::to_string(new_src) + "_" + std::to_string(new_dst));
                }
            }
        }
    }
    std::sort(coupling.begin(), coupling.end());
    coupling.erase(std::unique(coupling.begin(), coupling.end()), coupling.end());
    result["cx_coupling"] = coupling;

    json gate_errs = json::object();
    json gate_lens = json::object();
    if (subchip)
    {
        for (std::size_t new_idx = 0; new_idx < new_to_old.size(); ++new_idx)
        {
            IdxType old_idx = new_to_old[new_idx];
            if (old_idx >= 0 && old_idx < static_cast<IdxType>(subchip->single_qubit_errors.size()))
            {
                for (const auto &entry : subchip->single_qubit_errors[old_idx])
                {
                    gate_errs[entry.first + std::to_string(new_idx)] = entry.second;
                }
            }
            if (old_idx >= 0 && old_idx < static_cast<IdxType>(subchip->single_qubit_gate_lengths.size()))
            {
                for (const auto &entry : subchip->single_qubit_gate_lengths[old_idx])
                {
                    gate_lens[entry.first + std::to_string(new_idx)] = entry.second;
                }
            }
        }
        auto append_two_qubit = [&](const auto &source_map, json &dest) {
            for (const auto &entry : source_map)
            {
                IdxType new_ctrl = (entry.first.first >= 0 && entry.first.first < static_cast<IdxType>(old_to_new.size()))
                                       ? old_to_new[entry.first.first]
                                       : -1;
                IdxType new_tgt = (entry.first.second >= 0 && entry.first.second < static_cast<IdxType>(old_to_new.size()))
                                      ? old_to_new[entry.first.second]
                                      : -1;
                if (new_ctrl < 0 || new_tgt < 0)
                {
                    continue;
                }
                for (const auto &gate_entry : entry.second)
                {
                    const double value = gate_entry.second;
                    const std::string forward_key = gate_entry.first + std::to_string(new_ctrl) + "_" + std::to_string(new_tgt);
                    dest[forward_key] = value;
                    if (new_ctrl != new_tgt)
                    {
                        const std::string reverse_key = gate_entry.first + std::to_string(new_tgt) + "_" + std::to_string(new_ctrl);
                        if (!dest.contains(reverse_key))
                        {
                            dest[reverse_key] = value;
                        }
                    }
                }
            }
        };
        append_two_qubit(subchip->two_qubit_errors, gate_errs);
        append_two_qubit(subchip->two_qubit_gate_lengths, gate_lens);
    }
    result["gate_errs"] = gate_errs;
    result["gate_lens"] = gate_lens;

    auto append_optional_property = [&](const std::vector<std::optional<double>> &values,
                                        const char *key) {
        if (values.empty())
        {
            return;
        }
        json prop = json::object();
        bool any = false;
        for (std::size_t new_idx = 0; new_idx < new_to_old.size(); ++new_idx)
        {
            IdxType old_idx = new_to_old[new_idx];
            if (old_idx < 0 || old_idx >= static_cast<IdxType>(values.size()))
            {
                continue;
            }
            const auto &val = values[static_cast<std::size_t>(old_idx)];
            if (val.has_value())
            {
                prop[std::to_string(new_idx)] = *val;
                any = true;
            }
        }
        if (any)
        {
            result[key] = std::move(prop);
        }
    };
    if (subchip)
    {
        append_optional_property(subchip->t1, "T1");
        append_optional_property(subchip->t2, "T2");
        append_optional_property(subchip->freq, "freq");
        append_optional_property(subchip->readout_length, "readout_length");
        append_optional_property(subchip->prob_meas0_prep1, "prob_meas0_prep1");
        append_optional_property(subchip->prob_meas1_prep0, "prob_meas1_prep0");
    }

    nlohmann::json local_global = nlohmann::json::array();
    for (IdxType value : local_to_global)
    {
        local_global.push_back(value);
    }
    result["local_to_global"] = local_global;

    nlohmann::json logical_local = nlohmann::json::array();
    for (IdxType value : logical_mapping)
    {
        logical_local.push_back(value);
    }
    result["logical_to_local"] = logical_local;
    result["generated_by"] = "qasmtrans";

    return result;
}

void write_json_file(const nlohmann::json &data, const fs::path &output_path)
{
    if (output_path.has_parent_path())
    {
        std::error_code ec;
        fs::create_directories(output_path.parent_path(), ec);
        if (ec)
        {
            std::cerr << "Warning: failed to create directory '" << output_path.parent_path() << "' (" << ec.message() << ")" << std::endl;
        }
    }
    std::ofstream out(output_path);
    if (!out.is_open())
    {
        std::cerr << "Error: unable to open JSON output file '" << output_path << "'" << std::endl;
        return;
    }
    out << data.dump(2);
    out.close();
}

struct SubchipArtifact
{
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

struct PrunedSubchipData
{
    std::vector<IdxType> new_to_old;
    std::vector<IdxType> old_to_new;
    std::vector<QASMTrans::Gate> gates;
    std::vector<IdxType> logical_mapping;
    std::vector<IdxType> measurement_mapping;
    std::vector<IdxType> local_to_global;
    nlohmann::json device_json;
};

std::vector<IdxType> collect_used_nodes(const SubchipArtifact &artifact)
{
    std::vector<IdxType> used;
    if (!artifact.subchip)
    {
        return used;
    }

    IdxType local_size = artifact.subchip->chip_qubit_num;
    if (local_size <= 0)
    {
        return used;
    }

    std::vector<char> mark(static_cast<std::size_t>(local_size), 0);
    auto mark_idx = [&](IdxType idx)
    {
        if (idx >= 0 && idx < local_size)
        {
            mark[static_cast<std::size_t>(idx)] = 1;
        }
    };

    for (const auto &gate : artifact.local_gates)
    {
        mark_idx(gate.qubit);
        mark_idx(gate.ctrl);
        mark_idx(gate.extra);
    }
    for (IdxType value : artifact.logical_mapping)
    {
        mark_idx(value);
    }
    for (IdxType value : artifact.measurement_mapping)
    {
        mark_idx(value);
    }

    for (IdxType idx = 0; idx < local_size; ++idx)
    {
        if (mark[static_cast<std::size_t>(idx)])
        {
            used.push_back(idx);
        }
    }
    return used;
}

PrunedSubchipData prune_subchip_artifact(const SubchipArtifact &artifact,
                                         const std::string &subchip_name,
                                         IdxType mode)
{
    PrunedSubchipData result;

    if (!artifact.subchip)
    {
        result.device_json = build_subchip_json(nullptr, result.new_to_old, result.old_to_new,
                                                result.local_to_global, result.logical_mapping,
                                                subchip_name, mode);
        return result;
    }

    auto used_nodes = collect_used_nodes(artifact);
    if (used_nodes.empty())
    {
        used_nodes.reserve(static_cast<std::size_t>(artifact.subchip->chip_qubit_num));
        for (IdxType idx = 0; idx < artifact.subchip->chip_qubit_num; ++idx)
        {
            used_nodes.push_back(idx);
        }
    }

    std::sort(used_nodes.begin(), used_nodes.end());
    used_nodes.erase(std::unique(used_nodes.begin(), used_nodes.end()), used_nodes.end());

    result.new_to_old = used_nodes;
    result.old_to_new.assign(static_cast<std::size_t>(artifact.subchip->chip_qubit_num), -1);
    result.local_to_global.reserve(used_nodes.size());
    for (std::size_t new_idx = 0; new_idx < used_nodes.size(); ++new_idx)
    {
        IdxType old_idx = used_nodes[new_idx];
        if (old_idx >= 0 && old_idx < static_cast<IdxType>(result.old_to_new.size()))
        {
            result.old_to_new[static_cast<std::size_t>(old_idx)] = static_cast<IdxType>(new_idx);
        }
        if (old_idx >= 0 && old_idx < static_cast<IdxType>(artifact.local_to_global.size()))
        {
            result.local_to_global.push_back(artifact.local_to_global[static_cast<std::size_t>(old_idx)]);
        }
        else
        {
            result.local_to_global.push_back(old_idx);
        }
    }

    result.gates.reserve(artifact.local_gates.size());
    for (const auto &gate : artifact.local_gates)
    {
        Gate adjusted = gate;
        if (adjusted.qubit >= 0 && adjusted.qubit < static_cast<IdxType>(result.old_to_new.size()))
        {
            adjusted.qubit = result.old_to_new[static_cast<std::size_t>(adjusted.qubit)];
        }
        if (adjusted.ctrl >= 0 && adjusted.ctrl < static_cast<IdxType>(result.old_to_new.size()))
        {
            adjusted.ctrl = result.old_to_new[static_cast<std::size_t>(adjusted.ctrl)];
        }
        if (adjusted.extra >= 0 && adjusted.extra < static_cast<IdxType>(result.old_to_new.size()))
        {
            adjusted.extra = result.old_to_new[static_cast<std::size_t>(adjusted.extra)];
        }
        result.gates.push_back(adjusted);
    }

    result.logical_mapping = artifact.logical_mapping;
    for (auto &value : result.logical_mapping)
    {
        if (value >= 0 && value < static_cast<IdxType>(result.old_to_new.size()))
        {
            value = result.old_to_new[static_cast<std::size_t>(value)];
        }
        else
        {
            value = -1;
        }
    }

    result.measurement_mapping = artifact.measurement_mapping;
    for (auto &value : result.measurement_mapping)
    {
        if (value >= 0 && value < static_cast<IdxType>(result.old_to_new.size()))
        {
            value = result.old_to_new[static_cast<std::size_t>(value)];
        }
        else
        {
            value = -1;
        }
    }

    result.device_json = build_subchip_json(artifact.subchip,
                                            result.new_to_old,
                                            result.old_to_new,
                                            result.local_to_global,
                                            result.logical_mapping,
                                            subchip_name,
                                            mode);
    return result;
}

} // namespace

void print_help()
{
    // print the help function for all the options
    std::cout << "Usage: ./qasmtrans [options]" << std::endl;
    std::cout << "Option            Description" << std::endl;
    std::cout << "-i <path>         Input QASM file (repeat -i for multiple circuits)" << std::endl;
    std::cout << "-c <backend>      Path to backend configuration json file" << std::endl;
    std::cout << "-limited          Run the transpiler with limited physical qubits usage" << std::endl;
    std::cout << "-limited          Limit qubit usage to circuit than device. "
        << "It reduces qubit usage but may introduce extra routing cost or unable to route" << std::endl;
    std::cout << "-backend_list     Print the available device backends" << std::endl;
    std::cout << "-m <name>         Set the transpiler targeted device, default is ibmq" << std::endl;
    std::cout << "-v <0/1/2>        Set the output level, default is 0" << std::endl;
    std::cout << "-full_fidelity    Score Mapomatic candidates on the entire circuit instead of its critical path" << std::endl;
    std::cout << "-cp_mode <product|hybrid>  Choose scoring strategy (default product)" << std::endl;
    std::cout << "-mapomatic_limit <N>       Limit the number of candidate embeddings Mapomatic evaluates (default 1000)" << std::endl;
    std::cout << "--disable_mapomatic       Skip the calibration-aware Mapomatic pass" << std::endl;
    std::cout << "-o <path>         Set the output file, "
        << "default is data/output/transpiled_modename_filename.qasm" << std::endl;
    std::cout << "-h                print the help function" << std::endl;
}

int main(int argc, char **argv)
{
    bool run_with_limit = false;
    IdxType mode = 0;
    std::string mode_name = "ibmq";
    IdxType debug_level = 0;
    std::string output_path = "../data/output/";
    std::string backendpath;
    std::vector<std::string> input_files;
    bool disable_mapomatic = false;
    bool use_full_fidelity = false;
    CriticalPathHeuristicMode cp_mode = CriticalPathHeuristicMode::LogProduct;
    std::size_t mapomatic_limit = 1000;
    std::map<std::string, IdxType> machineQubits = {
        {"ibmq_toronto", 27},
        {"ibmq_jakarta", 7},
        {"ibmq_guadalupe", 16},
        {"ibm_seattle", 433},
        {"ibm_cairo", 27},
        {"ibm_brisbane", 127},
        {"dummy_ibmq12", 12},
        {"dummy_ibmq14", 14},
        {"dummy_ibmq15", 15},
        {"dummy_ibmq16", 16},
        {"dummy_ibmq30", 30},
        {"aspen_m3", 80},
        {"h1_2", 12},
        {"h1_1", 20}};
    if (argc == 1)
    {
        print_help();
        return 0;
    }
    else
    {
        if (cmdOptionExists(argv, argv + argc, "-h"))
        {
            print_help();
            return 0;
        }
        //! need a -m for different machine mode basis
        if (cmdOptionExists(argv, argv + argc, "-limited"))
        {
            run_with_limit = true;
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
                cp_mode = CriticalPathHeuristicMode::LogProduct;
            }
            else if (lowered == "hybrid")
            {
                cp_mode = CriticalPathHeuristicMode::Hybrid;
            }
            else
            {
                std::cerr << "Error: unknown -cp_mode value '" << cp_mode_value
                          << "'. Expected 'product' or 'hybrid'." << std::endl;
                return 1;
            }
        }
        if (cmdOptionExists(argv, argv + argc, "-mapomatic_limit"))
        {
            const char *opt = getCmdOption(argv, argv + argc, "-mapomatic_limit");
            if (!opt)
            {
                std::cerr << "Error: -mapomatic_limit requires a positive integer argument." << std::endl;
                return 1;
            }
            try
            {
                long long parsed = std::stoll(opt);
                if (parsed <= 0)
                {
                    std::cerr << "Error: -mapomatic_limit must be greater than zero (got " << parsed << ")." << std::endl;
                    return 1;
                }
                mapomatic_limit = static_cast<std::size_t>(parsed);
            }
            catch (const std::exception &)
            {
                std::cerr << "Error: failed to parse -mapomatic_limit argument '" << opt << "'." << std::endl;
                return 1;
            }
        }
        if (cmdOptionExists(argv, argv + argc, "-v"))
        {
            debug_level = IdxType(std::stoi(getCmdOption(argv, argv + argc, "-v")));
        }
        if (cmdOptionExists(argv, argv + argc, "-full_fidelity"))
        {
            use_full_fidelity = true;
        }
        if (cmdOptionExists(argv, argv + argc, "-o"))
        {
            output_path = std::string(getCmdOption(argv, argv + argc, "-o"));
        }
        int argi = 1;
        while (argi < argc)
        {
            std::string current = argv[argi];
            if (current == "-i")
            {
                if (argi + 1 >= argc)
                {
                    std::cerr << "Error: -i requires a following QASM file path." << std::endl;
                    return 1;
                }
                input_files.emplace_back(argv[argi + 1]);
                argi += 2;
                continue;
            }
            if (current == "--disable_mapomatic")
            {
                disable_mapomatic = true;
                ++argi;
                continue;
            }
            ++argi;
        }
        if (cmdOptionExists(argv, argv + argc, "-c"))
        {
            backendpath = std::string(getCmdOption(argv, argv + argc, "-c"));
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
            return 0;
        }
        if (input_files.empty())
        {
            std::cerr << "Error: missing input QASM file(s) (-i)." << std::endl;
            return 1;
        }
        if (backendpath.empty())
        {
            std::cerr << "Error: missing machine backend file via -c" << std::endl;
            return 1;
        }

        if (cmdOptionExists(argv, argv + argc, "-m"))
        {
            mode_name = std::string(getCmdOption(argv, argv + argc, "-m"));
            if (mode_name == "ibmq" || mode_name == "IBMQ")
            {
                mode = 0;
            }
            else if (mode_name == "ionq" || mode_name == "IonQ")
            {
                mode = 1;
            }
            else if (mode_name == "Quantinuum" || mode_name == "quantinuum")
            {
                mode = 2;
            }
            else if (mode_name == "Rigetti" || mode_name == "rigetti")
            {
                mode = 3;
            }
            else if (mode_name == "Quafu" || mode_name == "quafu")
            {
                mode = 4;
            }
            else
            {
                std::cout << "Invalid mode name, please check" << std::endl;
                return 0;
            }
        }

        std::vector<std::shared_ptr<Circuit>> circuits;
        std::vector<map<std::string, creg>> circuit_cregs;
        std::vector<IdxType> circuit_sizes;
        circuits.reserve(input_files.size());
        circuit_cregs.reserve(input_files.size());
        circuit_sizes.reserve(input_files.size());

        for (const auto &file : input_files)
        {
            qasm_parser parser(file.c_str());
            IdxType n_qubits = parser.num_qubits();
            auto circuit = make_shared<Circuit>(n_qubits);
            parser.loadin_circuit(circuit);
            map<string, creg> cregs = parser.get_list_cregs();
            circuits.push_back(circuit);
            circuit_cregs.push_back(cregs);
            circuit_sizes.push_back(n_qubits);
        }

        IdxType total_requested_qubits = 0;
        for (auto size : circuit_sizes)
        {
            total_requested_qubits += size;
        }
        if (total_requested_qubits == 0)
        {
            std::cerr << "Error: no qubits found in provided circuit(s)." << std::endl;
            return 1;
        }

        shared_ptr<Chip> chip = constructChip(total_requested_qubits, backendpath,
                                              run_with_limit, debug_level);
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

        std::vector<IdxType> partition_sizes = circuit_sizes;
        if (!partition_sizes.empty())
        {
            const double total_logical = static_cast<double>(total_requested_qubits);
            const IdxType device_capacity = chip->chip_qubit_num;
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
        }

        std::vector<std::vector<IdxType>> partitions;
        auto partition_start = std::chrono::steady_clock::now();
        bool partition_success = false;
        try
        {
            partitions = partition_chip(chip, partition_sizes);
            partition_success = true;
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
                std::cerr << "Error: contiguous fallback partitioning failed." << std::endl;
                return 1;
            }
        }
        auto partition_end = std::chrono::steady_clock::now();
        auto partition_ms = std::chrono::duration_cast<std::chrono::milliseconds>(partition_end - partition_start).count();
        std::cout << "[timing] partition_ms=" << partition_ms << std::endl;

        std::vector<Gate> combined_gates;
        combined_gates.reserve(1024);
        std::vector<IdxType> combined_mapping;
        combined_mapping.reserve(total_requested_qubits);
        std::map<std::string, creg> combined_cregs;
        std::vector<SubchipArtifact> subchip_artifacts;
        subchip_artifacts.reserve(circuits.size());

        try
        {
            for (std::size_t ci = 0; ci < circuits.size(); ++ci)
            {
                auto &circuit = circuits[ci];
                const auto &cregs = circuit_cregs[ci];

                std::vector<IdxType> local_to_global;
                auto subchip = make_subchip(chip, partitions[ci], local_to_global);

            if (debug_level > 0)
            {
                std::cout << "Circuit " << ci << ": allocated " << partitions[ci].size()
                          << " physical qubits." << std::endl;
            }
            if (debug_level > 1)
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

            transpiler(circuit, subchip, cregs, debug_level, mode, use_full_fidelity, cp_mode, disable_mapomatic, mapomatic_limit);

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
            subchip_artifacts.push_back(std::move(artifact));

            std::size_t classical_bits = 0;
            for (const auto &entry : cregs)
            {
                classical_bits += entry.second.width;
            }
            std::vector<IdxType> global_mapping(local_mapping.size(), -1);
            if (debug_level > 0)
            {
                std::cout << "  Logical qubits: " << local_mapping.size()
                          << ", classical bits: " << classical_bits << std::endl;
            }
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

            if (debug_level > 1)
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

            combined_gates.insert(combined_gates.end(), gates.begin(), gates.end());

            if (global_measurement.empty())
            {
                global_measurement = global_mapping;
            }
            combined_mapping.insert(combined_mapping.end(), global_measurement.begin(), global_measurement.end());

            std::ostringstream prefix_builder;
            prefix_builder << "circuit" << std::setw(2) << std::setfill('0') << ci << "_";
            const std::string prefix = prefix_builder.str();
            if (!subchip_artifacts.empty())
            {
                subchip_artifacts.back().circuit_prefix = prefix;
            }

            for (const auto &entry : cregs)
            {
                creg renamed = entry.second;
                renamed.name = prefix + entry.first;
                    combined_cregs.emplace(renamed.name, renamed);
                }
            }
        }
        catch (const std::exception &ex)
        {
            std::cerr << "Error while processing circuits: " << ex.what() << std::endl;
            return 1;
        }

        auto combined_circuit = std::make_shared<Circuit>(chip->chip_qubit_num);
        combined_circuit->set_gates(combined_gates);
        combined_circuit->set_creg(combined_cregs);
        combined_circuit->set_mapping(combined_mapping);
        std::string requested_output_path = output_path;

        if (debug_level > 0)
        {
            cout << "======== QASMTrans ========" << endl;
            cout << "Processed " << circuits.size() << " circuit(s) on backend: " << backendpath
                 << " (" << chip->chip_qubit_num << " physical qubits)" << endl;
            cout << "Combined logical qubits: " << total_requested_qubits << endl;
            cout << "Basis gate mode: " << mode_name << endl;
            cout << "Limit mode: " << (run_with_limit ? "True" : "False") << endl;
        }

        if (debug_level > 0)
        {
            const auto summary = compute_gate_summary(combined_gates, combined_circuit->num_qubits());
            std::cout << "[metrics] one_qubit_gates=" << summary.single_qubit
                      << " two_qubit_gates=" << summary.two_qubit
                      << " depth=" << summary.depth << std::endl;
        }

        std::string final_output_qasm = dumpQASM(combined_circuit, input_files.front().c_str(), requested_output_path, debug_level, mode);
        output_path = final_output_qasm;
        cout << "Saving output qasm to: " << final_output_qasm << endl;

        fs::path final_output_path(final_output_qasm);
        fs::path base_output_dir = final_output_path.has_parent_path() ? final_output_path.parent_path() : fs::current_path();
        std::string subchip_dir_name = final_output_path.stem().string() + "_subchips";
        fs::path subchip_dir = base_output_dir / subchip_dir_name;
        std::error_code subchip_ec;
        fs::create_directories(subchip_dir, subchip_ec);
        if (subchip_ec)
        {
            std::cerr << "Warning: failed to create subchip directory '" << subchip_dir << "' (" << subchip_ec.message() << ")" << std::endl;
        }

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

            auto pruned = prune_subchip_artifact(artifact, subchip_name, mode);
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
            dumpQASM(subchip_circuit, artifact.source_filename.c_str(), qasm_path.string(), debug_level > 1 ? debug_level : 0, mode);
        }

        auto overall_end = std::chrono::steady_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(overall_end - overall_start).count();
        std::cout << "[timing] total_ms=" << total_ms << std::endl;
        return 0;
    }
    std::cout << "Invalid Commend Line, Please Check" << std::endl;
    print_help();
        return 0;
    }
