#include <memory>
#include <string>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>

#include "../include/QASMTransPrimitives.hpp"
#include "../include/IR/chip.hpp"
#include "../include/parser/parser_util.hpp"
#include "../include/parser/qasm_parser.hpp"
#include "../include/circuit_passes/transpiler.hpp"
#include "../include/util/chip_partition.hpp"

using namespace QASMTrans;

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
    std::cout << "-full_fidelity    Use full-circuit fidelity heuristic instead of critical-path mode" << std::endl;
    std::cout << "-cp_mode <product|additive>  Choose critical-path aggregation strategy (default additive)" << std::endl;
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
    bool use_full_fidelity = false;
    CriticalPathHeuristicMode cp_mode = CriticalPathHeuristicMode::AdditiveAverage;
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
        if (cmdOptionExists(argv, argv + argc, "-full_fidelity"))
        {
            use_full_fidelity = true;
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
            else if (lowered == "additive" || lowered == "sum" || lowered == "average" || lowered == "avg")
            {
                cp_mode = CriticalPathHeuristicMode::AdditiveAverage;
            }
            else
            {
                std::cerr << "Error: unknown -cp_mode value '" << cp_mode_value
                          << "'. Expected 'product' or 'additive'." << std::endl;
                return 1;
            }
        }
        if (cmdOptionExists(argv, argv + argc, "-v"))
        {
            debug_level = IdxType(std::stoi(getCmdOption(argv, argv + argc, "-v")));
        }
        if (cmdOptionExists(argv, argv + argc, "-o"))
        {
            output_path = std::string(getCmdOption(argv, argv + argc, "-o"));
        }
        for (int argi = 1; argi < argc; ++argi)
        {
            std::string current = argv[argi];
            if (current == "-i")
            {
                int next = argi + 1;
                if (next >= argc || argv[next][0] == '-')
                {
                    std::cerr << "Error: -i requires a following QASM file path." << std::endl;
                    return 1;
                }
                while (next < argc && argv[next][0] != '-')
                {
                    input_files.emplace_back(argv[next]);
                    ++next;
                }
                argi = next - 1;
            }
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
        try
        {
            partitions = partition_chip(chip, partition_sizes);
        }
        catch (const std::exception &ex)
        {
            std::cerr << "Error during device partitioning: " << ex.what() << std::endl;
            return 1;
        }

        std::vector<Gate> combined_gates;
        combined_gates.reserve(1024);
        std::vector<IdxType> combined_mapping;
        combined_mapping.reserve(total_requested_qubits);
        std::map<std::string, creg> combined_cregs;

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

            transpiler(circuit, subchip, cregs, debug_level, mode, use_full_fidelity, cp_mode);

            std::vector<IdxType> local_mapping = circuit->get_mapping();
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

            std::vector<IdxType> measurement_mapping;
            for (const auto &entry : cregs)
            {
                const auto &qubit_indices = entry.second.qubit_indices;
                for (std::size_t pos = 0; pos < qubit_indices.size(); ++pos)
                {
                    IdxType logical_index = qubit_indices[pos];
                    IdxType physical = 0;
                    if (logical_index != UN_DEF && logical_index >= 0 &&
                        logical_index < static_cast<IdxType>(global_mapping.size()))
                    {
                        physical = global_mapping[static_cast<std::size_t>(logical_index)];
                    }
                    else if (!global_mapping.empty())
                    {
                        physical = global_mapping[pos % global_mapping.size()];
                    }
                    measurement_mapping.push_back(physical);
                }
            }
            if (measurement_mapping.empty())
            {
                measurement_mapping = global_mapping;
            }
            combined_mapping.insert(combined_mapping.end(), measurement_mapping.begin(), measurement_mapping.end());

            std::ostringstream prefix_builder;
            prefix_builder << "circuit" << std::setw(2) << std::setfill('0') << ci << "_";
            const std::string prefix = prefix_builder.str();

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

        if (debug_level > 0)
        {
            cout << "======== QASMTrans ========" << endl;
            cout << "Processed " << circuits.size() << " circuit(s) on backend: " << backendpath
                 << " (" << chip->chip_qubit_num << " physical qubits)" << endl;
            cout << "Combined logical qubits: " << total_requested_qubits << endl;
            cout << "Basis gate mode: " << mode_name << endl;
            cout << "Limit mode: " << (run_with_limit ? "True" : "False") << endl;
        }

        dumpQASM(combined_circuit, input_files.front().c_str(), output_path, debug_level, mode);
        cout << "Saving output qasm to: " << output_path << endl;
        return 0;
    }
    std::cout << "Invalid Commend Line, Please Check" << std::endl;
    print_help();
    return 0;
}
