#pragma once

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <utility>   // for std::make_pair
#include <algorithm> // for toLowerCase
#include <filesystem>
#include <system_error>
#include <cstring>

#include "QASMTransPrimitives.hpp"
#include "IR/gate.hpp"
#include "IR/circuit.hpp"

using namespace QASMTrans;

std::string toLowerCase(const std::string &str)
{
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });
    return result;
}
namespace
{
inline std::string mode_prefix(IdxType mode)
{
    switch (mode)
    {
    case 0:
        return "transpiled_IBMQ_";
    case 1:
        return "transpiled_IonQ_";
    case 2:
        return "transpiled_Quantinuum_";
    case 3:
        return "transpiled_Rigetti_";
    case 4:
        return "transpiled_Quafu_";
    default:
        return "transpiled_";
    }
}
} // namespace

// Function to write QASM file
std::string dumpQASM(std::shared_ptr<QASMTrans::Circuit> circuit, const char *filename, const std::string &requested_output_path, IdxType debug_level, IdxType mode)
{
    namespace fs = std::filesystem;
    std::string input_name = filename ? fs::path(filename).filename().string() : "circuit.qasm";

    const bool treat_as_directory = requested_output_path.empty() ||
                                    requested_output_path.back() == '/' ||
                                    requested_output_path.back() == '\\' ||
                                    fs::is_directory(fs::path(requested_output_path));

    fs::path final_path;
    if (treat_as_directory)
    {
        fs::path directory = requested_output_path.empty() ? fs::path("../data/output_qasm_file") : fs::path(requested_output_path);
        if (directory.filename() == ".")
        {
            directory = directory.parent_path();
        }
        std::error_code ec;
        fs::create_directories(directory, ec);
        if (ec)
        {
            std::cerr << "Warning: failed to create directory '" << directory << "' (" << ec.message() << ")" << std::endl;
        }
        final_path = directory / (mode_prefix(mode) + input_name);
    }
    else
    {
        final_path = fs::path(requested_output_path);
        fs::path parent = final_path.parent_path();
        if (!parent.empty())
        {
            std::error_code ec;
            fs::create_directories(parent, ec);
            if (ec)
            {
                std::cerr << "Warning: failed to create directory '" << parent << "' (" << ec.message() << ")" << std::endl;
            }
        }
    }

    IdxType n_qubits = IdxType(circuit->num_qubits());
    std::vector<QASMTrans::Gate> gate_info = circuit->get_gates();
    map<string, creg> cregs = circuit->get_cregs();
    std::map<std::string, IdxType> basis_gate_counts;

    std::ofstream qasm_file(final_path);
    if (!qasm_file.is_open())
    {
        std::cerr << "Error: unable to open QASM output file '" << final_path << "'" << std::endl;
        return final_path.string();
    }

    qasm_file << "OPENQASM 2.0;\n";
    qasm_file << "include \"qelib1.inc\";\n";
    qasm_file << "qreg q[" << n_qubits << "];\n";
    for (auto &creg : cregs)
    {
        qasm_file << "creg " << toLowerCase(creg.first) << "[" << creg.second.width << "];\n";
    }
    for (auto g : gate_info)
    {
        if (std::strcmp(QASMTrans::OP_NAMES[g.op_name], "MA") != 0)
        {
            if (!g.gateToString().empty())
            {
                qasm_file << toLowerCase(g.gateToString()) << "; \n";
                std::string gate_name = toLowerCase(g.name());
                basis_gate_counts[gate_name] += 1;
            }
        }
    }
    IdxType creg_index = 0;
    for (auto &creg : cregs)
    {
        for (auto &qubit : creg.second.qubit_indices)
        {
            qasm_file << "measure q[" << circuit->initial_mapping[creg_index] << "] -> " << toLowerCase(creg.first) << "[" << creg_index << "];\n";
            ++creg_index;
        }
    }
    qasm_file.close();

    if (debug_level > 0)
    {
        IdxType n_gates = circuit->num_gates();
        std::cout << "Transpiled circuit contains " << circuit->num_qubits() << " qubits and " << n_gates << " basis gates: ";
        for (auto &gate : basis_gate_counts)
        {
            std::cout << gate.first << ":" << gate.second << ", ";
        }
        std::cout << std::endl;
    }

    return final_path.string();
}
