#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <bitset>
#include <iomanip>
#include <cctype>

#include "../QASMTransPrimitives.hpp"

#include "../IR/gate.hpp"
#include "../IR/circuit.hpp"

#include "../dump_qasm.hpp"

#include "routing_mapping.hpp"
#include "decompose.hpp"
#include "remapping.hpp"

using namespace QASMTrans;
using namespace std;

inline std::string logical_label_for_gate(const Gate &gate)
{
    std::string name = OP_NAMES[gate.op_name];
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c)
                   { return static_cast<char>(std::tolower(c)); });
    std::ostringstream oss;
    oss << name;
    auto append_param = [&](const char *key, ValType value)
    {
        if (value != 0.0)
        {
            oss << "[" << key << "=" << std::setprecision(10) << value << "]";
        }
    };
    switch (gate.op_name)
    {
    case OP::RX:
    case OP::RY:
    case OP::RZ:
    case OP::RI:
    case OP::P:
    case OP::CRX:
    case OP::CRY:
    case OP::CRZ:
    case OP::CP:
    case OP::U:
    case OP::CU:
    case OP::RXX:
    case OP::RYY:
    case OP::RZZ:
        append_param("theta", gate.theta);
        append_param("phi", gate.phi);
        append_param("lambda", gate.lam);
        append_param("gamma", gate.gamma);
        break;
    default:
        break;
    }
    return oss.str();
}

inline void annotate_logical_gates(std::vector<Gate> &gates)
{
    IdxType next_id = 0;
    for (auto &gate : gates)
    {
        gate.set_logical_metadata(next_id++, logical_label_for_gate(gate));
    }
}

void transpiler(shared_ptr<Circuit> circuit, shared_ptr<Chip> chip, map<string, creg> list_cregs, IdxType debug_level, IdxType mode)
{
    circuit->set_creg(list_cregs);
    IdxType n_qubits = IdxType(circuit->num_qubits());
    IdxType chip_n_qubit = chip->chip_qubit_num;

    if (n_qubits > chip_n_qubit)
    {
        //std::cerr<<"Chip qubit number is smaller than the circuit."<<endl;
        //std::cerr<<"No transpilation has been performed."<<endl;
        throw std::logic_error{"Chip qubit number is smaller than the circuit. No transpilation has been performed."};
        std::exit(1);
    }

    //======================================== STEP-1: Initial Gate Decomposition =====================================
    cpu_timer initial_decompose_timer;
    initial_decompose_timer.start_timer();
    Decompose_three_to_two(circuit);
    initial_decompose_timer.stop_timer();
    double initial_decompose_time = initial_decompose_timer.measure();
    {
        std::vector<Gate> logical_stage_gates = circuit->get_gates();
        annotate_logical_gates(logical_stage_gates);
        circuit->set_gates(logical_stage_gates);
    }
    if (debug_level > 0)
        cout << "STEP-1. Initial gate decomposition time: " << (IdxType)initial_decompose_time << "ms" << endl;

    //======================================== STEP-2: Routing and Mapping ============================================
    cpu_timer routing_timer;
    routing_timer.start_timer();
    Routing(circuit, chip, debug_level);
    routing_timer.stop_timer();
    double routing_time = routing_timer.measure();
    if (debug_level > 0)
        cout << "STEP-2. Routing and mapping time: " << (IdxType)routing_time << "ms" << endl;
    if (debug_level > 1)
        cout << circuit->to_string() << endl;
    //======================================== STEP-3: Basis Gate Decomposition =======================================
    cpu_timer decompose_timer;
    decompose_timer.start_timer();
    Decompose(circuit, mode);
    decompose_timer.stop_timer();
    double decompose_time = decompose_timer.measure();
    if (debug_level > 0)
    {
        cout << "STEP-3. Basis gate decomposition time: " << (IdxType)decompose_time << "ms" << endl;
        cout << " total QASMTrans time: " << (IdxType)(initial_decompose_time + routing_time + decompose_time) << "ms" << endl;
    }
}
