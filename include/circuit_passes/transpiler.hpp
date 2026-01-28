#pragma once

#include <algorithm>
#include <bitset>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "../QASMTransPrimitives.hpp"

#include "../IR/gate.hpp"
#include "../IR/circuit.hpp"

#include "../cli_support.hpp"
#include "../dump_qasm.hpp"

#include "routing_mapping.hpp"
#include "fast_quality_routing.hpp"
#include "decompose.hpp"
#include "optimize_1q.hpp"
#include "optimize_2q.hpp"
#include "optimize_2q_synth.hpp"
#include "mapomatic.hpp"
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
    case OP::RZX:
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

inline void enforce_cx_direction(std::shared_ptr<Circuit> circuit,
                                 const std::shared_ptr<Chip> &chip,
                                 const std::unordered_set<std::string> &basis_gates)
{
    if (!circuit || !chip || chip->directed_edge_list.empty())
    {
        return;
    }
    if (basis_gates.empty())
    {
        return;
    }
    const bool has_h = basis_gates.find("h") != basis_gates.end();
    const bool has_rz = basis_gates.find("rz") != basis_gates.end();
    const bool has_sx = basis_gates.find("sx") != basis_gates.end();
    if (!has_h && !(has_rz && has_sx))
    {
        return;
    }
    std::vector<Gate> gates = circuit->get_gates();
    if (gates.empty())
    {
        return;
    }

    auto emit_h = [&](std::vector<Gate> &out, IdxType qubit)
    {
        if (has_h)
        {
            out.emplace_back(OP::H, qubit);
            return;
        }
        out.emplace_back(OP::RZ, qubit, -1, -1, 1, PI / 2);
        out.emplace_back(OP::SX, qubit);
        out.emplace_back(OP::RZ, qubit, -1, -1, 1, PI / 2);
    };

    std::vector<Gate> out;
    out.reserve(gates.size());
    for (const auto &gate : gates)
    {
        if (gate.op_name == OP::CX && !gate.has_custom_name())
        {
            const IdxType ctrl = gate.ctrl;
            const IdxType tgt = gate.qubit;
            bool allowed = false;
            if (ctrl >= 0 && ctrl < static_cast<IdxType>(chip->directed_edge_list.size()))
            {
                const auto &targets = chip->directed_edge_list[static_cast<std::size_t>(ctrl)];
                allowed = targets.find(tgt) != targets.end();
            }
            if (!allowed && ctrl >= 0 && tgt >= 0)
            {
                emit_h(out, ctrl);
                emit_h(out, tgt);
                Gate flipped = gate;
                flipped.ctrl = tgt;
                flipped.qubit = ctrl;
                out.push_back(flipped);
                emit_h(out, ctrl);
                emit_h(out, tgt);
                continue;
            }
        }
        out.push_back(gate);
    }
    circuit->set_gates(out);
}

enum class RoutingMode
{
    Sabre,
    FastQuality
};

void transpiler(shared_ptr<Circuit> circuit,
                shared_ptr<Chip> chip,
                map<string, creg> list_cregs,
                IdxType debug_level,
                IdxType mode,
                bool use_full_fidelity,
                CriticalPathHeuristicMode cp_mode,
                bool disable_mapomatic,
                std::size_t mapomatic_max_embeddings,
                bool enable_1q_opt,
                bool enable_2q_cancel,
                bool enable_commute_2q,
                bool enable_2q_synth,
                bool routing_decay,
                double routing_decay_increment,
                IdxType routing_decay_reset,
                IdxType routing_trials,
                bool enable_sabre_layout,
                RoutingMode routing_mode,
                std::size_t fast_quality_max_embeddings,
                const std::unordered_set<std::string> &basis_gates,
                std::optional<uint64_t> routing_seed = std::nullopt)
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
    auto format_ms = [](double ms) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << ms;
        return oss.str();
    };

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
        cout << "STEP-1. Initial gate decomposition time: " << format_ms(initial_decompose_time) << "ms" << endl;

    //======================================== STEP-2: Routing and Mapping ============================================
    cpu_timer routing_timer;
    routing_timer.start_timer();
    (void)enable_sabre_layout;
    if (routing_mode == RoutingMode::FastQuality)
    {
        fast_quality_routing(circuit, chip, debug_level, fast_quality_max_embeddings,
                             routing_decay, routing_decay_increment, routing_decay_reset,
                             routing_trials, routing_seed, false);
    }
    else
    {
        Routing(circuit, chip, debug_level, routing_decay, routing_decay_increment,
                routing_decay_reset, routing_trials, routing_seed);
    }
    {
        IdxType swap_count = 0;
        IdxType cx_count = 0;
        IdxType twoq_count = 0;
        for (const auto &gate : circuit->get_gates())
        {
            if (gate.op_name == OP::SWAP)
            {
                swap_count += 1;
            }
            if (gate.ctrl >= 0 && gate.qubit >= 0 && gate.extra < 0)
            {
                twoq_count += 1;
                if (gate.op_name == OP::CX)
                {
                    cx_count += 1;
                }
            }
        }
        circuit->set_routing_swap_count(swap_count);
        if (debug_level > 0)
        {
            std::cout << "STEP-2. Routing stats: twoq=" << twoq_count
                      << " cx=" << cx_count
                      << " swap=" << swap_count << std::endl;
        }
    }
    routing_timer.stop_timer();
    double routing_time = routing_timer.measure();
    if (debug_level > 0)
        cout << "STEP-2. Routing and mapping time: " << format_ms(routing_time) << "ms" << endl;
    if (debug_level > 1)
        cout << circuit->to_string() << endl;
    //======================================== STEP-3: Calibration-Aware Optimization =======================================
    double calib_time = 0.0;
    if (!disable_mapomatic)
    {
        cpu_timer calib_timer;
        calib_timer.start_timer();
        calibration_aware_optimization(circuit, chip, debug_level, use_full_fidelity, cp_mode, mapomatic_max_embeddings);
        calib_timer.stop_timer();
        calib_time = calib_timer.measure();
        if (debug_level > 0)
            cout << "STEP-3. Calibration-aware optimization time: " << format_ms(calib_time) << "ms" << endl;
        if (debug_level > 1)
            cout << circuit->to_string() << endl;
    }
    else if (debug_level > 0)
    {
        cout << "STEP-3. Calibration-aware optimization skipped (--disable_mapomatic)" << endl;
    }
    //======================================== STEP-4: Basis Gate Decomposition =======================================
    cpu_timer decompose_timer;
    decompose_timer.start_timer();
    Decompose(circuit, mode);
    bool can_2q_synth = enable_2q_synth;
#ifndef QASMTRANS_USE_EIGEN
    if (can_2q_synth)
    {
        if (debug_level > 0)
        {
            cout << "STEP-4. 2Q synthesis requested, but Eigen is not enabled at build time." << endl;
        }
        can_2q_synth = false;
    }
#endif

    const bool run_opt_loop = can_2q_synth || enable_commute_2q || enable_2q_cancel || enable_1q_opt;
    if (run_opt_loop)
    {
        const int kMaxOptLoops = 50;
        const int kBacktrackDepth = 5;
        QASMTrans::cli::GateSummary summary =
            QASMTrans::cli::compute_gate_summary(circuit->get_gates(), circuit->num_qubits());
        std::size_t best_depth = summary.depth;
        std::size_t best_size =
            static_cast<std::size_t>(summary.single_qubit + summary.two_qubit);
        std::vector<Gate> best_gates = circuit->get_gates();
        int since_best = 0;

        for (int iter = 0; iter < kMaxOptLoops; ++iter)
        {
            if (can_2q_synth)
            {
                QASMTrans::optimize::synthesize_2q_blocks(circuit, chip, &basis_gates);
            }
            if (enable_1q_opt)
            {
                QASMTrans::optimize::optimize_1q_gates_decomposition(circuit, chip, &basis_gates);
            }
            if (enable_commute_2q)
            {
                QASMTrans::optimize::commute_rz_through_2q(circuit);
                QASMTrans::optimize::commute_1q_through_2q(circuit);
                QASMTrans::optimize::commutative_cancel_2q(circuit);
            }
            if (enable_2q_cancel)
            {
                QASMTrans::optimize::cancel_adjacent_2q(circuit);
            }

            QASMTrans::cli::GateSummary next_summary =
                QASMTrans::cli::compute_gate_summary(circuit->get_gates(), circuit->num_qubits());
            std::size_t next_size =
                static_cast<std::size_t>(next_summary.single_qubit + next_summary.two_qubit);

            if (next_summary.depth < best_depth ||
                (next_summary.depth == best_depth && next_size < best_size))
            {
                best_depth = next_summary.depth;
                best_size = next_size;
                best_gates = circuit->get_gates();
                since_best = 0;
            }
            else if (next_summary.depth == best_depth && next_size == best_size)
            {
                break;
            }
            else
            {
                since_best += 1;
                if (since_best >= kBacktrackDepth)
                {
                    circuit->set_gates(best_gates);
                    break;
                }
            }
        }
    }

    enforce_cx_direction(circuit, chip, basis_gates);
    if (enable_1q_opt)
    {
        QASMTrans::optimize::optimize_1q_gates_decomposition(circuit, chip, &basis_gates);
    }
    decompose_timer.stop_timer();
    double decompose_time = decompose_timer.measure();
    if (debug_level > 0)
    {
        cout << "STEP-4. Basis gate decomposition time: " << format_ms(decompose_time) << "ms" << endl;
        cout << " total QASMTrans time: " << format_ms(initial_decompose_time + routing_time + decompose_time) << "ms" << endl;
    }

}
