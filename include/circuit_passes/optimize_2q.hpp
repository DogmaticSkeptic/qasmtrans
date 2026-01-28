#pragma once

#include <cmath>
#include <vector>

#include "../IR/circuit.hpp"
#include "../IR/gate.hpp"
#include "../QASMTransPrimitives.hpp"

namespace QASMTrans
{
    namespace optimize
    {
        namespace detail_2q
        {
            constexpr double kAngleEps = 1e-12;

            inline bool is_two_qubit_gate(const Gate &gate)
            {
                return gate.ctrl >= 0 && gate.qubit >= 0 && gate.extra < 0;
            }

            inline bool is_single_qubit_gate(const Gate &gate)
            {
                return gate.ctrl < 0 && gate.extra < 0 && gate.qubit >= 0;
            }

            inline bool is_single_qubit_rz(const Gate &gate)
            {
                return gate.ctrl < 0 && gate.extra < 0 && gate.qubit >= 0 && gate.op_name == OP::RZ;
            }

            inline bool is_symmetric_gate(OP op)
            {
                switch (op)
                {
                case OP::CZ:
                case OP::SWAP:
                case OP::ISWAP:
                case OP::ZZ:
                case OP::RXX:
                case OP::RYY:
                case OP::RZZ:
                    return true;
                default:
                    return false;
                }
            }

            inline bool pair_matches(const Gate &a, const Gate &b)
            {
                if (a.ctrl == b.ctrl && a.qubit == b.qubit)
                {
                    return true;
                }
                if (is_symmetric_gate(a.op_name))
                {
                    return a.ctrl == b.qubit && a.qubit == b.ctrl;
                }
                return false;
            }

            inline bool is_param_2q(OP op)
            {
                switch (op)
                {
                case OP::RXX:
                case OP::RYY:
                case OP::RZZ:
                case OP::RZX:
                case OP::ZZ:
                case OP::CRX:
                case OP::CRY:
                case OP::CRZ:
                case OP::CP:
                    return true;
                default:
                    return false;
                }
            }

            inline bool is_self_inverse_2q(OP op)
            {
                switch (op)
                {
                case OP::CX:
                case OP::CZ:
                case OP::SWAP:
                case OP::ECR:
                    return true;
                default:
                    return false;
                }
            }

            inline bool gate_overlaps(const Gate &gate, IdxType q0, IdxType q1)
            {
                if (gate.qubit == q0 || gate.qubit == q1)
                {
                    return true;
                }
                if (gate.ctrl == q0 || gate.ctrl == q1)
                {
                    return true;
                }
                if (gate.extra == q0 || gate.extra == q1)
                {
                    return true;
                }
                return false;
            }

            inline bool is_nonunitary_gate(const Gate &gate)
            {
                return gate.op_name == OP::M || gate.op_name == OP::MA || gate.op_name == OP::RESET;
            }

            inline bool is_z_rotation_1q(OP op)
            {
                switch (op)
                {
                case OP::Z:
                case OP::RZ:
                case OP::S:
                case OP::SDG:
                case OP::T:
                case OP::TDG:
                case OP::P:
                case OP::RI:
                    return true;
                default:
                    return false;
                }
            }

            inline bool is_x_rotation_1q(OP op)
            {
                switch (op)
                {
                case OP::X:
                case OP::RX:
                case OP::SX:
                    return true;
                default:
                    return false;
                }
            }

            inline bool commutes_with_twoq_1q(const Gate &twoq, const Gate &oneq)
            {
                if (!is_two_qubit_gate(twoq) || !is_single_qubit_gate(oneq))
                {
                    return false;
                }
                if (twoq.op_name == OP::CX)
                {
                    if (oneq.qubit == twoq.ctrl)
                    {
                        return is_z_rotation_1q(oneq.op_name);
                    }
                    if (oneq.qubit == twoq.qubit)
                    {
                        return is_x_rotation_1q(oneq.op_name);
                    }
                    return false;
                }
                if (twoq.op_name == OP::CZ || twoq.op_name == OP::ZZ || twoq.op_name == OP::RZZ)
                {
                    if (oneq.qubit == twoq.ctrl || oneq.qubit == twoq.qubit)
                    {
                        return is_z_rotation_1q(oneq.op_name);
                    }
                    return false;
                }
                if (twoq.op_name == OP::RXX)
                {
                    if (oneq.qubit == twoq.ctrl || oneq.qubit == twoq.qubit)
                    {
                        return is_x_rotation_1q(oneq.op_name);
                    }
                    return false;
                }
                return false;
            }

            inline bool near_zero(double value)
            {
                return std::abs(value) < kAngleEps;
            }

            inline bool rz_commutes_with_twoq(const Gate &rz, const Gate &twoq)
            {
                if (!is_single_qubit_rz(rz))
                {
                    return false;
                }
                if (!is_two_qubit_gate(twoq))
                {
                    return false;
                }
                switch (twoq.op_name)
                {
                case OP::CX:
                    return rz.qubit == twoq.ctrl;
                case OP::CZ:
                case OP::ZZ:
                case OP::RZZ:
                case OP::CP:
                    return (rz.qubit == twoq.ctrl) || (rz.qubit == twoq.qubit);
                default:
                    return false;
                }
            }
        } // namespace detail_2q

        inline void commute_rz_through_2q(std::shared_ptr<Circuit> circuit)
        {
            if (!circuit)
            {
                return;
            }
            std::vector<Gate> gates = circuit->get_gates();
            if (gates.size() < 2)
            {
                return;
            }

            std::size_t i = 0;
            while (i + 1 < gates.size())
            {
                Gate &cur = gates[i];
                Gate &next = gates[i + 1];
                if (!cur.has_custom_name() && detail_2q::rz_commutes_with_twoq(cur, next))
                {
                    std::swap(cur, next);
                    if (i + 1 < gates.size())
                    {
                        ++i;
                        continue;
                    }
                }
                ++i;
            }
            circuit->set_gates(gates);
        }

        inline void commute_1q_through_2q(std::shared_ptr<Circuit> circuit)
        {
            if (!circuit)
            {
                return;
            }
            std::vector<Gate> gates = circuit->get_gates();
            if (gates.size() < 2)
            {
                return;
            }

            constexpr int kMaxPasses = 1;
            for (int pass = 0; pass < kMaxPasses; ++pass)
            {
                bool changed = false;
                for (std::size_t i = 0; i + 1 < gates.size(); ++i)
                {
                    Gate &g0 = gates[i];
                    Gate &g1 = gates[i + 1];
                    if (g0.has_custom_name() || g1.has_custom_name())
                    {
                        continue;
                    }
                    if (detail_2q::is_nonunitary_gate(g0) || detail_2q::is_nonunitary_gate(g1))
                    {
                        continue;
                    }
                    // Move 1q gates rightwards across commuting 2q gates to merge 1q runs.
                    if (detail_2q::is_single_qubit_gate(g0) &&
                        detail_2q::is_two_qubit_gate(g1) &&
                        detail_2q::commutes_with_twoq_1q(g1, g0))
                    {
                        std::swap(g0, g1);
                        changed = true;
                    }
                }
                if (!changed)
                {
                    break;
                }
            }
            circuit->set_gates(gates);
        }

        inline void cancel_adjacent_2q(std::shared_ptr<Circuit> circuit)
        {
            if (!circuit)
            {
                return;
            }
            std::vector<Gate> gates = circuit->get_gates();
            if (gates.size() < 2)
            {
                return;
            }

            std::vector<Gate> out;
            out.reserve(gates.size());

            std::size_t idx = 0;
            while (idx < gates.size())
            {
                const Gate &gate = gates[idx];
                if (idx + 1 < gates.size())
                {
                    const Gate &next = gates[idx + 1];
                    if (!gate.has_custom_name() && !next.has_custom_name() &&
                        gate.op_name == next.op_name &&
                        detail_2q::is_two_qubit_gate(gate) &&
                        detail_2q::is_two_qubit_gate(next) &&
                        detail_2q::pair_matches(gate, next))
                    {
                        if (detail_2q::is_param_2q(gate.op_name))
                        {
                            double combined = gate.theta + next.theta;
                            if (!detail_2q::near_zero(combined))
                            {
                                Gate merged = gate;
                                merged.theta = combined;
                                out.push_back(merged);
                            }
                            idx += 2;
                            continue;
                        }
                        if (detail_2q::is_self_inverse_2q(gate.op_name))
                        {
                            idx += 2;
                            continue;
                        }
                    }
                }
                out.push_back(gate);
                ++idx;
            }

            circuit->set_gates(out);
        }

        inline void commutative_cancel_2q(std::shared_ptr<Circuit> circuit)
        {
            if (!circuit)
            {
                return;
            }
            std::vector<Gate> gates = circuit->get_gates();
            if (gates.size() < 2)
            {
                return;
            }

            std::vector<char> remove(gates.size(), 0);
            const std::size_t n = gates.size();
            for (std::size_t i = 0; i < n; ++i)
            {
                if (remove[i])
                {
                    continue;
                }
                const Gate &gate = gates[i];
                if (gate.has_custom_name())
                {
                    continue;
                }
                if (!detail_2q::is_two_qubit_gate(gate) || !detail_2q::is_self_inverse_2q(gate.op_name))
                {
                    continue;
                }
                const IdxType q0 = gate.ctrl;
                const IdxType q1 = gate.qubit;
                for (std::size_t j = i + 1; j < n; ++j)
                {
                    if (remove[j])
                    {
                        continue;
                    }
                    const Gate &cand = gates[j];
                    if (cand.has_custom_name())
                    {
                        if (detail_2q::gate_overlaps(cand, q0, q1))
                        {
                            break;
                        }
                        continue;
                    }
                    if (!detail_2q::gate_overlaps(cand, q0, q1))
                    {
                        continue;
                    }
                    if (detail_2q::is_nonunitary_gate(cand))
                    {
                        break;
                    }
                    if (detail_2q::is_two_qubit_gate(cand))
                    {
                        if (cand.op_name == gate.op_name && detail_2q::pair_matches(gate, cand))
                        {
                            remove[i] = 1;
                            remove[j] = 1;
                        }
                        break;
                    }
                    if (!detail_2q::commutes_with_twoq_1q(gate, cand))
                    {
                        break;
                    }
                }
            }

            std::vector<Gate> out;
            out.reserve(gates.size());
            for (std::size_t i = 0; i < gates.size(); ++i)
            {
                if (!remove[i])
                {
                    out.push_back(gates[i]);
                }
            }
            circuit->set_gates(out);
        }
    } // namespace optimize
} // namespace QASMTrans
