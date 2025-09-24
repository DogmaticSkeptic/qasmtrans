#pragma once

#include <algorithm>
#include <cctype>
#include <cstring>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../QASMTransPrimitives.hpp"

#include "../IR/gate.hpp"
#include "../IR/circuit.hpp"
#include "../IR/chip.hpp"
#include "../IR/graph.hpp"

#include "../nlomann/json.hpp"

#include <lemon/list_graph.h>
#include <lemon/vf2pp.h>

using namespace QASMTrans;
using namespace std;
using json = nlohmann::json;


namespace mapomatic_detail
{
    inline std::string to_lower_copy(const std::string &value)
    {
        std::string result;
        result.reserve(value.size());
        for (char ch : value)
        {
            result.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
        }
        return result;
    }

    inline IdxType remap_if_present(IdxType qubit, const std::unordered_map<IdxType, IdxType> &mapping)
    {
        auto it = mapping.find(qubit);
        return it != mapping.end() ? it->second : qubit;
    }

    inline double calibration_heuristic(const std::vector<Gate> &gates,
                                        const std::unordered_map<IdxType, IdxType> &physical_mapping,
                                        const shared_ptr<Chip> &chip)
    {
        if (!chip)
        {
            return 1.0;
        }

        double fidelity = 1.0;
        for (const auto &gate : gates)
        {
            if (strcmp(OP_NAMES[gate.op_name], "MA") == 0)
            {
                continue;
            }

            const std::string gate_name = to_lower_copy(OP_NAMES[gate.op_name]);
            double gate_error = 0.0;

            if (gate.ctrl >= 0 && gate.qubit >= 0)
            {
                const IdxType mapped_ctrl = remap_if_present(gate.ctrl, physical_mapping);
                const IdxType mapped_tgt = remap_if_present(gate.qubit, physical_mapping);

                auto pair_it = chip->two_qubit_errors.find({mapped_ctrl, mapped_tgt});
                if (pair_it == chip->two_qubit_errors.end())
                {
                    pair_it = chip->two_qubit_errors.find({mapped_tgt, mapped_ctrl});
                }

                if (pair_it != chip->two_qubit_errors.end())
                {
                    auto err_it = pair_it->second.find(gate_name);
                    if (err_it != pair_it->second.end())
                    {
                        gate_error = std::clamp(err_it->second, 0.0, 1.0);
                    }
                }
            }
            else if (gate.qubit >= 0)
            {
                const IdxType mapped_qubit = remap_if_present(gate.qubit, physical_mapping);
                if (mapped_qubit >= 0 && mapped_qubit < static_cast<IdxType>(chip->single_qubit_errors.size()))
                {
                    const auto &error_map = chip->single_qubit_errors[mapped_qubit];
                    auto err_it = error_map.find(gate_name);
                    if (err_it != error_map.end())
                    {
                        gate_error = std::clamp(err_it->second, 0.0, 1.0);
                    }
                }
            }

            fidelity *= (1.0 - gate_error);
        }

        return 1.0 - fidelity;
    }
}


void calibration_aware_optimization(shared_ptr<Circuit> circuit, 
                  shared_ptr<Chip> chip, 
                  IdxType debug_level)
{
    (void)debug_level;
    if (!circuit || !chip)
    {
        return;
    }

    std::vector<Gate> gates = circuit->get_gates();
    if (gates.empty())
    {
        return;
    }

    std::unordered_set<IdxType> used_qubits_set;
    used_qubits_set.reserve(gates.size() * 2);
    for (const auto &gate : gates)
    {
        if (gate.qubit >= 0)
        {
            used_qubits_set.insert(gate.qubit);
        }
        if (gate.ctrl >= 0)
        {
            used_qubits_set.insert(gate.ctrl);
        }
        if (gate.extra >= 0)
        {
            used_qubits_set.insert(gate.extra);
        }
    }

    if (used_qubits_set.empty())
    {
        return;
    }

    std::vector<IdxType> used_qubits(used_qubits_set.begin(), used_qubits_set.end());
    std::sort(used_qubits.begin(), used_qubits.end());

    lemon::ListGraph circuit_graph;
    lemon::ListGraph::NodeMap<IdxType> circuit_node_to_qubit(circuit_graph);
    std::unordered_map<IdxType, lemon::ListGraph::Node> qubit_to_circuit_node;
    qubit_to_circuit_node.reserve(used_qubits.size());

    for (IdxType qubit : used_qubits)
    {
        auto node = circuit_graph.addNode();
        circuit_node_to_qubit[node] = qubit;
        qubit_to_circuit_node.emplace(qubit, node);
    }

    std::set<std::pair<IdxType, IdxType>> circuit_edges;
    auto add_circuit_edge = [&](IdxType a, IdxType b)
    {
        if (a < 0 || b < 0 || a == b)
        {
            return;
        }
        auto it_a = qubit_to_circuit_node.find(a);
        auto it_b = qubit_to_circuit_node.find(b);
        if (it_a == qubit_to_circuit_node.end() || it_b == qubit_to_circuit_node.end())
        {
            return;
        }
        auto ordered = std::minmax(a, b);
        if (circuit_edges.insert(ordered).second)
        {
            circuit_graph.addEdge(it_a->second, it_b->second);
        }
    };

    for (const auto &gate : gates)
    {
        if (gate.ctrl >= 0 && gate.qubit >= 0)
        {
            add_circuit_edge(gate.ctrl, gate.qubit);
        }
        if (gate.ctrl >= 0 && gate.extra >= 0)
        {
            add_circuit_edge(gate.ctrl, gate.extra);
        }
        if (gate.qubit >= 0 && gate.extra >= 0)
        {
            add_circuit_edge(gate.qubit, gate.extra);
        }
    }

    lemon::ListGraph chip_graph;
    lemon::ListGraph::NodeMap<IdxType> chip_node_to_qubit(chip_graph);
    std::vector<lemon::ListGraph::Node> chip_nodes;
    chip_nodes.reserve(chip->chip_qubit_num);

    for (IdxType phys = 0; phys < chip->chip_qubit_num; ++phys)
    {
        auto node = chip_graph.addNode();
        chip_node_to_qubit[node] = phys;
        chip_nodes.push_back(node);
    }

    for (IdxType src = 0; src < static_cast<IdxType>(chip->edge_list.size()); ++src)
    {
        if (src >= static_cast<IdxType>(chip_nodes.size()))
        {
            break;
        }
        for (IdxType dst : chip->edge_list[src])
        {
            if (dst >= 0 && dst < static_cast<IdxType>(chip_nodes.size()) && src < dst)
            {
                chip_graph.addEdge(chip_nodes[src], chip_nodes[dst]);
            }
        }
    }

    if (lemon::countNodes(circuit_graph) == 0 || lemon::countNodes(chip_graph) == 0)
    {
        return;
    }

    lemon::ListGraph::NodeMap<lemon::ListGraph::Node> vf2_mapping(circuit_graph);
    lemon::ListGraph::NodeMap<int> circuit_labels(circuit_graph, 0);
    lemon::ListGraph::NodeMap<int> chip_labels(chip_graph, 0);

    constexpr std::size_t kMaxEmbeddings = 128;
    std::vector<std::unordered_map<IdxType, IdxType>> candidate_mappings;
    candidate_mappings.reserve(kMaxEmbeddings);

    lemon::Vf2pp<lemon::ListGraph, lemon::ListGraph,
                 lemon::ListGraph::NodeMap<lemon::ListGraph::Node>,
                 lemon::ListGraph::NodeMap<int>,
                 lemon::ListGraph::NodeMap<int>>
        vf2_algorithm(circuit_graph, chip_graph, vf2_mapping, circuit_labels, chip_labels);
    vf2_algorithm.mappingType(lemon::SUBGRAPH);

    std::size_t explored = 0;
    while (explored < kMaxEmbeddings && vf2_algorithm.find())
    {
        std::unordered_map<IdxType, IdxType> assignment;
        assignment.reserve(used_qubits.size());
        for (lemon::ListGraph::NodeIt node(circuit_graph); node != lemon::INVALID; ++node)
        {
            const IdxType original_qubit = circuit_node_to_qubit[node];
            const auto mapped_node = vf2_mapping[node];
            const IdxType mapped_qubit = chip_node_to_qubit[mapped_node];
            assignment.emplace(original_qubit, mapped_qubit);
        }
        candidate_mappings.push_back(std::move(assignment));
        ++explored;
    }

    if (candidate_mappings.empty())
    {
        std::unordered_map<IdxType, IdxType> identity;
        identity.reserve(used_qubits.size());
        for (IdxType qubit : used_qubits)
        {
            identity.emplace(qubit, qubit);
        }
        candidate_mappings.push_back(std::move(identity));
    }

    double best_score = std::numeric_limits<double>::max();
    std::size_t best_index = 0;
    for (std::size_t idx = 0; idx < candidate_mappings.size(); ++idx)
    {
        double score = mapomatic_detail::calibration_heuristic(gates, candidate_mappings[idx], chip);
        if (score < best_score)
        {
            best_score = score;
            best_index = idx;
        }
    }

    const auto &best_mapping = candidate_mappings[best_index];

    bool mapping_changed = false;
    std::vector<Gate> updated_gates = gates;
    for (auto &gate : updated_gates)
    {
        if (gate.qubit >= 0)
        {
            IdxType mapped = mapomatic_detail::remap_if_present(gate.qubit, best_mapping);
            mapping_changed = mapping_changed || mapped != gate.qubit;
            gate.qubit = mapped;
        }
        if (gate.ctrl >= 0)
        {
            IdxType mapped = mapomatic_detail::remap_if_present(gate.ctrl, best_mapping);
            mapping_changed = mapping_changed || mapped != gate.ctrl;
            gate.ctrl = mapped;
        }
        if (gate.extra >= 0)
        {
            IdxType mapped = mapomatic_detail::remap_if_present(gate.extra, best_mapping);
            mapping_changed = mapping_changed || mapped != gate.extra;
            gate.extra = mapped;
        }
    }

    if (mapping_changed)
    {
        circuit->set_gates(updated_gates);

        std::vector<IdxType> logical_to_physical = circuit->get_mapping();
        if (!logical_to_physical.empty())
        {
            for (auto &physical : logical_to_physical)
            {
                physical = mapomatic_detail::remap_if_present(physical, best_mapping);
            }
            circuit->set_mapping(logical_to_physical);
        }

        circuit->populate_connectivity();
    }
}
