#pragma once

#include <algorithm>
#include <chrono>
#include <deque>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../IR/chip.hpp"
#include "../IR/circuit.hpp"
#include "../IR/gate.hpp"
#include "routing_mapping.hpp"

#include <lemon/list_graph.h>
#include <lemon/vf2pp.h>

using namespace QASMTrans;

namespace QASMTrans
{
    inline void fast_quality_routing(std::shared_ptr<Circuit> circuit,
                                     const std::shared_ptr<Chip> &chip,
                                     IdxType debug_level,
                                     std::size_t max_embeddings,
                                     bool enable_decay,
                                     double decay_increment,
                                     IdxType decay_reset,
                                     IdxType routing_trials,
                                     std::optional<uint64_t> routing_seed,
                                     bool enable_sabre_layout)
    {
        if (!circuit || !chip)
        {
            return;
        }

        const IdxType logical_qubits = static_cast<IdxType>(circuit->num_qubits());
        if (logical_qubits <= 0)
        {
            return;
        }

        const IdxType physical_qubits = chip->qubit_num;
        if (physical_qubits < logical_qubits)
        {
            return;
        }

        std::vector<Gate> gates = circuit->get_gates();
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
        std::map<std::pair<IdxType, IdxType>, double> circuit_edge_weights;
        std::map<std::pair<IdxType, IdxType>, double> directed_edge_weights;
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
            circuit_edge_weights[ordered] += 1.0;
        };
        auto add_directed_edge = [&](IdxType src, IdxType dst)
        {
            if (src < 0 || dst < 0 || src == dst)
            {
                return;
            }
            directed_edge_weights[{src, dst}] += 1.0;
        };

        for (const auto &gate : gates)
        {
            if (gate.ctrl >= 0 && gate.qubit >= 0)
            {
                add_circuit_edge(gate.ctrl, gate.qubit);
                add_directed_edge(gate.ctrl, gate.qubit);
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
        if (circuit_edges.empty() && used_qubits.size() > 1)
        {
            for (std::size_t i = 1; i < used_qubits.size(); ++i)
            {
                add_circuit_edge(used_qubits[i - 1], used_qubits[i]);
            }
        }

        std::size_t embedding_limit = std::max<std::size_t>(1, max_embeddings);
        std::unordered_map<IdxType, IdxType> best_mapping;
        double best_cost = std::numeric_limits<double>::infinity();
        bool found_embedding = false;

        auto has_directed_edge = [&](IdxType src, IdxType dst) -> bool
        {
            if (src < 0 || dst < 0)
            {
                return false;
            }
            const auto &directed = chip->two_qubit_errors;
            if (!directed.empty())
            {
                if (directed.find({src, dst}) != directed.end())
                {
                    return true;
                }
                return false;
            }
            for (IdxType neighbor : chip->edge_list[src])
            {
                if (neighbor == dst)
                {
                    return true;
                }
            }
            return false;
        };

        const double direction_penalty = 2.0;

        auto run_vf2 = [&](lemon::ListGraph &chip_graph,
                           lemon::ListGraph::NodeMap<IdxType> &chip_node_to_qubit) -> std::size_t {
            if (lemon::countNodes(circuit_graph) == 0 || lemon::countNodes(chip_graph) == 0)
            {
                return 0;
            }
            lemon::ListGraph::NodeMap<lemon::ListGraph::Node> vf2_mapping(circuit_graph);
            lemon::ListGraph::NodeMap<int> circuit_labels(circuit_graph, 0);
            lemon::ListGraph::NodeMap<int> chip_labels(chip_graph, 0);

            lemon::Vf2pp<lemon::ListGraph, lemon::ListGraph,
                         lemon::ListGraph::NodeMap<lemon::ListGraph::Node>,
                         lemon::ListGraph::NodeMap<int>,
                         lemon::ListGraph::NodeMap<int>>
                vf2_algorithm(circuit_graph, chip_graph, vf2_mapping, circuit_labels, chip_labels);
            vf2_algorithm.mappingType(lemon::SUBGRAPH);

            const auto vf2_start = std::chrono::steady_clock::now();
            if (debug_level > 0)
            {
                std::cout << "Fast routing_mapping: VF2++ search begin (circuit nodes="
                          << used_qubits.size() << ", chip nodes=" << physical_qubits << ")."
                          << std::endl;
            }
            std::size_t explored = 0;
            while (explored < embedding_limit && vf2_algorithm.find())
            {
                std::unordered_map<IdxType, IdxType> assignment;
                assignment.reserve(used_qubits.size());
                for (lemon::ListGraph::NodeIt node(circuit_graph); node != lemon::INVALID; ++node)
                {
                    const IdxType logical = circuit_node_to_qubit[node];
                    const auto mapped_node = vf2_mapping[node];
                    const IdxType physical = chip_node_to_qubit[mapped_node];
                    assignment.emplace(logical, physical);
                }

                double cost = 0.0;
                for (const auto &edge : circuit_edges)
                {
                    auto it_a = assignment.find(edge.first);
                    auto it_b = assignment.find(edge.second);
                    if (it_a == assignment.end() || it_b == assignment.end())
                    {
                        continue;
                    }
                    IdxType pa = it_a->second;
                    IdxType pb = it_b->second;
                    if (pa < 0 || pb < 0 ||
                        pa >= static_cast<IdxType>(chip->distance_mat.size()) ||
                        pb >= static_cast<IdxType>(chip->distance_mat.size()))
                    {
                        cost = std::numeric_limits<double>::infinity();
                        break;
                    }
                    IdxType dist = chip->distance_mat[pa][pb];
                    if (dist == std::numeric_limits<IdxType>::max())
                    {
                        cost = std::numeric_limits<double>::infinity();
                        break;
                    }
                    double weight = 1.0;
                    auto weight_it = circuit_edge_weights.find(edge);
                    if (weight_it != circuit_edge_weights.end())
                    {
                        weight = weight_it->second;
                    }
                    cost += weight * static_cast<double>(dist);
                }
                for (const auto &entry : directed_edge_weights)
                {
                    const auto &edge = entry.first;
                    auto it_src = assignment.find(edge.first);
                    auto it_dst = assignment.find(edge.second);
                    if (it_src == assignment.end() || it_dst == assignment.end())
                    {
                        continue;
                    }
                    IdxType ps = it_src->second;
                    IdxType pd = it_dst->second;
                    if (ps >= 0 && pd >= 0 &&
                        ps < static_cast<IdxType>(chip->distance_mat.size()) &&
                        pd < static_cast<IdxType>(chip->distance_mat.size()))
                    {
                        if (chip->distance_mat[ps][pd] == 1 && !has_directed_edge(ps, pd))
                        {
                            cost += entry.second * direction_penalty;
                        }
                    }
                }

                if (cost < best_cost)
                {
                    best_cost = cost;
                    best_mapping = std::move(assignment);
                    found_embedding = true;
                }
                ++explored;
                if (debug_level > 1 && (explored % 25 == 0))
                {
                    const auto now = std::chrono::steady_clock::now();
                    const auto elapsed_ms =
                        std::chrono::duration_cast<std::chrono::milliseconds>(now - vf2_start).count();
                    std::cout << "  VF2++ embeddings explored=" << explored
                              << ", elapsed=" << elapsed_ms << " ms" << std::endl;
                }
            }
            if (debug_level > 0)
            {
                const auto now = std::chrono::steady_clock::now();
                const auto elapsed_ms =
                    std::chrono::duration_cast<std::chrono::milliseconds>(now - vf2_start).count();
                std::cout << "Fast routing_mapping: VF2++ search end, explored=" << explored
                          << ", found=" << (found_embedding ? "yes" : "no")
                          << ", elapsed=" << elapsed_ms << " ms" << std::endl;
            }
            return explored;
        };

        std::size_t explored = 0;
        {
            lemon::ListGraph chip_graph;
            lemon::ListGraph::NodeMap<IdxType> chip_node_to_qubit(chip_graph);
            std::vector<lemon::ListGraph::Node> chip_nodes;
            chip_nodes.reserve(static_cast<std::size_t>(physical_qubits));
            for (IdxType phys = 0; phys < physical_qubits; ++phys)
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
            explored = run_vf2(chip_graph, chip_node_to_qubit);
        }

        if (!found_embedding)
        {
            constexpr IdxType kMaxAugDistance = 2;
            lemon::ListGraph chip_graph;
            lemon::ListGraph::NodeMap<IdxType> chip_node_to_qubit(chip_graph);
            std::vector<lemon::ListGraph::Node> chip_nodes;
            chip_nodes.reserve(static_cast<std::size_t>(physical_qubits));
            for (IdxType phys = 0; phys < physical_qubits; ++phys)
            {
                auto node = chip_graph.addNode();
                chip_node_to_qubit[node] = phys;
                chip_nodes.push_back(node);
            }
            for (IdxType src = 0; src < physical_qubits; ++src)
            {
                for (IdxType dst = src + 1; dst < physical_qubits; ++dst)
                {
                    if (src >= static_cast<IdxType>(chip->distance_mat.size()) ||
                        dst >= static_cast<IdxType>(chip->distance_mat[src].size()))
                    {
                        continue;
                    }
                    IdxType dist = chip->distance_mat[src][dst];
                    if (dist <= kMaxAugDistance)
                    {
                        chip_graph.addEdge(chip_nodes[src], chip_nodes[dst]);
                    }
                }
            }
            explored = run_vf2(chip_graph, chip_node_to_qubit);
        }

        if (!found_embedding)
        {
            // Pre-layout seed: greedy walk for heaviest logical qubits.
            std::vector<double> logical_weight(static_cast<std::size_t>(logical_qubits), 0.0);
            for (const auto &entry : circuit_edge_weights)
            {
                const auto &edge = entry.first;
                double weight = entry.second;
                if (edge.first >= 0 && edge.first < logical_qubits)
                {
                    logical_weight[static_cast<std::size_t>(edge.first)] += weight;
                }
                if (edge.second >= 0 && edge.second < logical_qubits)
                {
                    logical_weight[static_cast<std::size_t>(edge.second)] += weight;
                }
            }

            std::vector<std::pair<IdxType, double>> logical_rank;
            logical_rank.reserve(static_cast<std::size_t>(logical_qubits));
            for (IdxType q = 0; q < logical_qubits; ++q)
            {
                logical_rank.emplace_back(q, logical_weight[static_cast<std::size_t>(q)]);
            }
            std::sort(logical_rank.begin(), logical_rank.end(),
                      [](const auto &a, const auto &b)
                      { return a.second > b.second; });

            std::vector<double> phys_score(static_cast<std::size_t>(physical_qubits), 0.0);
            for (IdxType p = 0; p < physical_qubits; ++p)
            {
                double degree = 0.0;
                if (p >= 0 && p < static_cast<IdxType>(chip->edge_list.size()))
                {
                    degree = static_cast<double>(chip->edge_list[static_cast<std::size_t>(p)].size());
                }
                double dist_sum = 0.0;
                if (p >= 0 && p < static_cast<IdxType>(chip->distance_mat.size()))
                {
                    for (IdxType d : chip->distance_mat[static_cast<std::size_t>(p)])
                    {
                        if (d < std::numeric_limits<IdxType>::max())
                        {
                            dist_sum += static_cast<double>(d);
                        }
                    }
                }
                phys_score[static_cast<std::size_t>(p)] = degree * 10.0 - dist_sum;
            }

            IdxType start = 0;
            double best = -std::numeric_limits<double>::infinity();
            for (IdxType p = 0; p < physical_qubits; ++p)
            {
                double score = phys_score[static_cast<std::size_t>(p)];
                if (score > best)
                {
                    best = score;
                    start = p;
                }
            }

            std::vector<IdxType> chosen;
            chosen.reserve(static_cast<std::size_t>(logical_qubits));
            std::unordered_set<IdxType> used_phys;
            used_phys.reserve(static_cast<std::size_t>(logical_qubits));
            chosen.push_back(start);
            used_phys.insert(start);
            while (chosen.size() < static_cast<std::size_t>(logical_qubits))
            {
                IdxType best_next = -1;
                double best_score = -std::numeric_limits<double>::infinity();
                for (IdxType cur : chosen)
                {
                    if (cur < 0 || cur >= static_cast<IdxType>(chip->edge_list.size()))
                    {
                        continue;
                    }
                    for (IdxType nb : chip->edge_list[static_cast<std::size_t>(cur)])
                    {
                        if (nb < 0 || nb >= physical_qubits || used_phys.count(nb))
                        {
                            continue;
                        }
                        double score = phys_score[static_cast<std::size_t>(nb)];
                        if (score > best_score)
                        {
                            best_score = score;
                            best_next = nb;
                        }
                    }
                }
                if (best_next < 0)
                {
                    for (IdxType p = 0; p < physical_qubits; ++p)
                    {
                        if (!used_phys.count(p))
                        {
                            best_next = p;
                            break;
                        }
                    }
                }
                if (best_next < 0)
                {
                    break;
                }
                chosen.push_back(best_next);
                used_phys.insert(best_next);
            }

            std::unordered_map<IdxType, IdxType> heuristic;
            heuristic.reserve(static_cast<std::size_t>(logical_qubits));
            const std::size_t limit = std::min(logical_rank.size(), chosen.size());
            for (std::size_t i = 0; i < limit; ++i)
            {
                heuristic.emplace(logical_rank[i].first, chosen[i]);
            }
            best_mapping = std::move(heuristic);
        }

        std::vector<IdxType> mapping(static_cast<std::size_t>(logical_qubits), -1);
        std::unordered_set<IdxType> used_phys;
        used_phys.reserve(best_mapping.size());
        for (const auto &entry : best_mapping)
        {
            if (entry.first >= 0 && entry.first < logical_qubits &&
                entry.second >= 0 && entry.second < physical_qubits)
            {
                mapping[static_cast<std::size_t>(entry.first)] = entry.second;
                used_phys.insert(entry.second);
            }
        }

        for (IdxType q = 0; q < logical_qubits; ++q)
        {
            if (mapping[static_cast<std::size_t>(q)] >= 0)
            {
                continue;
            }
            for (IdxType phys = 0; phys < physical_qubits; ++phys)
            {
                if (!used_phys.count(phys))
                {
                    mapping[static_cast<std::size_t>(q)] = phys;
                    used_phys.insert(phys);
                    break;
                }
            }
            if (mapping[static_cast<std::size_t>(q)] < 0)
            {
                mapping[static_cast<std::size_t>(q)] = q;
            }
        }

        auto compute_twoq_depth = [&](const std::vector<IdxType> &logical_to_phys) -> std::size_t
        {
            std::vector<std::size_t> phys_depth(static_cast<std::size_t>(physical_qubits), 0);
            std::size_t depth = 0;
            for (const auto &gate : gates)
            {
                if (gate.ctrl < 0 || gate.qubit < 0 || gate.extra >= 0)
                {
                    continue;
                }
                if (gate.ctrl >= logical_qubits || gate.qubit >= logical_qubits)
                {
                    continue;
                }
                IdxType p0 = logical_to_phys[static_cast<std::size_t>(gate.ctrl)];
                IdxType p1 = logical_to_phys[static_cast<std::size_t>(gate.qubit)];
                if (p0 < 0 || p1 < 0)
                {
                    continue;
                }
                std::size_t next = std::max(phys_depth[static_cast<std::size_t>(p0)],
                                            phys_depth[static_cast<std::size_t>(p1)]) +
                                   1;
                phys_depth[static_cast<std::size_t>(p0)] = next;
                phys_depth[static_cast<std::size_t>(p1)] = next;
                depth = std::max(depth, next);
            }
            return depth;
        };

        if (!circuit_edges.empty())
        {
            std::vector<std::vector<std::pair<IdxType, double>>> adjacency(
                static_cast<std::size_t>(logical_qubits));
            adjacency.reserve(static_cast<std::size_t>(logical_qubits));
            for (const auto &edge : circuit_edges)
            {
                double weight = 1.0;
                auto weight_it = circuit_edge_weights.find(edge);
                if (weight_it != circuit_edge_weights.end())
                {
                    weight = weight_it->second;
                }
                if (edge.first >= 0 && edge.first < logical_qubits &&
                    edge.second >= 0 && edge.second < logical_qubits)
                {
                    adjacency[static_cast<std::size_t>(edge.first)].push_back({edge.second, weight});
                    adjacency[static_cast<std::size_t>(edge.second)].push_back({edge.first, weight});
                }
            }

            auto distance_lookup = [&](IdxType a, IdxType b) -> double
            {
                if (a < 0 || b < 0 ||
                    a >= static_cast<IdxType>(chip->distance_mat.size()) ||
                    b >= static_cast<IdxType>(chip->distance_mat.size()))
                {
                    return std::numeric_limits<double>::infinity();
                }
                IdxType dist = chip->distance_mat[a][b];
                if (dist == std::numeric_limits<IdxType>::max())
                {
                    return std::numeric_limits<double>::infinity();
                }
                return static_cast<double>(dist);
            };

            auto total_cost = [&]() -> double
            {
                double cost = 0.0;
                for (const auto &edge : circuit_edges)
                {
                    IdxType la = edge.first;
                    IdxType lb = edge.second;
                    double weight = 1.0;
                    auto weight_it = circuit_edge_weights.find(edge);
                    if (weight_it != circuit_edge_weights.end())
                    {
                        weight = weight_it->second;
                    }
                    double dist = distance_lookup(mapping[static_cast<std::size_t>(la)],
                                                  mapping[static_cast<std::size_t>(lb)]);
                    cost += weight * dist;
                }
                return cost;
            };

            double current_cost = total_cost();
            std::size_t current_twoq_depth = compute_twoq_depth(mapping);
            const int max_passes = 2;
            const std::size_t max_neighbors = 4;
            int swaps_applied = 0;
            for (int pass = 0; pass < max_passes; ++pass)
            {
                bool improved = false;
                for (IdxType logical = 0; logical < logical_qubits; ++logical)
                {
                    auto &neighbors = adjacency[static_cast<std::size_t>(logical)];
                    if (neighbors.empty())
                    {
                        continue;
                    }
                    std::sort(neighbors.begin(), neighbors.end(),
                              [](const auto &a, const auto &b)
                              { return a.second > b.second; });
                    std::size_t count = 0;
                    for (const auto &entry : neighbors)
                    {
                        if (count++ >= max_neighbors)
                        {
                            break;
                        }
                        IdxType other = entry.first;
                        if (other < 0 || other >= logical_qubits || other == logical)
                        {
                            continue;
                        }
                        IdxType phys_a = mapping[static_cast<std::size_t>(logical)];
                        IdxType phys_b = mapping[static_cast<std::size_t>(other)];
                        if (phys_a == phys_b)
                        {
                            continue;
                        }

                        double delta = 0.0;
                        for (const auto &adj : adjacency[static_cast<std::size_t>(logical)])
                        {
                            IdxType neighbor = adj.first;
                            if (neighbor == other)
                            {
                                continue;
                            }
                            double weight = adj.second;
                            IdxType phys_neighbor = mapping[static_cast<std::size_t>(neighbor)];
                            double old_dist = distance_lookup(phys_a, phys_neighbor);
                            double new_dist = distance_lookup(phys_b, phys_neighbor);
                            delta += weight * (new_dist - old_dist);
                        }
                        for (const auto &adj : adjacency[static_cast<std::size_t>(other)])
                        {
                            IdxType neighbor = adj.first;
                            if (neighbor == logical)
                            {
                                continue;
                            }
                            double weight = adj.second;
                            IdxType phys_neighbor = mapping[static_cast<std::size_t>(neighbor)];
                            double old_dist = distance_lookup(phys_b, phys_neighbor);
                            double new_dist = distance_lookup(phys_a, phys_neighbor);
                            delta += weight * (new_dist - old_dist);
                        }
                        if (delta < -1e-9)
                        {
                            std::swap(mapping[static_cast<std::size_t>(logical)],
                                      mapping[static_cast<std::size_t>(other)]);
                            std::size_t candidate_depth = compute_twoq_depth(mapping);
                            if (candidate_depth < current_twoq_depth)
                            {
                                current_twoq_depth = candidate_depth;
                                current_cost += delta;
                                improved = true;
                                swaps_applied += 1;
                            }
                            else
                            {
                                std::swap(mapping[static_cast<std::size_t>(logical)],
                                          mapping[static_cast<std::size_t>(other)]);
                            }
                        }
                    }
                }
                if (!improved)
                {
                    break;
                }
            }

            if (debug_level > 0)
            {
                std::cout << "Fast routing_mapping refinement: swaps=" << swaps_applied
                          << " cost=" << current_cost << " twoq_depth=" << current_twoq_depth << std::endl;
            }
        }

        (void)enable_sabre_layout;

        if (debug_level > 0)
        {
            std::cout << "Fast routing_mapping: vf2_embeddings="
                      << (found_embedding ? explored : 0)
                      << " best_cost=" << best_cost
                      << " (routing via sabre heuristic)" << std::endl;
        }
        circuit->set_mapping(mapping);
        Routing(circuit, chip, debug_level, enable_decay, decay_increment,
                decay_reset, routing_trials, routing_seed);
    }
} // namespace QASMTrans
