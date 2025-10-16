#pragma once

#include <algorithm>
#include <cctype>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <array>
#include <chrono>
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

enum class CriticalPathHeuristicMode
{
    LogProduct,
    AdditiveAverage
};


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

    inline double lookup_gate_length(const shared_ptr<Chip> &chip, const Gate &gate)
    {
        if (!chip)
        {
            return 0.0;
        }

        const std::string gate_name = to_lower_copy(OP_NAMES[gate.op_name]);
        auto fetch_single = [&](IdxType qubit, const std::string &primary) -> double
        {
            if (qubit >= 0 &&
                qubit < static_cast<IdxType>(chip->single_qubit_gate_lengths.size()))
            {
                const auto &length_map = chip->single_qubit_gate_lengths[qubit];
                auto it = length_map.find(primary);
                if (it != length_map.end())
                {
                    return it->second;
                }
                static const std::array<const char *, 3> single_fallbacks = {"sx", "x", "id"};
                for (const char *fallback : single_fallbacks)
                {
                    if (!fallback || primary == fallback)
                    {
                        continue;
                    }
                    it = length_map.find(fallback);
                    if (it != length_map.end())
                    {
                        return it->second;
                    }
                }
            }
            return 0.0;
        };

        auto fetch_two_qubit = [&](IdxType ctrl, IdxType tgt,
                                   const std::string &primary) -> double
        {
            std::pair<IdxType, IdxType> key{ctrl, tgt};
            auto len_it = chip->two_qubit_gate_lengths.find(key);
            if (len_it == chip->two_qubit_gate_lengths.end())
            {
                key = {tgt, ctrl};
                len_it = chip->two_qubit_gate_lengths.find(key);
            }
            if (len_it != chip->two_qubit_gate_lengths.end())
            {
                auto gate_it = len_it->second.find(primary);
                if (gate_it != len_it->second.end())
                {
                    return gate_it->second;
                }
                static const std::array<const char *, 1> two_fallbacks = {"cx"};
                for (const char *fallback : two_fallbacks)
                {
                    if (!fallback || primary == fallback)
                    {
                        continue;
                    }
                    gate_it = len_it->second.find(fallback);
                    if (gate_it != len_it->second.end())
                    {
                        return gate_it->second;
                    }
                }
            }
            return 0.0;
        };

        if (gate.op_name == OP::SWAP)
        {
            if (gate.ctrl >= 0 && gate.qubit >= 0)
            {
                double cx_length =
                    fetch_two_qubit(gate.ctrl, gate.qubit, std::string("cx"));
                if (cx_length > 0.0)
                {
                    return 3.0 * cx_length;
                }
            }
            if (gate.ctrl >= 0 && gate.qubit >= 0)
            {
                double swap_length =
                    fetch_two_qubit(gate.ctrl, gate.qubit, gate_name);
                if (swap_length > 0.0)
                {
                    return swap_length;
                }
            }
            return 0.0;
        }

        if (gate.ctrl >= 0 && gate.qubit >= 0)
        {
            return fetch_two_qubit(gate.ctrl, gate.qubit, gate_name);
        }

        if (gate.qubit >= 0)
        {
            return fetch_single(gate.qubit, gate_name);
        }

        if (gate.extra >= 0)
        {
            return fetch_single(gate.extra, gate_name);
        }

        return 0.0;
    }

    struct CriticalPathResult
    {
        std::vector<IdxType> gate_indices;
        double latency = 0.0;
    };

    inline CriticalPathResult compute_critical_path(const std::vector<Gate> &gates,
                                                    const shared_ptr<Chip> &chip)
    {
        CriticalPathResult result;
        if (gates.empty())
        {
            return result;
        }

        IdxType qubit_capacity = chip ? chip->chip_qubit_num : 0;
        if (qubit_capacity <= 0)
        {
            for (const auto &gate : gates)
            {
                qubit_capacity = std::max(qubit_capacity, gate.qubit + 1);
                qubit_capacity = std::max(qubit_capacity, gate.ctrl + 1);
                qubit_capacity = std::max(qubit_capacity, gate.extra + 1);
            }
        }

        std::vector<double> ready_time(static_cast<std::size_t>(qubit_capacity), 0.0);
        std::vector<IdxType> last_gate(static_cast<std::size_t>(qubit_capacity), -1);
        std::vector<double> finish_times;
        finish_times.reserve(gates.size());
        std::vector<IdxType> predecessors;
        predecessors.reserve(gates.size());

        double max_finish_time = 0.0;
        IdxType max_finish_index = -1;

        auto ensure_capacity = [&](IdxType qubit)
        {
            if (qubit < 0)
            {
                return;
            }
            if (qubit >= static_cast<IdxType>(ready_time.size()))
            {
                ready_time.resize(static_cast<std::size_t>(qubit + 1), 0.0);
                last_gate.resize(static_cast<std::size_t>(qubit + 1), -1);
            }
        };

        for (IdxType idx = 0; idx < static_cast<IdxType>(gates.size()); ++idx)
        {
            const auto &gate = gates[static_cast<std::size_t>(idx)];
            std::vector<IdxType> touched_qubits;
            touched_qubits.reserve(3);
            if (gate.ctrl >= 0)
            {
                ensure_capacity(gate.ctrl);
                touched_qubits.push_back(gate.ctrl);
            }
            if (gate.qubit >= 0)
            {
                ensure_capacity(gate.qubit);
                touched_qubits.push_back(gate.qubit);
            }
            if (gate.extra >= 0)
            {
                ensure_capacity(gate.extra);
                touched_qubits.push_back(gate.extra);
            }

            std::sort(touched_qubits.begin(), touched_qubits.end());
            touched_qubits.erase(std::unique(touched_qubits.begin(), touched_qubits.end()), touched_qubits.end());

            double duration = lookup_gate_length(chip, gate);

            double start_time = 0.0;
            double predecessor_finish = -1.0;
            IdxType predecessor = -1;
            for (IdxType qubit : touched_qubits)
            {
                double ready = ready_time[static_cast<std::size_t>(qubit)];
                if (ready > start_time)
                {
                    start_time = ready;
                }
                IdxType last = last_gate[static_cast<std::size_t>(qubit)];
                if (last >= 0)
                {
                    double finish = finish_times[static_cast<std::size_t>(last)];
                    if (finish > predecessor_finish)
                    {
                        predecessor_finish = finish;
                        predecessor = last;
                    }
                }
            }

            double end_time = start_time + duration;
            finish_times.push_back(end_time);
            predecessors.push_back(predecessor);
            if (end_time >= max_finish_time)
            {
                max_finish_time = end_time;
                max_finish_index = idx;
            }

            for (IdxType qubit : touched_qubits)
            {
                ready_time[static_cast<std::size_t>(qubit)] = end_time;
                last_gate[static_cast<std::size_t>(qubit)] = idx;
            }
        }

        if (max_finish_index >= 0)
        {
            IdxType current = max_finish_index;
            while (current >= 0)
            {
                result.gate_indices.push_back(current);
                current = predecessors[static_cast<std::size_t>(current)];
            }
            std::reverse(result.gate_indices.begin(), result.gate_indices.end());
            result.latency = max_finish_time;
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
                                        const shared_ptr<Chip> &chip,
                                        const std::vector<IdxType> &critical_path,
                                        CriticalPathHeuristicMode mode)
    {
        if (!chip)
        {
            return mode == CriticalPathHeuristicMode::LogProduct ? 0.0 : -1.0;
        }

        constexpr double kMinComponent = 1e-15;
        double log_fidelity = 0.0;
        double additive_sum = 0.0;
        std::size_t processed_gate_count = 0;

        auto accumulate_gate = [&](const Gate &gate)
        {
            if (strcmp(OP_NAMES[gate.op_name], "MA") == 0)
            {
                return;
            }

            const std::string gate_name = to_lower_copy(OP_NAMES[gate.op_name]);
            double gate_error = 0.0;

            auto consider_two_qubit = [&](IdxType ctrl, IdxType tgt)
            {
                const IdxType mapped_ctrl = remap_if_present(ctrl, physical_mapping);
                const IdxType mapped_tgt = remap_if_present(tgt, physical_mapping);

                auto pair_it = chip->two_qubit_errors.find({mapped_ctrl, mapped_tgt});
                if (pair_it == chip->two_qubit_errors.end())
                {
                    pair_it = chip->two_qubit_errors.find({mapped_tgt, mapped_ctrl});
                }

                if (pair_it != chip->two_qubit_errors.end())
                {
                    const auto &error_map = pair_it->second;
                    auto err_it = error_map.find(gate_name);
                    if (err_it != error_map.end())
                    {
                        gate_error = std::clamp(err_it->second, 0.0, 1.0);
                        return;
                    }
                    static const std::array<const char *, 1> kTwoFallbacks = {"cx"};
                    for (const char *fallback : kTwoFallbacks)
                    {
                        if (!fallback || gate_name == fallback)
                        {
                            continue;
                        }
                        err_it = error_map.find(fallback);
                        if (err_it != error_map.end())
                        {
                            gate_error = std::clamp(err_it->second, 0.0, 1.0);
                            return;
                        }
                    }
                }
            };

            auto consider_single_qubit = [&](IdxType qubit)
            {
                const IdxType mapped_qubit = remap_if_present(qubit, physical_mapping);
                if (mapped_qubit >= 0 && mapped_qubit < static_cast<IdxType>(chip->single_qubit_errors.size()))
                {
                    const auto &error_map = chip->single_qubit_errors[mapped_qubit];
                    auto err_it = error_map.find(gate_name);
                    if (err_it != error_map.end())
                    {
                        gate_error = std::clamp(err_it->second, 0.0, 1.0);
                        return;
                    }
                    static const std::array<const char *, 3> kSingleFallbacks = {"sx", "x", "id"};
                    for (const char *fallback : kSingleFallbacks)
                    {
                        if (!fallback || gate_name == fallback)
                        {
                            continue;
                        }
                        err_it = error_map.find(fallback);
                        if (err_it != error_map.end())
                        {
                            gate_error = std::clamp(err_it->second, 0.0, 1.0);
                            return;
                        }
                    }
                }
            };

            if (gate.ctrl >= 0 && gate.qubit >= 0)
            {
                consider_two_qubit(gate.ctrl, gate.qubit);
            }
            else if (gate.qubit >= 0)
            {
                consider_single_qubit(gate.qubit);
            }

            double fidelity_component = 1.0 - gate_error;
            if (mode == CriticalPathHeuristicMode::LogProduct)
            {
                fidelity_component = std::max(fidelity_component, kMinComponent);
                log_fidelity += std::log(fidelity_component);
            }
            else
            {
                fidelity_component = std::clamp(fidelity_component, 0.0, 1.0);
                additive_sum += fidelity_component;
            }
            ++processed_gate_count;
        };

        if (!critical_path.empty())
        {
            for (IdxType index : critical_path)
            {
                if (index < 0 || static_cast<std::size_t>(index) >= gates.size())
                {
                    continue;
                }
                accumulate_gate(gates[static_cast<std::size_t>(index)]);
            }
        }
        else
        {
            for (const auto &gate : gates)
            {
                accumulate_gate(gate);
            }
        }

        if (processed_gate_count == 0)
        {
            return mode == CriticalPathHeuristicMode::LogProduct ? 0.0 : -1.0;
        }
        if (mode == CriticalPathHeuristicMode::LogProduct)
        {
            return -log_fidelity;
        }
        double average = additive_sum / static_cast<double>(processed_gate_count);
        average = std::clamp(average, 0.0, 1.0);
        return -average;
    }
}


void calibration_aware_optimization(shared_ptr<Circuit> circuit,
                  shared_ptr<Chip> chip,
                  IdxType debug_level,
                  bool use_full_fidelity,
                  CriticalPathHeuristicMode cp_mode)
{
    (void)debug_level;
    if (!circuit || !chip)
    {
        return;
    }

    double critical_path_prep_ms = 0.0;
    double embedding_time_ms = 0.0;
    double scoring_time_ms = 0.0;
    double apply_time_ms = 0.0;
    const auto critical_path_start = std::chrono::steady_clock::now();

    constexpr double kLatencyEpsilon = 1e-12;
    std::vector<Gate> gates = circuit->get_gates();
    std::vector<IdxType> critical_path_gates = circuit->get_critical_path();
    double critical_path_latency = circuit->get_critical_path_latency();
    bool critical_path_from_circuit = !critical_path_gates.empty() && critical_path_latency > kLatencyEpsilon;
    if (!critical_path_from_circuit)
    {
        auto computed = mapomatic_detail::compute_critical_path(gates, chip);
        critical_path_gates = std::move(computed.gate_indices);
        critical_path_latency = computed.latency;
    }
    if (critical_path_latency <= kLatencyEpsilon)
    {
        critical_path_gates.clear();
        critical_path_latency = 0.0;
    }
    std::vector<IdxType> empty_scoring_path;
    const std::vector<IdxType> &path_for_scoring = use_full_fidelity ? empty_scoring_path : critical_path_gates;
    if (gates.empty())
    {
        return;
    }

    critical_path_prep_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - critical_path_start).count();

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

    const auto embedding_start = std::chrono::steady_clock::now();

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
    bool vf2_found_embedding = false;

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
        vf2_found_embedding = true;
    }

    bool identity_fallback_used = false;
    if (candidate_mappings.empty())
    {
        std::unordered_map<IdxType, IdxType> identity;
        identity.reserve(used_qubits.size());
        for (IdxType qubit : used_qubits)
        {
            identity.emplace(qubit, qubit);
        }
        candidate_mappings.push_back(std::move(identity));
        identity_fallback_used = true;
    }

    embedding_time_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - embedding_start).count();

    struct ScoreSummary
    {
        double best_score = std::numeric_limits<double>::max();
        double worst_score = std::numeric_limits<double>::lowest();
        std::size_t best_index = 0;
    };

    std::vector<double> candidate_scores;
    candidate_scores.reserve(candidate_mappings.size());

    constexpr double kScoreTolerance = 1e-12;

    auto evaluate_scores = [&](const std::vector<IdxType> &path, bool log_candidates) -> ScoreSummary
    {
        candidate_scores.clear();
        candidate_scores.reserve(candidate_mappings.size());
        ScoreSummary summary;
        for (std::size_t idx = 0; idx < candidate_mappings.size(); ++idx)
        {
            double score = mapomatic_detail::calibration_heuristic(gates, candidate_mappings[idx], chip, path, cp_mode);
            if (log_candidates && idx < 5)
            {
                std::vector<std::pair<IdxType, IdxType>> mapping_pairs(candidate_mappings[idx].begin(), candidate_mappings[idx].end());
                std::sort(mapping_pairs.begin(), mapping_pairs.end());
                double display = (cp_mode == CriticalPathHeuristicMode::LogProduct)
                                     ? std::clamp(std::exp(-score), 0.0, 1.0)
                                     : std::clamp(-score, 0.0, 1.0);
                const char *score_label = (cp_mode == CriticalPathHeuristicMode::LogProduct) ? "log-score" : "avg-score";
                std::cout << std::scientific << std::setprecision(6);
                std::cout << "  Candidate " << idx << " " << score_label << ": " << score;
                std::cout << std::fixed << std::setprecision(6);
                std::cout << " fidelity: " << display << " mapping:";
                for (const auto &entry : mapping_pairs)
                {
                    std::cout << " " << entry.first << "->" << entry.second;
                }
                std::cout << std::endl;
            }
            candidate_scores.push_back(score);
            if (score < summary.best_score)
            {
                summary.best_score = score;
                summary.best_index = idx;
            }
            if (score > summary.worst_score)
            {
                summary.worst_score = score;
            }
        }
        return summary;
    };

    const bool initial_logging = debug_level > 1 && !use_full_fidelity;
    const auto scoring_start = std::chrono::steady_clock::now();
    ScoreSummary summary = evaluate_scores(path_for_scoring, initial_logging);
    bool forced_full_fidelity = false;
    if (!use_full_fidelity &&
        cp_mode == CriticalPathHeuristicMode::LogProduct &&
        std::abs(summary.worst_score - summary.best_score) < kScoreTolerance)
    {
        forced_full_fidelity = true;
        summary = evaluate_scores(empty_scoring_path, debug_level > 1);
        if (debug_level > 0)
        {
            std::cout << "Critical path scores were identical; falling back to full circuit fidelity for scoring." << std::endl;
        }
    }
    scoring_time_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - scoring_start).count();

    const auto &best_mapping = candidate_mappings[summary.best_index];

    const auto apply_start = std::chrono::steady_clock::now();

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
        auto recomputed = mapomatic_detail::compute_critical_path(updated_gates, chip);
        if (recomputed.latency <= kLatencyEpsilon)
        {
            recomputed.gate_indices.clear();
            recomputed.latency = 0.0;
        }
        circuit->set_critical_path(recomputed.gate_indices, recomputed.latency);
        if (!recomputed.gate_indices.empty() || recomputed.latency > 0.0)
        {
            critical_path_gates = std::move(recomputed.gate_indices);
            critical_path_latency = recomputed.latency;
        }

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
    else if (!critical_path_from_circuit &&
             (!critical_path_gates.empty() || critical_path_latency > 0.0))
    {
        circuit->set_critical_path(critical_path_gates, critical_path_latency);
    }

    apply_time_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - apply_start).count();

    if (debug_level > 0)
    {
        const bool used_full_circuit = use_full_fidelity || critical_path_gates.empty() || forced_full_fidelity;
        auto score_to_fidelity = [&](double score) -> double
        {
            if (cp_mode == CriticalPathHeuristicMode::LogProduct)
            {
                if (!std::isfinite(score))
                {
                    return 0.0;
                }
                double fidelity = std::exp(-score);
                if (!std::isfinite(fidelity))
                {
                    return 0.0;
                }
                return std::clamp(fidelity, 0.0, 1.0);
            }
            return std::clamp(-score, 0.0, 1.0);
        };
        double best_fidelity = score_to_fidelity(summary.best_score);
        double worst_fidelity = score_to_fidelity(summary.worst_score);
        std::cout << std::fixed << std::setprecision(6)
                  << "Calibration-aware heuristic "
                  << (used_full_circuit ? "(full circuit fidelity)" : "(critical path fidelity)")
                  << ": best " << best_fidelity << ", worst " << worst_fidelity << std::endl;
        double score_gap = summary.worst_score - summary.best_score;
        const char *spread_label = (cp_mode == CriticalPathHeuristicMode::LogProduct) ? "Log-fidelity spread" : "Average-fidelity spread";
        std::cout << spread_label << ": " << std::scientific << std::setprecision(3) << score_gap << std::fixed << std::setprecision(6) << std::endl;
        std::cout << "Critical path aggregation: "
                  << (cp_mode == CriticalPathHeuristicMode::LogProduct ? "log product" : "additive average") << std::endl;
        std::cout << std::fixed << std::setprecision(3)
                  << "Mapomatic timing (ms): critical-path prep " << critical_path_prep_ms
                  << ", embedding " << embedding_time_ms
                  << ", scoring " << scoring_time_ms
                  << ", apply " << apply_time_ms << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << std::fixed << std::setprecision(6)
                  << "Critical path latency: " << critical_path_latency << " s" << std::endl;
        std::cout << "Mapomatic evaluated " << candidate_mappings.size()
                  << " candidate mapping(s)"
                  << (identity_fallback_used ? " (identity fallback)" : vf2_found_embedding ? "" : " (VF2 search returned no embeddings)") << "." << std::endl;

        std::vector<double> dedup_scores;
        std::size_t unique_score_count = 0;
        if (!candidate_scores.empty())
        {
            dedup_scores = candidate_scores;
            std::sort(dedup_scores.begin(), dedup_scores.end());
            dedup_scores.erase(std::unique(dedup_scores.begin(), dedup_scores.end(),
                                           [&](double a, double b)
                                           { return std::abs(a - b) < kScoreTolerance; }),
                               dedup_scores.end());
            unique_score_count = dedup_scores.size();
            if (unique_score_count > 1 && debug_level > 1)
            {
                std::cout << "  Fidelity samples:";
                for (std::size_t i = 0; i < std::min<std::size_t>(dedup_scores.size(), static_cast<std::size_t>(5)); ++i)
                {
                    std::cout << " " << std::fixed << std::setprecision(6) << score_to_fidelity(dedup_scores[i]);
                }
                if (dedup_scores.size() > 5)
                {
                    std::cout << " ...";
                }
                std::cout << std::endl;
            }
        }
        if (unique_score_count <= 1)
        {
            std::cout << "Best and worst heuristic values match because all evaluated mappings scored identically." << std::endl;
        }
    }
}
