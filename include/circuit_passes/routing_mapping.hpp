#pragma once

#include <algorithm>
#include <cctype>
#include <chrono>
#include <deque>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "../QASMTransPrimitives.hpp"

#include "../IR/chip.hpp"
#include "../IR/circuit.hpp"
#include "../IR/gate.hpp"
#include "../IR/graph.hpp"

#include "../nlomann/json.hpp"

using namespace QASMTrans;
using namespace std;
using json = nlohmann::json;

constexpr double kSabreBasicWeight = 1.0;
constexpr double kSabreLookaheadWeight = 0.5;
constexpr IdxType kSabreLookaheadSize = 20;
constexpr double kSabreDirectionPenalty = 1.0;
constexpr double kSabreProgressWeight = 0.2;

// extract cx gates in json file for constructing graph
vector<pair<IdxType, IdxType>> extract_cx_pairs(const json &j) {
  vector<pair<IdxType, IdxType>> pairs;
  for (auto &item : j.items()) {
    string key = item.key();
    if (key.substr(0, 2) == "cx") {
      IdxType pos = key.find('_');
      if (pos != string::npos) {
        IdxType first = stoi(key.substr(2, pos - 2));
        IdxType second = stoi(key.substr(pos + 1));
        pairs.push_back(make_pair(first, second));
      }
    }
  }
  return pairs;
}

inline std::string to_lower_copy(const std::string &value) {
  std::string result;
  result.reserve(value.size());
  for (char ch : value) {
    result.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
  }
  return result;
}

inline double lookup_gate_length(const shared_ptr<Chip> &chip, const Gate &gate) {
  if (!chip) {
    return 0.0;
  }

  const std::string gate_name = to_lower_copy(OP_NAMES[gate.op_name]);

  auto fetch_single = [&](IdxType qubit) -> double {
    if (qubit >= 0 &&
        qubit < static_cast<IdxType>(chip->single_qubit_gate_lengths.size())) {
      const auto &length_map = chip->single_qubit_gate_lengths[qubit];
      auto it = length_map.find(gate_name);
      if (it != length_map.end()) {
        return it->second;
      }
    }
    return 0.0;
  };

  auto fetch_two_qubit = [&](IdxType ctrl, IdxType tgt,
                             const std::string &name) -> double {
    std::pair<IdxType, IdxType> key{ctrl, tgt};
    auto len_it = chip->two_qubit_gate_lengths.find(key);
    if (len_it == chip->two_qubit_gate_lengths.end()) {
      key = {tgt, ctrl};
      len_it = chip->two_qubit_gate_lengths.find(key);
    }
    if (len_it != chip->two_qubit_gate_lengths.end()) {
      auto gate_it = len_it->second.find(name);
      if (gate_it != len_it->second.end()) {
        return gate_it->second;
      }
    }
    return 0.0;
  };

  if (gate.op_name == OP::SWAP) {
    if (gate.ctrl >= 0 && gate.qubit >= 0) {
      double cx_length =
          fetch_two_qubit(gate.ctrl, gate.qubit, std::string("cx"));
      if (cx_length > 0.0) {
        return 3.0 * cx_length;
      }
    }
    // fall back to any explicit SWAP entry if present
    if (gate.ctrl >= 0 && gate.qubit >= 0) {
      double swap_length =
          fetch_two_qubit(gate.ctrl, gate.qubit, gate_name);
      if (swap_length > 0.0) {
        return swap_length;
      }
    }
    return 0.0;
  }

  if (gate.ctrl >= 0 && gate.qubit >= 0) {
    return fetch_two_qubit(gate.ctrl, gate.qubit, gate_name);
  }

  if (gate.qubit >= 0) {
    return fetch_single(gate.qubit);
  }

  if (gate.extra >= 0) {
    return fetch_single(gate.extra);
  }

  return 0.0;
}

void DAG_generator(IdxType qubit_num, vector<vector<IdxType>> &circuit,
                   vector<IdxType> &gate_state, vector<IdxType> &qubit_state,
                   vector<IdxType> &gate_dependency,
                   vector<vector<IdxType>> &following_gate_idx,
                   vector<IdxType> &first_layer_gates_idx) {
  IdxType gate_num = circuit.size();
  vector<IdxType> current_gate_idx(qubit_num, -1);
  following_gate_idx.resize(gate_num, vector<IdxType>(2, 0));
  gate_dependency.resize(gate_num, 0);
  for (IdxType i = 0; i < gate_num; i++) {
    vector<IdxType> gate = circuit[i];
    const auto valid_qubit = [&](IdxType q) {
      return q >= 0 && q < qubit_num;
    };

      if (valid_qubit(gate[0]) && current_gate_idx[gate[0]] == -1) {
        if (valid_qubit(gate[1]) && current_gate_idx[gate[1]] == -1) {
          first_layer_gates_idx.push_back(i);
          gate_state[i] = 2;
          if (valid_qubit(gate[0])) {
            qubit_state[gate[0]] = 1;
          }
          if (valid_qubit(gate[1])) {
            qubit_state[gate[1]] = 1;
          }
          gate_dependency[i] = 0;
        } else {
          gate_dependency[i] = 1;
        }
      }
      if (valid_qubit(gate[1]) && current_gate_idx[gate[1]] == -1 &&
          valid_qubit(gate[0]) && current_gate_idx[gate[0]] != -1) {
        gate_dependency[i] = 1;
      }
      for (IdxType j = 0; j < gate.size(); j++) {
        IdxType qubit = gate[j];
        if (!valid_qubit(qubit)) {
          continue;
        }
        if (current_gate_idx[qubit] != -1) {
          vector<IdxType> prior_gate = circuit[current_gate_idx[qubit]];
          IdxType qubit_idx;
          if (prior_gate[j] != qubit) {
            qubit_idx = 1 - j;
        } else {
          qubit_idx = j;
        }
        following_gate_idx[current_gate_idx[qubit]][qubit_idx] = i;
      }
      current_gate_idx[qubit] = i;
    }
  }
}
// #gate_state
// # 0 - not considered
// # 1 - in future gate queue
// # 2 - in current gate layer
// # 3 - executed

// #qubit_state
// # 0 - not occupied in current layer
// # 1 - occupied in current layer

void maintain_layer(vector<IdxType> &current_layer_gates_idx,
                    set<IdxType> &gate_execute_idx_list,
                    vector<vector<IdxType>> &circuit,
                    vector<IdxType> &gate_state,
                    vector<vector<IdxType>> &following_gate_idx,
                    vector<IdxType> &qubit_state,
                    vector<IdxType> &gate_dependency,
                    vector<IdxType> &updated_layer_gates_idx,
                    vector<IdxType> &future_layer_gates_idx, IdxType flag) {
  unordered_set<IdxType> updated_set;
  updated_layer_gates_idx.clear();
  auto valid_qubit = [&](IdxType q) {
    return q >= 0 && q < static_cast<IdxType>(qubit_state.size());
  };
  for (IdxType gate_idx : current_layer_gates_idx) {
    if (gate_execute_idx_list.count(gate_idx) > 0) {
      vector<IdxType> gate = circuit[gate_idx];
      gate_state[gate_idx] = 3;
      future_layer_gates_idx.erase(remove(future_layer_gates_idx.begin(),
                                          future_layer_gates_idx.end(),
                                          gate_idx),
                                       future_layer_gates_idx.end());

      if (valid_qubit(gate[0])) {
        qubit_state[gate[0]] = 0;
      }
      if (valid_qubit(gate[1])) {
        qubit_state[gate[1]] = 0;
      }
      vector<IdxType> following_gates = following_gate_idx[gate_idx];
      for (IdxType next_gate_idx : following_gates) {
        gate_dependency[next_gate_idx]--;
        if (gate_dependency[next_gate_idx] == 0) {
          updated_set.insert(next_gate_idx);
          gate_state[next_gate_idx] = 2;
          future_layer_gates_idx.erase(remove(future_layer_gates_idx.begin(),
                                              future_layer_gates_idx.end(),
                                              next_gate_idx),
                                       future_layer_gates_idx.end());
          if (valid_qubit(circuit[next_gate_idx][0])) {
            qubit_state[circuit[next_gate_idx][0]] = 1;
          }
          if (valid_qubit(circuit[next_gate_idx][1])) {
            qubit_state[circuit[next_gate_idx][1]] = 1;
          }
        }
      }
    } else {
      updated_set.insert(gate_idx);
    }
  }
  updated_layer_gates_idx.assign(updated_set.begin(), updated_set.end());
  if (!updated_layer_gates_idx.empty()) {
    IdxType start_gate = *min_element(updated_layer_gates_idx.begin(),
                                      updated_layer_gates_idx.end());
    for (IdxType gate_idx = start_gate;
         gate_idx < start_gate + 20 && gate_idx < circuit.size(); gate_idx++) {
      if (gate_state[gate_idx] == 0) {
        gate_state[gate_idx] = 1;
        if (flag != 0) {
          future_layer_gates_idx.push_back(gate_idx);
        }
      }
    }
  }
  if (flag == 0) {
    for (IdxType gate_idx = 0; gate_idx < circuit.size(); gate_idx++) {
      if (gate_state[gate_idx] == 1) {
        future_layer_gates_idx.push_back(gate_idx);
      }
    }
  }
  sort(updated_layer_gates_idx.begin(), updated_layer_gates_idx.end());
  sort(future_layer_gates_idx.begin(), future_layer_gates_idx.end());
}

double heuristic(const vector<IdxType> &new_mapping,
                 const vector<IdxType> &current_layer_gates_idx,
                 const vector<IdxType> &future_gates_idx,
                 const vector<vector<IdxType>> &distance_mat,
                 const vector<vector<IdxType>> &circuit,
                 const shared_ptr<Chip> &chip) {
  double cost = 0.0;
  double first_cost = 0.0;
  if (current_layer_gates_idx.empty()) {
    return 0;
  }
  const auto direction_allowed = [&](IdxType src, IdxType dst) -> bool {
    if (!chip || chip->directed_edge_list.empty()) {
      return true;
    }
    if (src < 0 || dst < 0 ||
        src >= static_cast<IdxType>(chip->directed_edge_list.size())) {
      return true;
    }
    const auto &targets = chip->directed_edge_list[static_cast<size_t>(src)];
    return targets.find(dst) != targets.end();
  };
  for (IdxType gate_idx : current_layer_gates_idx) {
    if (gate_idx < 0 || gate_idx >= static_cast<IdxType>(circuit.size())) {
      continue;
    }
    const auto &gate = circuit[gate_idx];
    if (gate.size() < 2) {
      continue;
    }
    IdxType ctrl = gate[0];
    IdxType tgt = gate[1];
    if (ctrl < 0 || ctrl >= static_cast<IdxType>(new_mapping.size()) ||
        tgt < 0 || tgt >= static_cast<IdxType>(new_mapping.size())) {
      continue;
    }
    IdxType mapped_ctrl = new_mapping[ctrl];
    IdxType mapped_tgt = new_mapping[tgt];
    if (mapped_ctrl < 0 ||
        mapped_ctrl >= static_cast<IdxType>(distance_mat.size()) ||
        mapped_tgt < 0 ||
        mapped_tgt >= static_cast<IdxType>(distance_mat[mapped_ctrl].size())) {
      continue;
    }
    double dist = distance_mat[mapped_ctrl][mapped_tgt];
    first_cost += dist;
    if (dist == 1 && !direction_allowed(mapped_ctrl, mapped_tgt)) {
      first_cost += kSabreDirectionPenalty;
    }
  }
  first_cost /= std::max<IdxType>(1, current_layer_gates_idx.size());
  if (future_gates_idx.empty()) {
    cost = first_cost;
    return cost;
  }
  double second_cost = 0.0;
  IdxType lookahead_limit = static_cast<IdxType>(
      std::min<std::size_t>(future_gates_idx.size(),
                            static_cast<std::size_t>(kSabreLookaheadSize)));
  for (IdxType i = 0; i < lookahead_limit; ++i) {
    IdxType gate_idx = future_gates_idx[static_cast<std::size_t>(i)];
    if (gate_idx < 0 || gate_idx >= static_cast<IdxType>(circuit.size())) {
      continue;
    }
    const auto &gate = circuit[gate_idx];
    if (gate.size() < 2) {
      continue;
    }
    IdxType ctrl = gate[0];
    IdxType tgt = gate[1];
    if (ctrl < 0 || ctrl >= static_cast<IdxType>(new_mapping.size()) ||
        tgt < 0 || tgt >= static_cast<IdxType>(new_mapping.size())) {
      continue;
    }
    IdxType mapped_ctrl = new_mapping[ctrl];
    IdxType mapped_tgt = new_mapping[tgt];
    if (mapped_ctrl < 0 ||
        mapped_ctrl >= static_cast<IdxType>(distance_mat.size()) ||
        mapped_tgt < 0 ||
        mapped_tgt >= static_cast<IdxType>(distance_mat[mapped_ctrl].size())) {
      continue;
    }
    double dist = distance_mat[mapped_ctrl][mapped_tgt];
    second_cost += dist;
    if (dist == 1 && !direction_allowed(mapped_ctrl, mapped_tgt)) {
      second_cost += kSabreDirectionPenalty;
    }
  }
  second_cost /= std::max<IdxType>(1, lookahead_limit);
  cost = kSabreBasicWeight * first_cost +
         kSabreLookaheadWeight * second_cost;
  return cost;
}

vector<IdxType> find_reverse_mapping(const vector<IdxType> &mapping,
                                     IdxType qubit_num) {
  vector<IdxType> reverse_mapping(qubit_num, -1);
  for (IdxType l_qubit = 0; l_qubit < mapping.size(); ++l_qubit) {
    IdxType p_qubit = mapping[l_qubit];
    if (p_qubit >= 0 && p_qubit < qubit_num) {
      reverse_mapping[p_qubit] = l_qubit;
    } else {
      // Handle the error in an appropriate way for your program.
      // cerr << "Invalid qubit index: " << p_qubit << endl;
    }
  }
  return reverse_mapping;
}

vector<IdxType> pick_one_movement(vector<IdxType> &mapping,
                                  const vector<IdxType> &current_layer,
                                  const vector<IdxType> &future_layer,
                                  const vector<vector<IdxType>> &distance_mat,
                                  IdxType qubit_num,
                                  const vector<vector<IdxType>> &circuit,
                                  shared_ptr<Chip> chip,
                                  std::vector<double> *decay,
                                  double decay_increment,
                                  IdxType decay_reset,
                                  IdxType &decay_steps,
                                  bool enable_decay,
                                  IdxType debug_level,
                                  const std::unordered_set<IdxType> *critical_gate_set = nullptr,
                                  const std::deque<std::pair<IdxType, IdxType>> *recent_swaps = nullptr) {
  vector<IdxType> l2p_mapping = mapping;
  vector<IdxType> key_p_qubits;
  for (IdxType gate_idx : current_layer) {
    if (gate_idx < 0 || gate_idx >= static_cast<IdxType>(circuit.size())) {
      continue;
    }
    vector<IdxType> gate = circuit[gate_idx];
    if (gate.size() < 2) {
      continue;
    }
    if (gate[0] >= 0 && gate[0] < static_cast<IdxType>(l2p_mapping.size())) {
      key_p_qubits.push_back(l2p_mapping[gate[0]]);
    }
    if (gate[1] >= 0 && gate[1] < static_cast<IdxType>(l2p_mapping.size())) {
      key_p_qubits.push_back(l2p_mapping[gate[1]]);
    }
  }
  if (key_p_qubits.empty()) {
    return {-1, -1};
  }
  vector<vector<IdxType>> possible_pairs;
  for (IdxType p_qubit : key_p_qubits) {
    if (!chip || p_qubit < 0 ||
        p_qubit >= static_cast<IdxType>(chip->edge_list.size())) {
      continue;
    }
    for (IdxType p_qubit_target : chip->edge_list[p_qubit]) {
      if (p_qubit_target < 0 ||
          p_qubit_target >= static_cast<IdxType>(chip->edge_list.size())) {
        continue;
      }
      possible_pairs.push_back({p_qubit, p_qubit_target});
    }
  }
  if (possible_pairs.empty()) {
    return {-1, -1};
  }
  const auto direction_allowed = [&](IdxType src, IdxType dst) -> bool {
    if (!chip || chip->directed_edge_list.empty()) {
      return true;
    }
    if (src < 0 || dst < 0 ||
        src >= static_cast<IdxType>(chip->directed_edge_list.size())) {
      return true;
    }
    const auto &targets = chip->directed_edge_list[static_cast<size_t>(src)];
    return targets.find(dst) != targets.end();
  };

  auto weighted_score = [&](const vector<IdxType> &layer,
                            const vector<IdxType> &l2p) -> double {
    if (layer.empty()) {
      return 0.0;
    }
    double acc = 0.0;
    double denom = 0.0;
    for (IdxType gate_idx : layer) {
      if (gate_idx < 0 || gate_idx >= static_cast<IdxType>(circuit.size())) {
        continue;
      }
      const auto &gate = circuit[gate_idx];
      if (gate.size() < 2) {
        continue;
      }
      IdxType ctrl = gate[0];
      IdxType tgt = gate[1];
      if (ctrl < 0 || tgt < 0 ||
          ctrl >= static_cast<IdxType>(l2p.size()) ||
          tgt >= static_cast<IdxType>(l2p.size())) {
        continue;
      }
      IdxType phys_ctrl = l2p[ctrl];
      IdxType phys_tgt = l2p[tgt];
      if (phys_ctrl < 0 || phys_tgt < 0 ||
          phys_ctrl >= static_cast<IdxType>(distance_mat.size()) ||
          phys_tgt >= static_cast<IdxType>(distance_mat[phys_ctrl].size())) {
        continue;
      }
      double weight = 1.0;
      if (critical_gate_set && critical_gate_set->count(gate_idx) > 0) {
        weight = 1.5;
      }
      double dist = static_cast<double>(distance_mat[phys_ctrl][phys_tgt]);
      acc += weight * dist;
      denom += weight;
      if (distance_mat[phys_ctrl][phys_tgt] == 1 &&
          !direction_allowed(phys_ctrl, phys_tgt)) {
        acc += weight * kSabreDirectionPenalty;
      }
    }
    if (denom <= 0.0) {
      return 0.0;
    }
    return acc / denom;
  };

  auto count_executable = [&](const vector<IdxType> &l2p) -> IdxType {
    IdxType count = 0;
    for (IdxType gate_idx : current_layer) {
      if (gate_idx < 0 || gate_idx >= static_cast<IdxType>(circuit.size())) {
        continue;
      }
      const auto &gate = circuit[gate_idx];
      if (gate.size() < 2) {
        continue;
      }
      IdxType ctrl = gate[0];
      IdxType tgt = gate[1];
      if (ctrl < 0 || tgt < 0 ||
          ctrl >= static_cast<IdxType>(l2p.size()) ||
          tgt >= static_cast<IdxType>(l2p.size())) {
        continue;
      }
      IdxType phys_ctrl = l2p[ctrl];
      IdxType phys_tgt = l2p[tgt];
      if (phys_ctrl < 0 || phys_tgt < 0 ||
          phys_ctrl >= static_cast<IdxType>(distance_mat.size()) ||
          phys_tgt >= static_cast<IdxType>(distance_mat[phys_ctrl].size())) {
        continue;
      }
      if (distance_mat[phys_ctrl][phys_tgt] == 1 &&
          direction_allowed(phys_ctrl, phys_tgt)) {
        count += 1;
      }
    }
    return count;
  };

  IdxType baseline_exec = count_executable(mapping);
  vector<double> score(possible_pairs.size(), 0.0);
  const auto is_tabu = [&](IdxType a, IdxType b) -> bool {
    if (!recent_swaps || recent_swaps->empty()) {
      return false;
    }
    auto ordered = std::minmax(a, b);
    for (const auto &p : *recent_swaps) {
      auto pord = std::minmax(p.first, p.second);
      if (pord.first == ordered.first && pord.second == ordered.second) {
        return true;
      }
    }
    return false;
  };
  const double tabu_penalty = 1000.0;
  for (size_t pair_idx = 0; pair_idx < possible_pairs.size(); ++pair_idx) {
    vector<IdxType> pair = possible_pairs[pair_idx];
    vector<IdxType> p2l_mapping = find_reverse_mapping(l2p_mapping, qubit_num);
    swap(p2l_mapping[pair[0]], p2l_mapping[pair[1]]);
    vector<IdxType> temp_l2p_mapping =
        find_reverse_mapping(p2l_mapping, qubit_num);
    double base_score = kSabreBasicWeight * weighted_score(current_layer, temp_l2p_mapping) +
                        kSabreLookaheadWeight * weighted_score(future_layer, temp_l2p_mapping);
    IdxType exec_after = count_executable(temp_l2p_mapping);
    double progress_bonus =
        static_cast<double>(exec_after - baseline_exec) * kSabreProgressWeight;
    base_score -= progress_bonus;
    if (enable_decay && decay != nullptr && !decay->empty()) {
      IdxType p0 = pair[0];
      IdxType p1 = pair[1];
      if (p0 >= 0 && p1 >= 0 &&
          p0 < static_cast<IdxType>(decay->size()) &&
          p1 < static_cast<IdxType>(decay->size())) {
        base_score *= std::max((*decay)[p0], (*decay)[p1]);
      }
    }
    if (is_tabu(pair[0], pair[1])) {
      base_score += tabu_penalty;
    }
    score[pair_idx] = base_score;
  }
  size_t best_move_idx =
      distance(score.begin(), min_element(score.begin(), score.end()));
  vector<IdxType> pair = possible_pairs[best_move_idx];
  if (debug_level > 1) {
    static IdxType pick_calls = 0;
    pick_calls += 1;
    if (pick_calls % 200 == 0) {
      vector<IdxType> p2l_mapping_dbg = find_reverse_mapping(l2p_mapping, qubit_num);
      swap(p2l_mapping_dbg[pair[0]], p2l_mapping_dbg[pair[1]]);
      vector<IdxType> temp_l2p_dbg = find_reverse_mapping(p2l_mapping_dbg, qubit_num);
      IdxType exec_after_dbg = count_executable(temp_l2p_dbg);
      cout << "[sabre] pick_one swap=(" << pair[0] << "," << pair[1]
           << ") score=" << fixed << setprecision(3) << score[best_move_idx]
           << " baseline_exec=" << baseline_exec
           << " exec_after=" << exec_after_dbg
           << " possible_pairs=" << possible_pairs.size()
           << endl;
    }
  }
  vector<IdxType> p2l_mapping = find_reverse_mapping(l2p_mapping, qubit_num);
  swap(p2l_mapping[pair[0]], p2l_mapping[pair[1]]);
  vector<IdxType> new_mapping = find_reverse_mapping(p2l_mapping, qubit_num);
  mapping = new_mapping;
  if (enable_decay && decay != nullptr && !decay->empty()) {
    decay_steps += 1;
    if (decay_reset > 0 && decay_steps >= decay_reset) {
      std::fill(decay->begin(), decay->end(), 1.0);
      decay_steps = 0;
    } else {
      if (pair[0] >= 0 && pair[0] < static_cast<IdxType>(decay->size())) {
        (*decay)[pair[0]] += decay_increment;
      }
      if (pair[1] >= 0 && pair[1] < static_cast<IdxType>(decay->size())) {
        (*decay)[pair[1]] += decay_increment;
      }
    }
  }
  return pair;
}

vector<IdxType> force_route_pair(vector<IdxType> &mapping,
                                 const vector<IdxType> &current_layer,
                                 const vector<vector<IdxType>> &distance_mat,
                                 const vector<vector<IdxType>> &circuit,
                                 IdxType qubit_num,
                                 shared_ptr<Chip> chip,
                                 IdxType debug_level) {
  if (!chip || current_layer.empty()) {
    return {-1, -1};
  }
  IdxType best_gate = -1;
  IdxType best_dist = std::numeric_limits<IdxType>::max();
  IdxType best_ctrl = -1;
  IdxType best_tgt = -1;
  for (IdxType gate_idx : current_layer) {
    if (gate_idx < 0 || gate_idx >= static_cast<IdxType>(circuit.size())) {
      continue;
    }
    const auto &gate = circuit[gate_idx];
    if (gate.size() < 2) {
      continue;
    }
    IdxType ctrl = gate[0];
    IdxType tgt = gate[1];
    if (ctrl < 0 || tgt < 0 || ctrl >= static_cast<IdxType>(mapping.size()) ||
        tgt >= static_cast<IdxType>(mapping.size())) {
      continue;
    }
    IdxType phys_ctrl = mapping[ctrl];
    IdxType phys_tgt = mapping[tgt];
    if (phys_ctrl < 0 || phys_tgt < 0 ||
        phys_ctrl >= static_cast<IdxType>(distance_mat.size()) ||
        phys_tgt >= static_cast<IdxType>(distance_mat[phys_ctrl].size())) {
      continue;
    }
    IdxType dist = distance_mat[phys_ctrl][phys_tgt];
    if (dist < best_dist) {
      best_dist = dist;
      best_gate = gate_idx;
      best_ctrl = ctrl;
      best_tgt = tgt;
    }
  }
  if (best_gate < 0 || best_dist <= 1) {
    return {-1, -1};
  }

  IdxType phys_ctrl = mapping[best_ctrl];
  IdxType phys_tgt = mapping[best_tgt];
  vector<IdxType> best_pair = {-1, -1};
  IdxType best_new_dist = best_dist;

  if (phys_ctrl >= 0 && phys_ctrl < static_cast<IdxType>(chip->edge_list.size())) {
    for (IdxType neighbor : chip->edge_list[phys_ctrl]) {
      if (neighbor < 0 || neighbor >= static_cast<IdxType>(distance_mat.size())) {
        continue;
      }
      IdxType new_dist = distance_mat[neighbor][phys_tgt];
      if (new_dist < best_new_dist) {
        best_new_dist = new_dist;
        best_pair = {phys_ctrl, neighbor};
      }
    }
  }

  if (phys_tgt >= 0 && phys_tgt < static_cast<IdxType>(chip->edge_list.size())) {
    for (IdxType neighbor : chip->edge_list[phys_tgt]) {
      if (neighbor < 0 || neighbor >= static_cast<IdxType>(distance_mat.size())) {
        continue;
      }
      IdxType new_dist = distance_mat[phys_ctrl][neighbor];
      if (new_dist < best_new_dist) {
        best_new_dist = new_dist;
        best_pair = {phys_tgt, neighbor};
      }
    }
  }

  if (best_pair[0] < 0 || best_pair[1] < 0) {
    if (phys_ctrl >= 0 && phys_ctrl < static_cast<IdxType>(chip->edge_list.size()) &&
        !chip->edge_list[phys_ctrl].empty()) {
      best_pair = {phys_ctrl, chip->edge_list[phys_ctrl][0]};
    }
  }

  if (best_pair[0] < 0 || best_pair[1] < 0) {
    return {-1, -1};
  }
  if (debug_level > 1) {
    cout << "[sabre] force_route gate=" << best_gate
         << " dist=" << best_dist
         << " pair=(" << best_pair[0] << "," << best_pair[1] << ")"
         << " new_dist=" << best_new_dist
         << endl;
  }

  vector<IdxType> p2l_mapping = find_reverse_mapping(mapping, qubit_num);
  swap(p2l_mapping[best_pair[0]], p2l_mapping[best_pair[1]]);
  mapping = find_reverse_mapping(p2l_mapping, qubit_num);
  return best_pair;
}

set<IdxType>
find_executable_gates(const vector<IdxType> &mapping,
                      const vector<IdxType> &current_layer,
                      const vector<vector<IdxType>> &circuit,
                      const vector<vector<IdxType>> &distance_mat,
                      const shared_ptr<Chip> &chip) {
  (void)chip;
  set<IdxType> executable_gates;
  // Pre-allocate memory using .reserve() for the worst-case scenario where
  // every gate is executable. executable_gates.reserve(current_layer.size());
  for (IdxType gate_idx : current_layer) {
    if (gate_idx < 0 || gate_idx >= static_cast<IdxType>(circuit.size())) {
      continue;
    }
    const auto &gate = circuit[gate_idx];
    if (gate.size() < 2) {
      continue;
    }
    IdxType ctrl_qubit = gate[0];
    IdxType tgt_qubit = gate[1];
    if (ctrl_qubit < 0 || ctrl_qubit >= static_cast<IdxType>(mapping.size()) ||
        tgt_qubit < 0 || tgt_qubit >= static_cast<IdxType>(mapping.size())) {
      continue;
    }
    IdxType mapped_gate_zero = mapping[ctrl_qubit];
    IdxType mapped_gate_one = mapping[tgt_qubit];
    if (mapped_gate_zero < 0 ||
        mapped_gate_zero >= static_cast<IdxType>(distance_mat.size()) ||
        mapped_gate_one < 0 ||
        mapped_gate_one >=
            static_cast<IdxType>(distance_mat[mapped_gate_zero].size())) {
      continue;
    }
    if (distance_mat[mapped_gate_zero][mapped_gate_one] == 1) {
      executable_gates.insert(gate_idx);
    }
  }
  // executable_gates.shrink_to_fit(); // Shrink the allocated memory to fit the
  // actual usage.
  return executable_gates;
}

vector<pair<IdxType, IdxType>> sortWithSwaps(vector<IdxType> &lst) {
  vector<IdxType> sorted_lst;
  vector<IdxType> temp_lst;
  for (IdxType x : lst) {
    if (x != -1) {
      temp_lst.push_back(x);
    }
  }
  for (IdxType x : lst) {
    if (x != -1) {
      sorted_lst.push_back(x);
    }
  }
  sort(sorted_lst.begin(), sorted_lst.end());
  vector<pair<IdxType, IdxType>> swaps;
  for (IdxType i = 0; i < lst.size(); i++) {
    if (lst[i] != -1 && lst[i] != sorted_lst[i]) {
      IdxType j = find(lst.begin() + i, lst.end(), sorted_lst[i]) - lst.begin();
      swap(lst[i], lst[j]);
      swaps.push_back(make_pair(lst[i], lst[j]));
    }
  }
  return swaps;
}

IdxType one_round_optimization(
    vector<IdxType> &initial_mapping, vector<Gate> circuit_gate,
    vector<vector<IdxType>> distance_mat, vector<Gate> gate_info,
    shared_ptr<Chip> chip, vector<vector<IdxType>> gate_qubit,
    vector<Gate> &return_circuit, IdxType debug_level,
    std::vector<IdxType> *critical_path_indices = nullptr,
    double *critical_path_latency = nullptr,
    bool enable_decay = false,
    double decay_increment = 0.001,
    IdxType decay_reset = 5) {
  IdxType swap_num = 0;
  vector<IdxType> mapping = initial_mapping;

  IdxType physical_qubit_count =
      chip ? chip->chip_qubit_num : static_cast<IdxType>(mapping.size());
  if (physical_qubit_count <= 0) {
    physical_qubit_count = static_cast<IdxType>(mapping.size());
  }
  std::vector<double> qubit_ready_time(physical_qubit_count, 0.0);
  std::vector<IdxType> qubit_last_gate(physical_qubit_count, -1);
  std::vector<double> gate_finish_times;
  std::vector<IdxType> gate_predecessor;
  double max_finish_time = 0.0;
  IdxType max_finish_index = -1;
  std::vector<double> decay;
  IdxType decay_steps = 0;
  if (enable_decay) {
    decay.assign(static_cast<size_t>(physical_qubit_count), 1.0);
  }

  auto record_gate = [&](const Gate &gate) {
    std::vector<IdxType> touched_qubits;
    if (gate.ctrl >= 0) {
      touched_qubits.push_back(gate.ctrl);
    }
    if (gate.qubit >= 0) {
      touched_qubits.push_back(gate.qubit);
    }
    if (gate.extra >= 0) {
      touched_qubits.push_back(gate.extra);
    }
    std::sort(touched_qubits.begin(), touched_qubits.end());
    touched_qubits.erase(
        std::unique(touched_qubits.begin(), touched_qubits.end()),
        touched_qubits.end());

    double duration = lookup_gate_length(chip, gate);
    double start_time = 0.0;
    IdxType predecessor = -1;
    double predecessor_finish = -1.0;
    for (IdxType phys_qubit : touched_qubits) {
      if (phys_qubit < 0 ||
          phys_qubit >= static_cast<IdxType>(qubit_ready_time.size())) {
        continue;
      }
      double ready = qubit_ready_time[phys_qubit];
      if (ready > start_time) {
        start_time = ready;
      }
      IdxType last_gate_idx = qubit_last_gate[phys_qubit];
      if (last_gate_idx >= 0) {
        double finish = gate_finish_times[last_gate_idx];
        if (finish > predecessor_finish) {
          predecessor_finish = finish;
          predecessor = last_gate_idx;
        }
      }
    }

    double end_time = start_time + duration;
    IdxType gate_index = static_cast<IdxType>(return_circuit.size());
    return_circuit.push_back(gate);
    gate_finish_times.push_back(end_time);
    gate_predecessor.push_back(predecessor);
    if (end_time >= max_finish_time) {
      max_finish_time = end_time;
      max_finish_index = gate_index;
    }
    for (IdxType phys_qubit : touched_qubits) {
      if (phys_qubit < 0 ||
          phys_qubit >= static_cast<IdxType>(qubit_ready_time.size())) {
        continue;
      }
      qubit_ready_time[phys_qubit] = end_time;
      qubit_last_gate[phys_qubit] = gate_index;
    }
  };

  //^find all single qubit dependency
  IdxType executed_gates_num = 0;
  IdxType gate_num = circuit_gate.size();
  vector<vector<IdxType>> circuit(gate_num, vector<IdxType>(2, 0));
  for (IdxType i = 0; i < gate_num; i++) {
    circuit[i][0] = circuit_gate[i].ctrl;
    circuit[i][1] = circuit_gate[i].qubit;
  }
  IdxType qubit_num = distance_mat.size();
  std::unordered_set<IdxType> critical_gate_set;
  if (!circuit.empty() && qubit_num > 0) {
    std::vector<IdxType> depth_per_qubit(qubit_num, 0);
    std::vector<IdxType> gate_depth(circuit.size(), 0);
    IdxType max_depth = 0;
    for (IdxType i = 0; i < static_cast<IdxType>(circuit.size()); ++i) {
      const auto &gate = circuit[i];
      if (gate.size() < 2) {
        continue;
      }
      IdxType ctrl = gate[0];
      IdxType tgt = gate[1];
      if (ctrl < 0 || tgt < 0 ||
          ctrl >= static_cast<IdxType>(depth_per_qubit.size()) ||
          tgt >= static_cast<IdxType>(depth_per_qubit.size())) {
        continue;
      }
      IdxType d = std::max(depth_per_qubit[ctrl], depth_per_qubit[tgt]) + 1;
      gate_depth[i] = d;
      depth_per_qubit[ctrl] = d;
      depth_per_qubit[tgt] = d;
      if (d > max_depth) {
        max_depth = d;
      }
    }
    if (max_depth > 0) {
      IdxType threshold = max_depth > 2 ? max_depth - 2 : max_depth - 1;
      for (IdxType i = 0; i < static_cast<IdxType>(gate_depth.size()); ++i) {
        if (gate_depth[i] >= threshold) {
          critical_gate_set.insert(i);
        }
      }
    }
  }
  vector<IdxType> gate_state(gate_num, 0);
  vector<IdxType> gate_dependency(gate_num, 2);
  vector<IdxType> qubit_state(qubit_num, 0);
  vector<vector<IdxType>> following_gates_idx;
  vector<IdxType> first_layer_gates_idx;
  DAG_generator(qubit_num, circuit, gate_state, qubit_state, gate_dependency,
                following_gates_idx, first_layer_gates_idx);
  vector<IdxType> current_layer;
  for (IdxType gate_idx : first_layer_gates_idx) {
    current_layer.push_back(gate_idx);
  }
  vector<IdxType> future_layer;
  set<IdxType> gate_execute_idx_list;
  vector<IdxType> updated_layer_gates_idx;
  maintain_layer(current_layer, gate_execute_idx_list, circuit, gate_state,
                 following_gates_idx, qubit_state, gate_dependency,
                 updated_layer_gates_idx, future_layer, 0);
  current_layer = updated_layer_gates_idx;
  IdxType layer_index = 0;
  IdxType single_gate_count = 0;
  vector<IdxType> num_single_before;
  vector<Gate> single_gate_info;
  for (IdxType i = 0; i < gate_info.size(); i++) {
    // cout << OP_NAMES[gate_info[i].op_name] << " (" << gate_info[i].ctrl << ",
    // " << gate_info[i].qubit <<")"<< endl;
    if (gate_info[i].ctrl == -1 &&
        strcmp(OP_NAMES[gate_info[i].op_name], "MA") != 0) {
      single_gate_count++;
      single_gate_info.push_back(gate_info[i]);
    } else {
      num_single_before.push_back(single_gate_count);
    }
  }
  map<IdxType, vector<IdxType>> qubit_to_gate_indices;
  vector<vector<IdxType>> dependency_vector(gate_num);

  IdxType two_qubit_gate_index = 0;
  IdxType single_qubit_index = 0;

  for (IdxType i = 0; i < gate_info.size(); ++i) {
    const auto &gate = gate_info[i];
    // If it's a two-qubit gate, we check whether the involved qubits were
    // touched by a single-qubit gate before
    if (gate.ctrl != -1) {
      // We look at all the previous single-qubit gates involving the control or
      // target qubits
      if (qubit_to_gate_indices.count(gate.ctrl)) {
        for (const auto &idx : qubit_to_gate_indices[gate.ctrl]) {
          dependency_vector[two_qubit_gate_index].push_back(idx);
        }
        qubit_to_gate_indices.erase(gate.ctrl);
      }
      if (qubit_to_gate_indices.count(gate.qubit)) {
        for (const auto &idx : qubit_to_gate_indices[gate.qubit]) {
          dependency_vector[two_qubit_gate_index].push_back(idx);
        }
        qubit_to_gate_indices.erase(gate.qubit);
      }
      two_qubit_gate_index++;
    }
    // We assume any other gate is a single-qubit gate
    else {
      if (strcmp(OP_NAMES[gate_info[i].op_name], "MA") != 0) {
        qubit_to_gate_indices[gate.qubit].push_back(single_qubit_index++);
      }
    }
  }
  double total_maIdxTypeainlayer_time = 0;
  double total_pickone_time = 0;
  IdxType no_progress_steps = 0;
  IdxType attempt_limit =
      std::max<IdxType>(10, qubit_num * 10);
  set<IdxType> visited_gate;
  IdxType cur = 0;
  IdxType iter_count = 0;
  const auto loop_start = std::chrono::steady_clock::now();
  std::deque<std::pair<IdxType, IdxType>> recent_swaps;
  const std::size_t recent_limit = 8;
  while (executed_gates_num < gate_num) {
    set<IdxType> execute_gates_idx =
        find_executable_gates(mapping, current_layer, circuit, distance_mat, chip);
    if (debug_level > 1 && (iter_count % 200 == 0)) {
      const auto now = std::chrono::steady_clock::now();
      const auto elapsed_ms =
          std::chrono::duration_cast<std::chrono::milliseconds>(now - loop_start).count();
      cout << "[sabre] iter=" << iter_count
           << " exec=" << executed_gates_num << "/" << gate_num
           << " cur=" << current_layer.size()
           << " fut=" << future_layer.size()
           << " ready=" << execute_gates_idx.size()
           << " swaps=" << swap_num
           << " no_progress=" << no_progress_steps
           << " elapsed=" << elapsed_ms << " ms"
           << endl;
    }
    // cout << current_layer.size()<<endl;
    for (IdxType ee : execute_gates_idx) {
      vector<IdxType> cur_index_vector = dependency_vector[ee];
      for (IdxType cur_index : cur_index_vector) {
        //^ push back all the single qubit gate
        Gate cur_gate = single_gate_info[cur_index];
        IdxType q_qubit = mapping[single_gate_info[cur_index].qubit];
        cur_gate.qubit = q_qubit;
        // Gate new_gate = Gate(cur_gate.op_name,
        // mapping[single_gate_info[single_gate_index].qubit], -1,
        // cur_gate.theta); new_gate.set_gm(cur_gate.gm_real, cur_gate.gm_imag,
        // 2);
        record_gate(cur_gate);
        visited_gate.insert(cur_index);
      }
      Gate cur_gate = circuit_gate[ee];
      IdxType q_qubit = mapping[cur_gate.qubit];
      IdxType c_qubit = mapping[cur_gate.ctrl];
      cur_gate.qubit = q_qubit;
      cur_gate.ctrl = c_qubit;
      // Gate new_gate = Gate(cur_gate.op_name, mapping[cur_gate.qubit],
      // mapping[cur_gate.ctrl], cur_gate.theta);
      // new_gate.set_gm(cur_gate.gm_real, cur_gate.gm_imag, 4);
      record_gate(cur_gate);
    }
    if (!execute_gates_idx.empty()) {
      cpu_timer trans_timer;
      trans_timer.start_timer();
      maintain_layer(current_layer, execute_gates_idx, circuit, gate_state,
                     following_gates_idx, qubit_state, gate_dependency,
                     updated_layer_gates_idx, future_layer, 1);
      trans_timer.stop_timer();
      total_maIdxTypeainlayer_time += trans_timer.measure();

      current_layer = updated_layer_gates_idx;
      executed_gates_num += execute_gates_idx.size();
      no_progress_steps = 0;
      if (enable_decay && !decay.empty()) {
        std::fill(decay.begin(), decay.end(), 1.0);
        decay_steps = 0;
      }
    } else {
      if (debug_level > 1 && (no_progress_steps % 200 == 0)) {
        cout << "[sabre] no executable gates; no_progress_steps="
             << no_progress_steps << " cur=" << current_layer.size()
             << " fut=" << future_layer.size() << endl;
      }
      if (debug_level > 1 && (no_progress_steps % 500 == 0)) {
        cout << "[sabre] stuck snapshot: current_layer gates=" << current_layer.size() << endl;
        for (IdxType gate_idx : current_layer) {
          if (gate_idx < 0 || gate_idx >= gate_num) {
            continue;
          }
          const auto &g = circuit_gate[gate_idx];
          IdxType lq = g.qubit;
          IdxType lc = g.ctrl;
          IdxType pq = (lq >= 0 && lq < static_cast<IdxType>(mapping.size())) ? mapping[lq] : -1;
          IdxType pc = (lc >= 0 && lc < static_cast<IdxType>(mapping.size())) ? mapping[lc] : -1;
          IdxType dist = -1;
          if (pq >= 0 && pc >= 0 &&
              pq < static_cast<IdxType>(distance_mat.size()) &&
              pc < static_cast<IdxType>(distance_mat[pq].size())) {
            dist = distance_mat[pq][pc];
          }
          cout << "  gate " << gate_idx << " (" << lc << "," << lq
               << ") phys=(" << pc << "," << pq << ") dist=" << dist << endl;
        }
      }
      cpu_timer trans_timer;
      trans_timer.start_timer();
      vector<IdxType> pair;
      if (attempt_limit > 0 && no_progress_steps >= attempt_limit) {
        pair = force_route_pair(mapping, current_layer, distance_mat, circuit,
                                qubit_num, chip, debug_level);
        no_progress_steps = 0;
      } else {
        pair =
            pick_one_movement(mapping, current_layer, future_layer, distance_mat,
                              qubit_num, circuit, chip, &decay, decay_increment,
                              decay_reset, decay_steps, enable_decay, debug_level,
                              critical_gate_set.empty() ? nullptr : &critical_gate_set,
                              &recent_swaps);
        no_progress_steps += 1;
      }
      trans_timer.stop_timer();
      total_pickone_time += trans_timer.measure();
      if (pair.size() < 2 || pair[0] < 0 || pair[1] < 0) {
        if (debug_level > 1) {
          cout << "No valid swap candidate found. Aborting further routing iterations." << endl;
        }
        executed_gates_num = gate_num;
        break;
      }
      // cout << "swap " << pair[0] << " " << pair[1] << endl;
      // all_gate_output.push_back({pair[0], pair[1]});
      Gate SWAPG = Gate(OP::SWAP, IdxType(pair[1]), IdxType(pair[0]));
      record_gate(SWAPG);
      // all_gate_type.push_back(1);
      swap_num += 1;
      auto ordered = std::minmax(pair[0], pair[1]);
      recent_swaps.emplace_back(ordered.first, ordered.second);
      if (recent_swaps.size() > recent_limit) {
        recent_swaps.pop_front();
      }
    }
    layer_index += 1;
    iter_count += 1;
    // executed_gates_num += gate_num;
  }
  single_gate_count = 0;
  for (IdxType i = 0; i < single_gate_info.size(); i++) {
    if (single_gate_info[i].ctrl == -1 &&
        strcmp(OP_NAMES[single_gate_info[i].op_name], "MA") != 0 &&
        visited_gate.find(i) == visited_gate.end()) {
      Gate cur_gate = single_gate_info[i];
      IdxType q_qubit = mapping[single_gate_info[i].qubit];
      cur_gate.qubit = q_qubit;
      record_gate(cur_gate);
    }
  }
  if (critical_path_indices != nullptr || critical_path_latency != nullptr) {
    std::vector<IdxType> path_indices;
    IdxType current = max_finish_index;
    while (current >= 0) {
      path_indices.push_back(current);
      current = gate_predecessor[current];
    }
    std::reverse(path_indices.begin(), path_indices.end());
    if (critical_path_indices != nullptr) {
      *critical_path_indices = path_indices;
    }
    if (critical_path_latency != nullptr) {
      *critical_path_latency = max_finish_time;
    }
  }
  initial_mapping = mapping;
  if (debug_level > 1) {
    cout << "total maIdxTypeainlayer time is: " << fixed << setprecision(1)
         << total_maIdxTypeainlayer_time << endl;
    cout << "total pick one swap time is: " << fixed << setprecision(1)
         << total_pickone_time << endl;
  }
  return swap_num;
}

void Routing(shared_ptr<Circuit> circuit, shared_ptr<Chip> chip,
             IdxType debug_level, bool enable_decay = false,
             double decay_increment = 0.001, IdxType decay_reset = 5,
             IdxType routing_trials = 1,
             std::optional<uint64_t> seed = std::nullopt) {
  IdxType n_qubits = IdxType(circuit->num_qubits());
  vector<Gate> gate_info = circuit->get_gates();

  vector<Gate> cx_gates;
  for (IdxType i = 0; i < gate_info.size(); i++) {
    Gate gate = gate_info[i];
    if (gate.ctrl != -1 && strcmp(OP_NAMES[gate.op_name], "MA") != 0) {
      cx_gates.push_back(gate);
    }
  }
  struct RoutingTrialResult {
    std::vector<Gate> gates;
    std::vector<IdxType> mapping;
    std::vector<IdxType> critical_path;
    double critical_path_latency = 0.0;
    IdxType swap_count = 0;
  };

  auto count_swaps = [](const std::vector<Gate> &gates) -> IdxType {
    IdxType count = 0;
    for (const auto &gate : gates) {
      if (gate.op_name == OP::SWAP) {
        count += 1;
      }
    }
    return count;
  };

  auto run_trial = [&](std::vector<IdxType> initial_mapping,
                       IdxType trial_debug) -> RoutingTrialResult {
    RoutingTrialResult result;
    std::vector<Gate> cx_gates_local = cx_gates;
    std::vector<vector<IdxType>> all_gates_index;
    std::vector<Gate> return_circuit;

    if (trial_debug > 1)
      cout << "******* 1st round sabre optimization *******" << endl;
    one_round_optimization(initial_mapping, cx_gates_local, chip->distance_mat,
                           gate_info, chip, all_gates_index, return_circuit,
                           trial_debug, nullptr, nullptr, enable_decay,
                           decay_increment, decay_reset);

    if (trial_debug > 1)
      cout << "******* 2nd round sabre optimization *******" << endl;
    std::reverse(cx_gates_local.begin(), cx_gates_local.end());
    std::vector<vector<IdxType>> reverse_gate_qubit;
    for (IdxType i = all_gates_index.size(); i > 0; --i) {
      reverse_gate_qubit.push_back(all_gates_index[i - 1]);
    }
    return_circuit.clear();
    one_round_optimization(initial_mapping, cx_gates_local, chip->distance_mat,
                           gate_info, chip, reverse_gate_qubit, return_circuit,
                           trial_debug, nullptr, nullptr, enable_decay,
                           decay_increment, decay_reset);

    if (trial_debug > 1)
      cout << "******* 3rd round sabre optimization *******" << endl;
    return_circuit.clear();
    std::reverse(cx_gates_local.begin(), cx_gates_local.end());
    if (trial_debug > 1) {
      cout << "initial mapping is:";
      for (IdxType ini : initial_mapping) {
        cout << ini << " ";
      }
      cout << endl;
    }
    std::vector<IdxType> critical_path_indices;
    double critical_path_latency = 0.0;
    one_round_optimization(
        initial_mapping, cx_gates_local, chip->distance_mat, gate_info, chip,
        all_gates_index, return_circuit, trial_debug, &critical_path_indices,
        &critical_path_latency, enable_decay, decay_increment, decay_reset);

    result.gates = return_circuit;
    result.mapping = initial_mapping;
    result.critical_path = critical_path_indices;
    result.critical_path_latency = critical_path_latency;
    result.swap_count = count_swaps(return_circuit);
    return result;
  };

  IdxType trials = routing_trials < 1 ? 1 : routing_trials;
  mt19937 g;
  if (seed.has_value()) {
    g.seed(static_cast<mt19937::result_type>(*seed));
  } else {
    random_device rd;
    g.seed(rd());
  }

  RoutingTrialResult best_result;
  bool have_best = false;

  const auto is_valid_mapping = [&](const std::vector<IdxType> &mapping) {
    if (mapping.size() < static_cast<std::size_t>(n_qubits)) {
      return false;
    }
    std::unordered_set<IdxType> used;
    used.reserve(mapping.size());
    for (IdxType i = 0; i < n_qubits; ++i) {
      IdxType phys = mapping[static_cast<std::size_t>(i)];
      if (phys < 0 || phys >= chip->chip_qubit_num) {
        return false;
      }
      if (!used.insert(phys).second) {
        return false;
      }
    }
    return true;
  };

  std::vector<IdxType> preset_mapping = circuit->get_mapping();
  if (is_valid_mapping(preset_mapping)) {
    auto result = run_trial(preset_mapping, debug_level);
    best_result = std::move(result);
    have_best = true;
  }
  IdxType layout_routing_trials = routing_trials < 1 ? 1 : routing_trials;
  if (layout_routing_trials > 3) {
    layout_routing_trials = 3;
  }

  for (IdxType trial = 0; trial < trials; ++trial) {
    std::vector<IdxType> initial_mapping(n_qubits, 0);
    iota(initial_mapping.begin(), initial_mapping.end(), 0);
    if (seed.has_value()) {
      g.seed(static_cast<mt19937::result_type>(*seed + trial));
    }
    shuffle(initial_mapping.begin(), initial_mapping.end(), g);

    IdxType trial_debug = (trials == 1) ? debug_level : 0;
    auto result = run_trial(initial_mapping, trial_debug);
    if (!have_best || result.swap_count < best_result.swap_count ||
        (result.swap_count == best_result.swap_count &&
         result.gates.size() < best_result.gates.size())) {
      best_result = std::move(result);
      have_best = true;
    }
  }

  if (have_best) {
    circuit->set_mapping(best_result.mapping);
    circuit->set_gates(best_result.gates);
    circuit->set_critical_path(best_result.critical_path,
                               best_result.critical_path_latency);
    circuit->set_routing_swap_count(best_result.swap_count);
  }
}

void SabreLayout(shared_ptr<Circuit> circuit, shared_ptr<Chip> chip,
                 IdxType debug_level, bool enable_decay = false,
                 double decay_increment = 0.001, IdxType decay_reset = 5,
                 IdxType routing_trials = 1,
                 std::optional<uint64_t> seed = std::nullopt,
                 IdxType layout_iterations = 2,
                 IdxType layout_trials = 1) {
  if (!circuit || !chip) {
    return;
  }
  IdxType kLayoutIterations = layout_iterations > 0 ? layout_iterations : 2;
  IdxType trials = layout_trials > 0 ? layout_trials : 1;
  std::vector<IdxType> best_mapping = circuit->get_mapping();
  IdxType best_swaps = std::numeric_limits<IdxType>::max();
  bool have_best = false;

  const std::vector<Gate> base_gates = circuit->get_gates();
  std::vector<Gate> reverse_gates = base_gates;
  std::reverse(reverse_gates.begin(), reverse_gates.end());

  IdxType layout_routing_trials = routing_trials < 1 ? 1 : routing_trials;
  if (layout_routing_trials > 3) {
    layout_routing_trials = 3;
  }

  std::vector<IdxType> dense_mapping;
  bool have_dense_mapping = false;
  {
    const IdxType logical_qubits = static_cast<IdxType>(circuit->num_qubits());
    const IdxType physical_qubits = chip ? chip->chip_qubit_num : 0;
    if (logical_qubits > 0 && physical_qubits >= logical_qubits) {
      std::vector<double> logical_weight(static_cast<std::size_t>(logical_qubits), 0.0);
      for (const auto &gate : base_gates) {
        if (gate.ctrl >= 0 && gate.qubit >= 0) {
          if (gate.ctrl < logical_qubits) {
            logical_weight[static_cast<std::size_t>(gate.ctrl)] += 1.0;
          }
          if (gate.qubit < logical_qubits) {
            logical_weight[static_cast<std::size_t>(gate.qubit)] += 1.0;
          }
        }
      }
      std::vector<std::pair<IdxType, double>> logical_rank;
      logical_rank.reserve(static_cast<std::size_t>(logical_qubits));
      for (IdxType q = 0; q < logical_qubits; ++q) {
        logical_rank.emplace_back(q, logical_weight[static_cast<std::size_t>(q)]);
      }
      std::sort(logical_rank.begin(), logical_rank.end(),
                [](const auto &a, const auto &b) { return a.second > b.second; });

      std::vector<std::pair<IdxType, std::pair<double, double>>> phys_rank;
      phys_rank.reserve(static_cast<std::size_t>(physical_qubits));
      for (IdxType p = 0; p < physical_qubits; ++p) {
        double degree = 0.0;
        if (p >= 0 && p < static_cast<IdxType>(chip->edge_list.size())) {
          degree = static_cast<double>(chip->edge_list[static_cast<std::size_t>(p)].size());
        }
        double dist_sum = 0.0;
        if (p >= 0 && p < static_cast<IdxType>(chip->distance_mat.size())) {
          for (IdxType d : chip->distance_mat[static_cast<std::size_t>(p)]) {
            if (d < std::numeric_limits<IdxType>::max()) {
              dist_sum += static_cast<double>(d);
            }
          }
        }
        phys_rank.emplace_back(p, std::make_pair(degree, -dist_sum));
      }
      std::sort(phys_rank.begin(), phys_rank.end(),
                [](const auto &a, const auto &b) {
                  if (a.second.first != b.second.first) {
                    return a.second.first > b.second.first;
                  }
                  return a.second.second > b.second.second;
                });

      dense_mapping.assign(static_cast<std::size_t>(logical_qubits), -1);
      for (IdxType i = 0; i < logical_qubits; ++i) {
        dense_mapping[static_cast<std::size_t>(logical_rank[static_cast<std::size_t>(i)].first)] =
            phys_rank[static_cast<std::size_t>(i)].first;
      }
      have_dense_mapping = true;
    }
  }

  for (IdxType trial = 0; trial < trials; ++trial) {
    std::vector<IdxType> trial_mapping = circuit->get_mapping();
    if (trial_mapping.empty()) {
      trial_mapping.resize(circuit->num_qubits(), 0);
      std::iota(trial_mapping.begin(), trial_mapping.end(), 0);
    }
    if (trial == 0 && have_dense_mapping) {
      trial_mapping = dense_mapping;
    } else {
      std::mt19937 g;
      if (seed.has_value()) {
        g.seed(static_cast<mt19937::result_type>(*seed + trial));
      } else {
        std::random_device rd;
        g.seed(rd());
      }
      std::shuffle(trial_mapping.begin(), trial_mapping.end(), g);
    }

    for (int iter = 0; iter < kLayoutIterations; ++iter) {
      std::optional<uint64_t> iter_seed = seed;
      if (seed.has_value()) {
        iter_seed = static_cast<uint64_t>(*seed + static_cast<uint64_t>(trial * kLayoutIterations + iter));
      }

      // Forward pass.
      auto forward_circuit = std::make_shared<Circuit>(circuit->num_qubits());
      forward_circuit->set_gates(base_gates);
      forward_circuit->set_mapping(trial_mapping);
      forward_circuit->set_creg(circuit->get_cregs());
      Routing(forward_circuit, chip, 0, enable_decay, decay_increment,
              decay_reset, layout_routing_trials, iter_seed);
      trial_mapping = forward_circuit->get_mapping();

      // Backward pass on reversed circuit.
      auto backward_circuit = std::make_shared<Circuit>(circuit->num_qubits());
      backward_circuit->set_gates(reverse_gates);
      backward_circuit->set_mapping(trial_mapping);
      backward_circuit->set_creg(circuit->get_cregs());
      Routing(backward_circuit, chip, 0, enable_decay, decay_increment,
              decay_reset, layout_routing_trials, iter_seed);
      trial_mapping = backward_circuit->get_mapping();
    }

    auto eval_circuit = std::make_shared<Circuit>(circuit->num_qubits());
    eval_circuit->set_gates(base_gates);
    eval_circuit->set_mapping(trial_mapping);
    eval_circuit->set_creg(circuit->get_cregs());
    Routing(eval_circuit, chip, 0, enable_decay, decay_increment,
            decay_reset, layout_routing_trials, seed);
    IdxType swap_count = eval_circuit->get_routing_swap_count();

    if (!have_best || swap_count < best_swaps) {
      best_swaps = swap_count;
      best_mapping = trial_mapping;
      have_best = true;
      if (best_swaps == 0) {
        break;
      }
    }
  }

  if (have_best) {
    circuit->set_mapping(best_mapping);
  }
  if (debug_level > 0) {
    std::cout << "STEP-2. Sabre layout stage selected initial mapping." << std::endl;
  }
}
