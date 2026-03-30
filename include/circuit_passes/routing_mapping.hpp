#pragma once

#include <algorithm>
#include <iomanip>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <string>
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

// extract cx gates in json file for constructing graph
vector<pair<IdxType, IdxType>> extract_cx_pairs(const json &j)
{
    vector<pair<IdxType, IdxType>> pairs;
    for (auto &item : j.items())
    {
        string key = item.key();
        if (key.substr(0, 2) == "cx")
        {
            IdxType pos = key.find('_');
            if (pos != string::npos)
            {
                IdxType first = stoi(key.substr(2, pos - 2));
                IdxType second = stoi(key.substr(pos + 1));
                pairs.push_back(make_pair(first, second));
            }
        }
    }
    return pairs;
}

void DAG_generator(IdxType qubit_num, vector<vector<IdxType>> &circuit, vector<IdxType> &gate_state, vector<IdxType> &qubit_state, vector<IdxType> &gate_dependency, vector<vector<IdxType>> &following_gate_idx, vector<IdxType> &first_layer_gates_idx)
{
    IdxType gate_num = circuit.size();
    vector<IdxType> current_gate_idx(qubit_num, -1);
    following_gate_idx.resize(gate_num, vector<IdxType>(2, 0));
    gate_dependency.resize(gate_num, 0);
    for (IdxType i = 0; i < gate_num; i++)
    {
        vector<IdxType> gate = circuit[i];
        if (current_gate_idx[gate[0]] == -1)
        {
            if (current_gate_idx[gate[1]] == -1)
            {
                first_layer_gates_idx.push_back(i);
                gate_state[i] = 2;
                qubit_state[gate[0]] = 1;
                qubit_state[gate[1]] = 1;
                gate_dependency[i] = 0;
            }
            else
            {
                gate_dependency[i] = 1;
            }
        }
        if (current_gate_idx[gate[1]] == -1 && current_gate_idx[gate[0]] != -1)
        {
            gate_dependency[i] = 1;
        }
        for (IdxType j = 0; j < gate.size(); j++)
        {
            IdxType qubit = gate[j];
            if (current_gate_idx[qubit] != -1)
            {
                vector<IdxType> prior_gate = circuit[current_gate_idx[qubit]];
                IdxType qubit_idx;
                if (prior_gate[j] != qubit)
                {
                    qubit_idx = 1 - j;
                }
                else
                {
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

void maintain_layer(vector<IdxType> &current_layer_gates_idx, set<IdxType> &gate_execute_idx_list, vector<vector<IdxType>> &circuit, vector<IdxType> &gate_state, vector<vector<IdxType>> &following_gate_idx, vector<IdxType> &qubit_state, vector<IdxType> &gate_dependency, vector<IdxType> &updated_layer_gates_idx, vector<IdxType> &future_layer_gates_idx, IdxType flag)
{
    unordered_set<IdxType> updated_set;
    updated_layer_gates_idx.clear();
    for (IdxType gate_idx : current_layer_gates_idx)
    {
        if (gate_execute_idx_list.count(gate_idx) > 0)
        {
            vector<IdxType> gate = circuit[gate_idx];
            gate_state[gate_idx] = 3;
            future_layer_gates_idx.erase(remove(future_layer_gates_idx.begin(), future_layer_gates_idx.end(), gate_idx), future_layer_gates_idx.end());

            qubit_state[gate[0]] = 0;
            qubit_state[gate[1]] = 0;
            vector<IdxType> following_gates = following_gate_idx[gate_idx];
            for (IdxType next_gate_idx : following_gates)
            {
                gate_dependency[next_gate_idx]--;
                if (gate_dependency[next_gate_idx] == 0)
                {
                    updated_set.insert(next_gate_idx);
                    gate_state[next_gate_idx] = 2;
                    future_layer_gates_idx.erase(remove(future_layer_gates_idx.begin(), future_layer_gates_idx.end(), next_gate_idx), future_layer_gates_idx.end());
                    qubit_state[circuit[next_gate_idx][0]] = 1;
                    qubit_state[circuit[next_gate_idx][1]] = 1;
                }
            }
        }
        else
        {
            updated_set.insert(gate_idx);
        }
    }
    updated_layer_gates_idx.assign(updated_set.begin(), updated_set.end());
    if (!updated_layer_gates_idx.empty())
    {
        IdxType start_gate = *min_element(updated_layer_gates_idx.begin(), updated_layer_gates_idx.end());
        for (IdxType gate_idx = start_gate; gate_idx < start_gate + 20 && gate_idx < circuit.size(); gate_idx++)
        {
            if (gate_state[gate_idx] == 0)
            {
                gate_state[gate_idx] = 1;
                if (flag != 0)
                {
                    future_layer_gates_idx.push_back(gate_idx);
                }
            }
        }
    }
    if (flag == 0)
    {
        for (IdxType gate_idx = 0; gate_idx < circuit.size(); gate_idx++)
        {
            if (gate_state[gate_idx] == 1)
            {
                future_layer_gates_idx.push_back(gate_idx);
            }
        }
    }
    sort(updated_layer_gates_idx.begin(), updated_layer_gates_idx.end());
    sort(future_layer_gates_idx.begin(), future_layer_gates_idx.end());
}

double heuristic(const vector<IdxType> &new_mapping, const vector<IdxType> &current_layer_gates_idx, const vector<IdxType> &future_gates_idx, const vector<vector<IdxType>> &distance_mat, const vector<vector<IdxType>> &circuit)
{
    double cost = 0.0;
    double first_cost = 0.0;
    if (current_layer_gates_idx.empty())
    {
        return 0;
    }
    for (IdxType gate_idx : current_layer_gates_idx)
    {
        vector<IdxType> gate = circuit[gate_idx];
        first_cost += distance_mat[new_mapping[gate[0]]][new_mapping[gate[1]]];
    }
    first_cost /= current_layer_gates_idx.size();
    if (future_gates_idx.empty())
    {
        cost = first_cost;
        return cost;
    }
    double second_cost = 0.0;
    for (IdxType gate_idx : future_gates_idx)
    {
        vector<IdxType> gate = circuit[gate_idx];
        second_cost += distance_mat[new_mapping[gate[0]]][new_mapping[gate[1]]];
    }
    second_cost /= future_gates_idx.size();
    cost = first_cost + 0.5 * second_cost;
    return cost;
}

vector<IdxType> find_reverse_mapping(const vector<IdxType> &mapping, IdxType qubit_num)
{
    vector<IdxType> reverse_mapping(qubit_num, -1);
    for (IdxType l_qubit = 0; l_qubit < mapping.size(); ++l_qubit)
    {
        IdxType p_qubit = mapping[l_qubit];
        if (p_qubit >= 0 && p_qubit < qubit_num)
        {
            reverse_mapping[p_qubit] = l_qubit;
        }
    }
    return reverse_mapping;
}

vector<IdxType> pick_one_movement(vector<IdxType> &mapping, const vector<IdxType> &current_layer, const vector<IdxType> &future_layer, const vector<vector<IdxType>> &distance_mat, IdxType qubit_num, const vector<vector<IdxType>> &circuit, shared_ptr<Chip> chip)
{
    vector<IdxType> l2p_mapping = mapping;
    vector<IdxType> key_p_qubits;
    for (IdxType gate_idx : current_layer)
    {
        vector<IdxType> gate = circuit[gate_idx];
        key_p_qubits.push_back(l2p_mapping[gate[0]]);
        key_p_qubits.push_back(l2p_mapping[gate[1]]);
    }
    vector<vector<IdxType>> possible_pairs;
    for (IdxType p_qubit : key_p_qubits)
    {
        for (IdxType p_qubit_target : chip->edge_list[p_qubit])
        {
            possible_pairs.push_back({p_qubit, p_qubit_target});
        }
    }
    vector<double> score(possible_pairs.size(), 0.0);
    for (size_t pair_idx = 0; pair_idx < possible_pairs.size(); ++pair_idx)
    {
        vector<IdxType> pair = possible_pairs[pair_idx];
        vector<IdxType> p2l_mapping = find_reverse_mapping(l2p_mapping, qubit_num);
        swap(p2l_mapping[pair[0]], p2l_mapping[pair[1]]);
        vector<IdxType> temp_l2p_mapping = find_reverse_mapping(p2l_mapping, qubit_num);
        score[pair_idx] = heuristic(temp_l2p_mapping, current_layer, future_layer, distance_mat, circuit);
    }
    size_t best_move_idx = distance(score.begin(), min_element(score.begin(), score.end()));
    vector<IdxType> pair = possible_pairs[best_move_idx];
    vector<IdxType> p2l_mapping = find_reverse_mapping(l2p_mapping, qubit_num);
    swap(p2l_mapping[pair[0]], p2l_mapping[pair[1]]);
    vector<IdxType> new_mapping = find_reverse_mapping(p2l_mapping, qubit_num);
    mapping = new_mapping;
    return pair;
}

set<IdxType> find_executable_gates(const vector<IdxType> &mapping, const vector<IdxType> &current_layer,
                                   const vector<vector<IdxType>> &circuit, const vector<vector<IdxType>> &distance_mat)
{
    set<IdxType> executable_gates;
    for (IdxType gate_idx : current_layer)
    {
        IdxType mapped_gate_zero = mapping[circuit[gate_idx][0]];
        IdxType mapped_gate_one = mapping[circuit[gate_idx][1]];
        if (distance_mat[mapped_gate_zero][mapped_gate_one] == 1)
        {
            executable_gates.insert(gate_idx);
        }
    }
    return executable_gates;
}

vector<pair<IdxType, IdxType>> sortWithSwaps(vector<IdxType> &lst)
{
    vector<IdxType> sorted_lst;
    vector<IdxType> temp_lst;
    for (IdxType x : lst)
    {
        if (x != -1)
        {
            temp_lst.push_back(x);
        }
    }
    for (IdxType x : lst)
    {
        if (x != -1)
        {
            sorted_lst.push_back(x);
        }
    }
    sort(sorted_lst.begin(), sorted_lst.end());
    vector<pair<IdxType, IdxType>> swaps;
    for (IdxType i = 0; i < lst.size(); i++)
    {
        if (lst[i] != -1 && lst[i] != sorted_lst[i])
        {
            IdxType j = find(lst.begin() + i, lst.end(), sorted_lst[i]) - lst.begin();
            swap(lst[i], lst[j]);
            swaps.push_back(make_pair(lst[i], lst[j]));
        }
    }
    return swaps;
}

IdxType one_round_optimization(vector<IdxType> &initial_mapping, vector<Gate> circuit_gate, vector<vector<IdxType>> distance_mat,
                               vector<Gate> gate_info, shared_ptr<Chip> chip, vector<vector<IdxType>> gate_qubit, vector<Gate> &return_circuit, IdxType debug_level)
{
    IdxType swap_num = 0;
    vector<IdxType> mapping = initial_mapping;

    IdxType executed_gates_num = 0;
    IdxType gate_num = circuit_gate.size();
    vector<vector<IdxType>> circuit(gate_num, vector<IdxType>(2, 0));
    for (IdxType i = 0; i < gate_num; i++)
    {
        circuit[i][0] = circuit_gate[i].ctrl;
        circuit[i][1] = circuit_gate[i].qubit;
    }
    IdxType qubit_num = distance_mat.size();
    vector<IdxType> gate_state(gate_num, 0);
    vector<IdxType> gate_dependency(gate_num, 2);
    vector<IdxType> qubit_state(qubit_num, 0);
    vector<vector<IdxType>> following_gates_idx;
    vector<IdxType> first_layer_gates_idx;
    DAG_generator(qubit_num, circuit, gate_state, qubit_state, gate_dependency, following_gates_idx, first_layer_gates_idx);
    vector<IdxType> current_layer;
    for (IdxType gate_idx : first_layer_gates_idx)
    {
        current_layer.push_back(gate_idx);
    }
    vector<IdxType> future_layer;
    set<IdxType> gate_execute_idx_list;
    vector<IdxType> updated_layer_gates_idx;
    maintain_layer(current_layer, gate_execute_idx_list, circuit, gate_state, following_gates_idx, qubit_state, gate_dependency, updated_layer_gates_idx, future_layer, 0);
    current_layer = updated_layer_gates_idx;
    IdxType single_gate_count = 0;
    vector<IdxType> num_single_before;
    vector<Gate> single_gate_info;
    for (IdxType i = 0; i < gate_info.size(); i++)
    {
        if (gate_info[i].ctrl == -1 && strcmp(OP_NAMES[gate_info[i].op_name], "MA") != 0)
        {
            single_gate_count++;
            single_gate_info.push_back(gate_info[i]);
        }
        else
        {
            num_single_before.push_back(single_gate_count);
        }
    }
    map<IdxType, vector<IdxType>> qubit_to_gate_indices;
    vector<vector<IdxType>> dependency_vector(gate_num);

    IdxType two_qubit_gate_index = 0;
    IdxType single_qubit_index = 0;

    for (IdxType i = 0; i < gate_info.size(); ++i)
    {
        const auto &gate = gate_info[i];
        if (gate.ctrl != -1)
        {
            if (qubit_to_gate_indices.count(gate.ctrl))
            {
                for (const auto &idx : qubit_to_gate_indices[gate.ctrl])
                {
                    dependency_vector[two_qubit_gate_index].push_back(idx);
                }
                qubit_to_gate_indices.erase(gate.ctrl);
            }
            if (qubit_to_gate_indices.count(gate.qubit))
            {
                for (const auto &idx : qubit_to_gate_indices[gate.qubit])
                {
                    dependency_vector[two_qubit_gate_index].push_back(idx);
                }
                qubit_to_gate_indices.erase(gate.qubit);
            }
            two_qubit_gate_index++;
        }
        else
        {
            if (strcmp(OP_NAMES[gate_info[i].op_name], "MA") != 0)
            {
                qubit_to_gate_indices[gate.qubit].push_back(single_qubit_index++);
            }
        }
    }
    double total_maintainlayer_time = 0;
    double total_pickone_time = 0;
    set<IdxType> visited_gate;
    while (executed_gates_num < gate_num)
    {
        set<IdxType> execute_gates_idx = find_executable_gates(mapping, current_layer, circuit, distance_mat);
        for (IdxType ee : execute_gates_idx)
        {
            vector<IdxType> cur_index_vector = dependency_vector[ee];
            for (IdxType cur_index : cur_index_vector)
            {
                Gate cur_gate = single_gate_info[cur_index];
                IdxType q_qubit = mapping[single_gate_info[cur_index].qubit];
                cur_gate.qubit = q_qubit;
                return_circuit.push_back(cur_gate);
                visited_gate.insert(cur_index);
            }
            Gate cur_gate = circuit_gate[ee];
            IdxType q_qubit = mapping[cur_gate.qubit];
            IdxType c_qubit = mapping[cur_gate.ctrl];
            cur_gate.qubit = q_qubit;
            cur_gate.ctrl = c_qubit;
            return_circuit.push_back(cur_gate);
        }
        if (!execute_gates_idx.empty())
        {
            cpu_timer trans_timer;
            trans_timer.start_timer();
            maintain_layer(current_layer, execute_gates_idx, circuit, gate_state, following_gates_idx, qubit_state, gate_dependency, updated_layer_gates_idx, future_layer, 1);
            trans_timer.stop_timer();
            total_maintainlayer_time += trans_timer.measure();

            current_layer = updated_layer_gates_idx;
            executed_gates_num += execute_gates_idx.size();
        }
        else
        {
            cpu_timer trans_timer;
            trans_timer.start_timer();
            vector<IdxType> pair = pick_one_movement(mapping, current_layer, future_layer, distance_mat, qubit_num, circuit, chip);
            trans_timer.stop_timer();
            total_pickone_time += trans_timer.measure();
            Gate SWAPG = Gate(OP::SWAP, IdxType(pair[1]), IdxType(pair[0]));
            return_circuit.push_back(SWAPG);
            swap_num += 1;
        }
    }
    single_gate_count = 0;
    for (IdxType i = 0; i < single_gate_info.size(); i++)
    {
        if (single_gate_info[i].ctrl == -1 && strcmp(OP_NAMES[single_gate_info[i].op_name], "MA") != 0 && visited_gate.find(i) == visited_gate.end())
        {
            Gate cur_gate = single_gate_info[i];
            IdxType q_qubit = mapping[single_gate_info[i].qubit];
            cur_gate.qubit = q_qubit;
            return_circuit.push_back(cur_gate);
        }
    }
    initial_mapping = mapping;
    if (debug_level > 1)
    {
        cout << "total maintainlayer time is: " << fixed << setprecision(1)
             << total_maintainlayer_time << endl;
        cout << "total pick one swap time is: " << fixed << setprecision(1) << total_pickone_time << endl;
    }
    return swap_num;
}

void Routing(shared_ptr<Circuit> circuit, shared_ptr<Chip> chip,
             IdxType debug_level, bool enable_decay = false,
             double decay_increment = 0.001, IdxType decay_reset = 5,
             std::optional<uint64_t> seed = std::nullopt)
{
    (void)enable_decay;
    (void)decay_increment;
    (void)decay_reset;

    IdxType n_qubits = IdxType(circuit->num_qubits());
    vector<Gate> gate_info = circuit->get_gates();

    vector<Gate> cx_gates;
    for (IdxType i = 0; i < gate_info.size(); i++)
    {
        Gate gate = gate_info[i];
        if (gate.ctrl != -1 && strcmp(OP_NAMES[gate.op_name], "MA") != 0)
        {
            cx_gates.push_back(gate);
        }
    }

    mt19937 g;
    if (seed.has_value())
    {
        g.seed(static_cast<mt19937::result_type>(*seed));
    }
    else
    {
        random_device rd;
        g.seed(rd());
    }

    vector<Gate> best_return_circuit;
    vector<IdxType> best_mapping;
    IdxType best_swap_num = -1;
    const vector<IdxType> requested_initial_mapping = circuit->get_mapping();
    bool has_explicit_initial_mapping = requested_initial_mapping.size() == static_cast<size_t>(n_qubits);
    if (has_explicit_initial_mapping)
    {
        std::vector<IdxType> seen(n_qubits, 0);
        for (IdxType physical : requested_initial_mapping)
        {
            if (physical < 0 || physical >= n_qubits || seen[physical] != 0)
            {
                has_explicit_initial_mapping = false;
                break;
            }
            seen[physical] = 1;
        }
    }

    vector<vector<IdxType>> all_gates_index;
    for (IdxType trial = 0; trial < 1; ++trial)
    {
        vector<IdxType> initial_mapping;
        if (has_explicit_initial_mapping)
        {
            initial_mapping = requested_initial_mapping;
        }
        else
        {
            initial_mapping.resize(n_qubits, 0);
            iota(initial_mapping.begin(), initial_mapping.end(), 0);
            if (seed.has_value())
            {
                g.seed(static_cast<mt19937::result_type>(*seed + static_cast<uint64_t>(trial)));
            }
            shuffle(initial_mapping.begin(), initial_mapping.end(), g);
        }

        if (has_explicit_initial_mapping)
        {
            if (debug_level > 1)
                cout << "******* fixed-layout routing *******" << endl;
            vector<Gate> return_circuit;
            IdxType swap_num = one_round_optimization(initial_mapping, cx_gates, chip->distance_mat, gate_info, chip, all_gates_index, return_circuit, debug_level);
            best_swap_num = swap_num;
            best_mapping = initial_mapping;
            best_return_circuit = return_circuit;
            break;
        }

        if (debug_level > 1)
            cout << "******* 1st round sabre optimization *******" << endl;
        vector<Gate> return_circuit;
        return_circuit.clear();
        IdxType swap_num = one_round_optimization(initial_mapping, cx_gates, chip->distance_mat, gate_info, chip, all_gates_index, return_circuit, debug_level);

        if (debug_level > 1)
            cout << "******* 2nd round sabre optimization *******" << endl;
        reverse(cx_gates.begin(), cx_gates.end());
        vector<vector<IdxType>> reverse_gate_qubit;
        for (IdxType i = all_gates_index.size(); i > 0; --i)
        {
            reverse_gate_qubit.push_back(all_gates_index[i - 1]);
        }
        return_circuit.clear();
        swap_num = one_round_optimization(initial_mapping, cx_gates, chip->distance_mat, gate_info, chip, reverse_gate_qubit, return_circuit, debug_level);

        if (debug_level > 1)
            cout << "******* 3rd round sabre optimization *******" << endl;
        return_circuit.clear();
        reverse(cx_gates.begin(), cx_gates.end());
        if (debug_level > 1)
        {
            cout << "initial mapping is:";
            for (IdxType ini : initial_mapping)
            {
                cout << ini << " ";
            }
            cout << endl;
        }
        swap_num = one_round_optimization(initial_mapping, cx_gates, chip->distance_mat, gate_info, chip, all_gates_index, return_circuit, debug_level);

        if (best_swap_num < 0 || swap_num < best_swap_num ||
            (swap_num == best_swap_num && return_circuit.size() < best_return_circuit.size()))
        {
            best_swap_num = swap_num;
            best_mapping = initial_mapping;
            best_return_circuit = return_circuit;
        }
    }

    circuit->clear_critical_path();
    circuit->set_mapping(best_mapping);
    circuit->set_routing_swap_count(best_swap_num < 0 ? 0 : best_swap_num);
    circuit->set_gates(best_return_circuit);
}
