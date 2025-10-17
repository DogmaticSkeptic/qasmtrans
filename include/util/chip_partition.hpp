#pragma once

#include <algorithm>
#include <functional>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../QASMTransPrimitives.hpp"
#include "../IR/chip.hpp"

namespace QASMTrans
{

inline std::vector<std::vector<IdxType>> partition_chip(const std::shared_ptr<Chip> &chip,
                                                        const std::vector<IdxType> &partition_sizes)
{
    if (!chip)
    {
        throw std::invalid_argument("partition_chip: null chip pointer");
    }

    if (partition_sizes.empty())
    {
        return {};
    }

    std::vector<IdxType> available_nodes(chip->chip_qubit_num);
    std::iota(available_nodes.begin(), available_nodes.end(), 0);

    IdxType total_requested = 0;
    for (auto size : partition_sizes)
    {
        if (size <= 0)
        {
            throw std::invalid_argument("partition_chip: partition sizes must be positive");
        }
        total_requested += size;
    }
    if (total_requested > static_cast<IdxType>(available_nodes.size()))
    {
        throw std::invalid_argument("partition_chip: requested qubits exceed device capacity");
    }

    std::vector<IdxType> degrees(available_nodes.size(), 0);
    for (IdxType node = 0; node < static_cast<IdxType>(chip->edge_list.size()); ++node)
    {
        degrees[node] = static_cast<IdxType>(chip->edge_list[node].size());
    }

    struct PartitionRequest
    {
        IdxType size;
        std::size_t original_index;
    };

    std::vector<PartitionRequest> requests;
    requests.reserve(partition_sizes.size());
    for (std::size_t i = 0; i < partition_sizes.size(); ++i)
    {
        requests.push_back(PartitionRequest{partition_sizes[i], i});
    }

    std::sort(requests.begin(), requests.end(),
              [](const PartitionRequest &a, const PartitionRequest &b)
              {
                  if (a.size != b.size)
                  {
                      return a.size > b.size;
                  }
                  return a.original_index < b.original_index;
              });

    std::vector<std::vector<IdxType>> result(partition_sizes.size());
    std::vector<char> used(available_nodes.size(), 0);

    auto bfs_collect = [&](IdxType seed, IdxType needed, const std::vector<char> &used_state) -> std::vector<IdxType>
    {
        if (used_state[seed])
        {
            return {};
        }

        std::vector<char> visited(available_nodes.size(), 0);
        std::queue<IdxType> q;
        std::vector<IdxType> collected;
        q.push(seed);
        visited[seed] = 1;
        collected.push_back(seed);

        while (!q.empty() && static_cast<IdxType>(collected.size()) < needed)
        {
            IdxType current = q.front();
            q.pop();
            if (current < 0 || current >= static_cast<IdxType>(chip->edge_list.size()))
            {
                continue;
            }

            std::vector<IdxType> neighbours;
            neighbours.reserve(chip->edge_list[current].size());
            for (IdxType neighbour : chip->edge_list[current])
            {
                if (neighbour < 0 || neighbour >= static_cast<IdxType>(available_nodes.size()))
                {
                    continue;
                }
                if (used_state[neighbour] || visited[neighbour])
                {
                    continue;
                }
                neighbours.push_back(neighbour);
            }
            std::sort(neighbours.begin(), neighbours.end(),
                      [&](IdxType lhs, IdxType rhs)
                      {
                          if (degrees[lhs] != degrees[rhs])
                          {
                              return degrees[lhs] < degrees[rhs];
                          }
                          return lhs < rhs;
                      });

            for (IdxType neighbour : neighbours)
            {
                if (static_cast<IdxType>(collected.size()) >= needed)
                {
                    break;
                }
                visited[neighbour] = 1;
                q.push(neighbour);
                collected.push_back(neighbour);
            }
        }

        if (static_cast<IdxType>(collected.size()) == needed)
        {
            return collected;
        }
        return {};
    };

    std::function<bool(std::size_t)> assign_partition = [&](std::size_t request_index) -> bool
    {
        if (request_index >= requests.size())
        {
            return true;
        }

        const auto &req = requests[request_index];
        std::vector<IdxType> candidate_seeds;
        candidate_seeds.reserve(available_nodes.size());
        for (IdxType node = 0; node < static_cast<IdxType>(available_nodes.size()); ++node)
        {
            if (!used[node])
            {
                candidate_seeds.push_back(node);
            }
        }
        std::sort(candidate_seeds.begin(), candidate_seeds.end(),
                  [&](IdxType lhs, IdxType rhs)
                  {
                      if (degrees[lhs] != degrees[rhs])
                      {
                          return degrees[lhs] < degrees[rhs];
                      }
                      return lhs < rhs;
                  });

        for (IdxType seed : candidate_seeds)
        {
            auto candidate_nodes = bfs_collect(seed, req.size, used);
            if (candidate_nodes.empty())
            {
                continue;
            }

            for (IdxType node : candidate_nodes)
            {
                used[node] = 1;
            }

            result[req.original_index] = candidate_nodes;

            if (assign_partition(request_index + 1))
            {
                return true;
            }

            for (IdxType node : candidate_nodes)
            {
                used[node] = 0;
            }
        }

        return false;
    };

    if (!assign_partition(0))
    {
        throw std::runtime_error("partition_chip: unable to allocate contiguous region for circuit");
    }

    return result;
}

inline std::shared_ptr<Chip> make_subchip(const std::shared_ptr<Chip> &chip,
                                          const std::vector<IdxType> &nodes,
                                          std::vector<IdxType> &local_to_global)
{
    if (!chip)
    {
        throw std::invalid_argument("make_subchip: null chip pointer");
    }
    if (nodes.empty())
    {
        throw std::invalid_argument("make_subchip: empty node set");
    }

    const std::size_t local_n = nodes.size();
    local_to_global = nodes;

    std::unordered_map<IdxType, IdxType> global_to_local;
    global_to_local.reserve(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        global_to_local[nodes[i]] = static_cast<IdxType>(i);
    }

    std::vector<std::vector<IdxType>> adj(local_n, std::vector<IdxType>(local_n, 0));
    std::vector<std::vector<IdxType>> edge_list(local_n);
    for (std::size_t i = 0; i < nodes.size(); ++i)
    {
        IdxType global_src = nodes[i];
        if (global_src < 0 || global_src >= static_cast<IdxType>(chip->edge_list.size()))
        {
            continue;
        }
        for (IdxType neighbour : chip->edge_list[global_src])
        {
            auto it = global_to_local.find(neighbour);
            if (it == global_to_local.end())
            {
                continue;
            }
            IdxType local_dst = it->second;
            adj[i][local_dst] = 1;
        }
    }
    for (std::size_t i = 0; i < local_n; ++i)
    {
        for (std::size_t j = 0; j < local_n; ++j)
        {
            if (adj[i][j] == 1)
            {
                edge_list[i].push_back(static_cast<IdxType>(j));
            }
        }
    }

    std::vector<std::vector<IdxType>> distance_mat = floyd(static_cast<IdxType>(local_n), adj);
    auto subchip = std::make_shared<Chip>(static_cast<IdxType>(local_n), adj, edge_list, distance_mat);
    subchip->chip_qubit_num = static_cast<IdxType>(local_n);

    subchip->single_qubit_errors.assign(local_n, {});
    subchip->single_qubit_gate_lengths.assign(local_n, {});
    for (std::size_t i = 0; i < local_n; ++i)
    {
        IdxType global = nodes[i];
        if (global < 0 || global >= static_cast<IdxType>(chip->single_qubit_errors.size()))
        {
            continue;
        }
        subchip->single_qubit_errors[i] = chip->single_qubit_errors[global];
        subchip->single_qubit_gate_lengths[i] = chip->single_qubit_gate_lengths[global];
    }

    for (const auto &entry : chip->two_qubit_errors)
    {
        auto it_ctrl = global_to_local.find(entry.first.first);
        auto it_tgt = global_to_local.find(entry.first.second);
        if (it_ctrl == global_to_local.end() || it_tgt == global_to_local.end())
        {
            continue;
        }
        subchip->two_qubit_errors[{it_ctrl->second, it_tgt->second}] = entry.second;
    }
    for (const auto &entry : chip->two_qubit_gate_lengths)
    {
        auto it_ctrl = global_to_local.find(entry.first.first);
        auto it_tgt = global_to_local.find(entry.first.second);
        if (it_ctrl == global_to_local.end() || it_tgt == global_to_local.end())
        {
            continue;
        }
        subchip->two_qubit_gate_lengths[{it_ctrl->second, it_tgt->second}] = entry.second;
    }

    return subchip;
}

} // namespace QASMTrans
