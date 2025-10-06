#pragma once

#include <algorithm>
#include <fstream>
#include <map>
#include <stdexcept>
#include <set>
#include <string>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>
#include <filesystem>
#include <cmath>
#include <limits>
#include <cstddef>

#include "QASMTransPrimitives.hpp"
#include "IR/gate.hpp"
#include "IR/circuit.hpp"
#include "nlomann/json.hpp"

namespace QASMTrans
{
    namespace pulses
    {
        using json = nlohmann::json;

        struct PulseDefinition
        {
            std::string id;
            std::string gate;
            std::vector<IdxType> qubits;
            std::string shape;
            std::string waveform_type;
            ValType width = 0.0;
            ValType amplitude = 0.0;
            std::string note;
            std::vector<ValType> samples_i;
            std::vector<ValType> samples_q;
            bool is_virtual = false;
            std::map<std::string, ValType> parameters;
        };

        struct PulseTemplateLibrary
        {
            std::string name;
            std::string version;
            std::map<std::string, std::vector<PulseDefinition>> definitions;
        };

        struct BackendTimingData
        {
            std::string name;
            std::string version;
            std::unordered_map<std::string, ValType> gate_lengths;
        };

        struct GateTimingInfo
        {
            std::string gate;
            std::vector<IdxType> qubits;
            ValType start = 0.0;
            ValType duration = 0.0;
            ValType finish = 0.0;
            std::ptrdiff_t predecessor = -1;
        };

        struct GateContribution
        {
            ValType total_duration = 0.0;
            size_t count = 0;
        };

        struct CriticalPathResult
        {
            ValType total_duration = 0.0;
            std::vector<size_t> gate_indices;
            std::unordered_map<std::string, GateContribution> contributions;
        };

        inline std::string toLower(std::string value)
        {
            std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
                           { return static_cast<char>(std::tolower(c)); });
            return value;
        }

        inline std::string joinQubits(const std::vector<IdxType> &qubits, const std::string &delimiter)
        {
            std::ostringstream oss;
            for (size_t i = 0; i < qubits.size(); ++i)
            {
                if (i > 0)
                {
                    oss << delimiter;
                }
                oss << qubits[i];
            }
            return oss.str();
        }

        inline std::string makePulseKey(const std::string &gate, const std::vector<IdxType> &qubits)
        {
            std::ostringstream oss;
            oss << gate << ":";
            oss << joinQubits(qubits, "_");
            return oss.str();
        }

        inline std::string makePulseIdentifier(const std::string &gate, const std::vector<IdxType> &qubits)
        {
            std::ostringstream oss;
            oss << gate;
            for (auto q : qubits)
            {
                oss << "_q" << q;
            }
            return oss.str();
        }

        inline std::vector<IdxType> extractGateQubits(const Gate &gate)
        {
            std::vector<IdxType> qubits;
            if (gate.ctrl >= 0)
            {
                qubits.push_back(gate.ctrl);
            }
            if (gate.qubit >= 0)
            {
                qubits.push_back(gate.qubit);
            }
            if (gate.extra >= 0)
            {
                qubits.push_back(gate.extra);
            }
            return qubits;
        }

        inline ValType gateParameterValue(const Gate &gate, const std::string &name)
        {
            if (name == "theta")
            {
                return gate.theta;
            }
            if (name == "phi")
            {
                return gate.phi;
            }
            if (name == "lambda" || name == "lam")
            {
                return gate.lam;
            }
            if (name == "gamma")
            {
                return gate.gamma;
            }
            return std::numeric_limits<ValType>::quiet_NaN();
        }

        inline const PulseDefinition *selectPulseDefinition(const std::vector<PulseDefinition> &candidates,
                                                           const Gate &gate,
                                                           const std::string &gate_name)
        {
            if (candidates.empty())
            {
                return nullptr;
            }
            if (candidates.size() == 1)
            {
                return &candidates.front();
            }
            constexpr ValType PARAM_TOL = 1e-6;
            const PulseDefinition *fallback = nullptr;
            for (const auto &candidate : candidates)
            {
                if (candidate.parameters.empty())
                {
                    if (fallback == nullptr)
                    {
                        fallback = &candidate;
                    }
                    continue;
                }
                bool matched = true;
                for (const auto &param : candidate.parameters)
                {
                    ValType gate_value = gateParameterValue(gate, param.first);
                    if (std::isnan(gate_value) || std::abs(gate_value - param.second) > PARAM_TOL)
                    {
                        matched = false;
                        break;
                    }
                }
                if (matched)
                {
                    return &candidate;
                }
            }
            if (fallback != nullptr)
            {
                return fallback;
            }
            return nullptr;
        }

        inline PulseTemplateLibrary loadPulseTemplate(const std::string &template_path)
        {
            PulseTemplateLibrary library;
            std::ifstream input(template_path);
            if (!input.is_open())
            {
                throw std::logic_error("Unable to open pulse template at " + template_path);
            }
            json document = json::parse(input, nullptr, true, true);
            library.name = document.value("name", std::string{});
            library.version = document.value("version", std::string{});
            if (!document.contains("pulse_definitions") || !document["pulse_definitions"].is_array())
            {
                throw std::logic_error("Pulse template " + template_path + " is missing 'pulse_definitions' array");
            }
            for (const auto &entry : document["pulse_definitions"])
            {
                PulseDefinition definition;
                definition.gate = toLower(entry.value("gate", std::string{}));
                definition.qubits = entry.value("qubits", std::vector<IdxType>{});
                definition.shape = entry.value("shape", std::string{});
                definition.waveform_type = entry.value("waveform_type", std::string{});
                definition.width = entry.value("width", ValType{0.0});
                definition.amplitude = entry.value("amplitude", ValType{0.0});
                definition.note = entry.value("note", std::string{});
                definition.id = entry.value("id", std::string{});
                definition.is_virtual = entry.value("virtual", false);
                if (entry.contains("parameters") && entry["parameters"].is_object())
                {
                    for (auto it = entry["parameters"].begin(); it != entry["parameters"].end(); ++it)
                    {
                        if (it.value().is_number())
                        {
                            definition.parameters.emplace(it.key(), it.value().get<ValType>());
                        }
                    }
                }
                if (entry.contains("samples_i") && entry["samples_i"].is_array())
                {
                    definition.samples_i = entry["samples_i"].get<std::vector<ValType>>();
                }
                if (entry.contains("samples_q") && entry["samples_q"].is_array())
                {
                    definition.samples_q = entry["samples_q"].get<std::vector<ValType>>();
                }
                if (definition.gate.empty() || definition.qubits.empty())
                {
                    continue;
                }
                if (definition.id.empty())
                {
                    definition.id = makePulseIdentifier(definition.gate, definition.qubits);
                }
                if (definition.waveform_type.empty())
                {
                    definition.waveform_type = definition.shape;
                }
                const std::string key = makePulseKey(definition.gate, definition.qubits);
                library.definitions[key].push_back(definition);
            }
            return library;
        }

        inline BackendTimingData loadBackendTiming(const std::string &backend_path)
        {
            BackendTimingData data;
            std::ifstream input(backend_path);
            if (!input.is_open())
            {
                throw std::logic_error("Unable to open backend config at " + backend_path);
            }
            json document = json::parse(input, nullptr, true, true);
            data.name = document.value("name", std::string{});
            data.version = document.value("version", std::string{});
            if (document.contains("gate_lens") && document["gate_lens"].is_object())
            {
                for (auto it = document["gate_lens"].begin(); it != document["gate_lens"].end(); ++it)
                {
                    data.gate_lengths.emplace(it.key(), it.value().get<ValType>());
                }
            }
            return data;
        }

    } // namespace pulses

    inline pulses::CriticalPathResult computeCriticalPath(const std::vector<pulses::GateTimingInfo> &timings)
    {
        using namespace pulses;
        CriticalPathResult result;
        if (timings.empty())
        {
            return result;
        }

        ValType max_finish = 0.0;
        std::ptrdiff_t max_index = -1;
        for (size_t idx = 0; idx < timings.size(); ++idx)
        {
            if (timings[idx].finish >= max_finish)
            {
                max_finish = timings[idx].finish;
                max_index = static_cast<std::ptrdiff_t>(idx);
            }
        }

        result.total_duration = max_finish;
        std::vector<size_t> path_indices;
        while (max_index >= 0)
        {
            size_t current = static_cast<size_t>(max_index);
            path_indices.push_back(current);
            max_index = timings[current].predecessor;
        }
        std::reverse(path_indices.begin(), path_indices.end());
        result.gate_indices = path_indices;

        for (size_t idx : path_indices)
        {
            const GateTimingInfo &info = timings[idx];
            GateContribution &entry = result.contributions[info.gate];
            entry.total_duration += info.duration;
            entry.count += 1;
        }

        return result;
    }

    inline void dumpPulses(std::shared_ptr<Circuit> circuit,
                           const char *input_filename,
                           const std::string &backend_path,
                           const std::string &pulse_template_path,
                           const std::string &output_path,
                           IdxType debug_level)
    {
        using namespace pulses;
        PulseTemplateLibrary library = loadPulseTemplate(pulse_template_path);
        BackendTimingData backend = loadBackendTiming(backend_path);

        const std::map<std::string, std::vector<PulseDefinition>> &definitions = library.definitions;
        std::set<std::string> used_ids;
        std::vector<json> schedule;
        struct QubitAvailability
        {
            ValType finish = 0.0;
            std::ptrdiff_t gate_index = -1;
        };
        std::unordered_map<IdxType, QubitAvailability> availability;
        const std::vector<Gate> gates = circuit->get_gates();
        std::vector<GateTimingInfo> gate_timings;
        gate_timings.reserve(gates.size());
        size_t gate_index = 0;
        for (const auto &gate : gates)
        {
            std::string gate_name = toLower(OP_NAMES[gate.op_name]);
            if (gate_name.empty())
            {
                continue;
            }
            if (gate_name == "ma" || gate_name == "m")
            {
                continue;
            }
            std::vector<IdxType> qubits = extractGateQubits(gate);
            if (qubits.empty())
            {
                continue;
            }
            const std::string key = makePulseKey(gate_name, qubits);
            auto def_it = definitions.find(key);
            if (def_it == definitions.end())
            {
                throw std::logic_error(
                    "Pulse template missing definition for gate '" + gate_name + "' on qubits [" +
                    joinQubits(qubits, ",") + "]");
            }
            const PulseDefinition *definition_ptr = selectPulseDefinition(def_it->second, gate, gate_name);
            if (definition_ptr == nullptr)
            {
                std::ostringstream oss;
                oss << "Pulse template missing matching definition for gate '" << gate_name
                    << "' on qubits [" << joinQubits(qubits, ",") << "]";
                if (gate.theta != 0.0)
                {
                    oss << " (theta=" << gate.theta << ")";
                }
                if (gate.phi != 0.0)
                {
                    oss << " (phi=" << gate.phi << ")";
                }
                if (gate.lam != 0.0)
                {
                    oss << " (lambda=" << gate.lam << ")";
                }
                if (gate.gamma != 0.0)
                {
                    oss << " (gamma=" << gate.gamma << ")";
                }
                throw std::logic_error(oss.str());
            }
            const PulseDefinition &definition = *definition_ptr;
            if (!definition.is_virtual && toLower(definition.waveform_type) == "arbitrary")
            {
                if (definition.samples_i.empty() && definition.samples_q.empty())
                {
                    throw std::logic_error("Arbitrary waveform for gate '" + gate_name + "' on qubits [" +
                                           joinQubits(qubits, ",") + "] is missing samples");
                }
            }
            ValType start_time = 0.0;
            std::ptrdiff_t predecessor = -1;
            for (auto q : qubits)
            {
                auto avail_it = availability.find(q);
                if (avail_it != availability.end())
                {
                    if (avail_it->second.finish >= start_time)
                    {
                        start_time = avail_it->second.finish;
                        predecessor = avail_it->second.gate_index;
                    }
                }
            }
            ValType duration = std::max(definition.width, ValType{0.0});
            ValType finish = start_time + duration;

            size_t current_gate_index = gate_index;

            json entry;
            entry["index"] = gate_index++;
            entry["gate"] = gate_name;
            entry["qubits"] = qubits;
            entry["pulse_id"] = definition.id;
            entry["start_time"] = start_time;
            entry["duration"] = duration;
            if (definition.is_virtual)
            {
                entry["virtual"] = true;
            }
            json parameters = json::object();
            if (gate.theta != 0.0)
            {
                parameters["theta"] = gate.theta;
            }
            if (gate.phi != 0.0)
            {
                parameters["phi"] = gate.phi;
            }
            if (gate.lam != 0.0)
            {
                parameters["lambda"] = gate.lam;
            }
            if (gate.gamma != 0.0)
            {
                parameters["gamma"] = gate.gamma;
            }
            if (!parameters.empty())
            {
                entry["parameters"] = parameters;
            }
            schedule.push_back(entry);
            for (auto q : qubits)
            {
                availability[q] = QubitAvailability{finish, static_cast<std::ptrdiff_t>(current_gate_index)};
            }
            used_ids.insert(definition.id);
            gate_timings.push_back(GateTimingInfo{gate_name, qubits, start_time, duration, finish, predecessor});
        }

        json pulse_library = json::array();
        for (const auto &pair : definitions)
        {
            for (const auto &definition : pair.second)
            {
                if (used_ids.find(definition.id) == used_ids.end())
                {
                    continue;
                }
                json entry;
                entry["id"] = definition.id;
                entry["gate"] = definition.gate;
                entry["qubits"] = definition.qubits;
                if (!definition.shape.empty())
                {
                    entry["shape"] = definition.shape;
                }
                if (!definition.waveform_type.empty())
                {
                    entry["waveform_type"] = definition.waveform_type;
                }
                entry["width"] = definition.width;
                entry["amplitude"] = definition.amplitude;
                if (!definition.note.empty())
                {
                    entry["note"] = definition.note;
                }
                if (definition.is_virtual)
                {
                    entry["virtual"] = true;
                }
                if (!definition.parameters.empty())
                {
                    json param_obj = json::object();
                    for (const auto &param : definition.parameters)
                    {
                        param_obj[param.first] = param.second;
                    }
                    entry["parameters"] = param_obj;
                }
                if (!definition.samples_i.empty())
                {
                    entry["samples_i"] = definition.samples_i;
                }
                if (!definition.samples_q.empty())
                {
                    entry["samples_q"] = definition.samples_q;
                }
                pulse_library.push_back(entry);
            }
        }

        json backend_info = json::object();
        backend_info["name"] = !backend.name.empty() ? backend.name : library.name;
        backend_info["version"] = !backend.version.empty() ? backend.version : library.version;
        backend_info["input_qasm"] = std::filesystem::path(input_filename).filename().string();
        backend_info["backend_config"] = backend_path;
        backend_info["pulse_template"] = pulse_template_path;
        backend_info["num_qubits"] = circuit->num_qubits();
        backend_info["total_pulses"] = schedule.size();
        ValType total_duration = 0.0;
        for (const auto &pair : availability)
        {
            total_duration = std::max(total_duration, pair.second.finish);
        }
        backend_info["total_duration"] = total_duration;
        json output;
        output["backend"] = backend_info;
        output["pulse_library"] = pulse_library;
        output["schedule"] = schedule;

        CriticalPathResult critical_path_result = computeCriticalPath(gate_timings);
        json critical_summary = json::object();
        critical_summary["total_duration"] = critical_path_result.total_duration;
        critical_summary["gate_count"] = critical_path_result.gate_indices.size();
        json path_indices = json::array();
        for (size_t idx : critical_path_result.gate_indices)
        {
            path_indices.push_back(idx);
        }
        critical_summary["path_indices"] = path_indices;

        std::vector<std::pair<std::string, GateContribution>> contributions(
            critical_path_result.contributions.begin(), critical_path_result.contributions.end());
        std::sort(contributions.begin(), contributions.end(), [](const auto &lhs, const auto &rhs)
                  {
                      if (std::abs(lhs.second.total_duration - rhs.second.total_duration) > 1e-9)
                      {
                          return lhs.second.total_duration > rhs.second.total_duration;
                      }
                      if (lhs.second.count != rhs.second.count)
                      {
                          return lhs.second.count > rhs.second.count;
                      }
                      return lhs.first < rhs.first;
                  });
        json top_gates = json::array();
        size_t gate_limit = std::min<size_t>(10, contributions.size());
        for (size_t i = 0; i < gate_limit; ++i)
        {
            json entry = json::object();
            entry["gate"] = contributions[i].first;
            entry["count"] = contributions[i].second.count;
            entry["total_duration"] = contributions[i].second.total_duration;
            top_gates.push_back(entry);
        }
        critical_summary["top_gates"] = top_gates;
        output["critical_path"] = critical_summary;

        std::filesystem::path out_path(output_path);
        if (out_path.has_parent_path())
        {
            std::filesystem::create_directories(out_path.parent_path());
        }
        std::ofstream out_file(out_path);
        if (!out_file.is_open())
        {
            throw std::logic_error("Unable to open output pulse file " + output_path);
        }
        out_file << output.dump(2) << std::endl;
        out_file.close();

        if (debug_level > 0)
        {
            std::cout << "Saved pulse schedule with " << schedule.size() << " entries to " << output_path << std::endl;
        }
    }
}
