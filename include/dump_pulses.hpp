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
            ValType width = 0.0;
            ValType amplitude = 0.0;
            std::string note;
        };

        struct PulseTemplateLibrary
        {
            std::string name;
            std::string version;
            std::map<std::string, PulseDefinition> definitions;
        };

        struct BackendTimingData
        {
            std::string name;
            std::string version;
            std::unordered_map<std::string, ValType> gate_lengths;
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

        inline std::string backendLensKey(const std::string &gate, const std::vector<IdxType> &qubits)
        {
            std::ostringstream oss;
            oss << gate;
            if (!qubits.empty())
            {
                oss << qubits[0];
                for (size_t i = 1; i < qubits.size(); ++i)
                {
                    oss << '_' << qubits[i];
                }
            }
            return oss.str();
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
                definition.width = entry.value("width", ValType{0.0});
                definition.amplitude = entry.value("amplitude", ValType{0.0});
                definition.note = entry.value("note", std::string{});
                definition.id = entry.value("id", std::string{});
                if (definition.gate.empty() || definition.qubits.empty())
                {
                    continue;
                }
                if (definition.id.empty())
                {
                    definition.id = makePulseIdentifier(definition.gate, definition.qubits);
                }
                const std::string key = makePulseKey(definition.gate, definition.qubits);
                library.definitions.insert({key, definition});
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

        inline PulseDefinition makeFallbackDefinition(const std::string &gate,
                                                      const std::vector<IdxType> &qubits,
                                                      const BackendTimingData &backend)
        {
            PulseDefinition fallback;
            fallback.gate = gate;
            fallback.qubits = qubits;
            fallback.shape = "gaussian";
            const std::string lens_key = backendLensKey(gate, qubits);
            auto length_it = backend.gate_lengths.find(lens_key);
            fallback.width = (length_it != backend.gate_lengths.end()) ? length_it->second : ValType{0.0};
            if (gate == "rz")
            {
                fallback.amplitude = 0.0;
            }
            else if (qubits.size() > 1)
            {
                fallback.amplitude = 0.2;
            }
            else
            {
                fallback.amplitude = 0.1;
            }
            fallback.note = "auto-generated fallback for missing pulse template";
            fallback.id = makePulseIdentifier(gate, qubits);
            return fallback;
        }
    } // namespace pulses

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

        std::map<std::string, PulseDefinition> definitions = library.definitions;
        std::set<std::string> used_keys;
        std::vector<json> schedule;
        std::unordered_map<IdxType, ValType> availability;
        std::vector<std::string> missing_definitions;

        const std::vector<Gate> gates = circuit->get_gates();
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
                PulseDefinition fallback = makeFallbackDefinition(gate_name, qubits, backend);
                definitions.insert({key, fallback});
                def_it = definitions.find(key);
                missing_definitions.push_back(key);
            }
            const PulseDefinition &definition = def_it->second;
            ValType start_time = 0.0;
            for (auto q : qubits)
            {
                auto avail_it = availability.find(q);
                if (avail_it != availability.end())
                {
                    start_time = std::max(start_time, avail_it->second);
                }
            }
            ValType duration = std::max(definition.width, ValType{0.0});
            ValType finish = start_time + duration;
            for (auto q : qubits)
            {
                availability[q] = finish;
            }

            json entry;
            entry["index"] = gate_index++;
            entry["gate"] = gate_name;
            entry["qubits"] = qubits;
            entry["pulse_id"] = definition.id;
            entry["start_time"] = start_time;
            entry["duration"] = duration;
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
            used_keys.insert(key);
        }

        json pulse_library = json::array();
        for (const auto &key : used_keys)
        {
            const auto &definition = definitions.at(key);
            json entry;
            entry["id"] = definition.id;
            entry["gate"] = definition.gate;
            entry["qubits"] = definition.qubits;
            if (!definition.shape.empty())
            {
                entry["shape"] = definition.shape;
            }
            entry["width"] = definition.width;
            entry["amplitude"] = definition.amplitude;
            if (!definition.note.empty())
            {
                entry["note"] = definition.note;
            }
            pulse_library.push_back(entry);
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
            total_duration = std::max(total_duration, pair.second);
        }
        backend_info["total_duration"] = total_duration;
        if (!missing_definitions.empty())
        {
            backend_info["auto_generated_pulses"] = missing_definitions;
        }

        json output;
        output["backend"] = backend_info;
        output["pulse_library"] = pulse_library;
        output["schedule"] = schedule;

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
