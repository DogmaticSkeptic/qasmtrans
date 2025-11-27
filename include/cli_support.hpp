#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include <unordered_set>
#include <unordered_map>

#include "QASMTransPrimitives.hpp"
#include "IR/chip.hpp"
#include "IR/gate.hpp"
#include "nlomann/json.hpp"

namespace QASMTrans
{
    extern std::unordered_set<std::string> g_device_basis_gates;
    extern std::unordered_map<std::string, std::string> g_merged_gate_aliases;

    namespace cli
    {
        // Shared helpers for CLI orchestration
        std::string derive_pulse_output_path(const std::string &qasm_output_path);

        struct GateSummary
        {
            std::size_t single_qubit = 0;
            std::size_t two_qubit = 0;
            std::size_t depth = 0;
        };

        GateSummary compute_gate_summary(const std::vector<Gate> &gates, IdxType initial_capacity);

        struct PrunedSubchipData
        {
            std::vector<IdxType> new_to_old;
            std::vector<IdxType> old_to_new;
            std::vector<QASMTrans::Gate> gates;
            std::vector<IdxType> logical_mapping;
            std::vector<IdxType> measurement_mapping;
            std::vector<IdxType> local_to_global;
            nlohmann::json device_json;
        };

        void ingest_backend_metadata(const std::string &backendpath,
                                     std::unordered_set<std::string> &device_basis_gates,
                                     std::unordered_map<std::string, std::string> &merged_gate_aliases);

        PrunedSubchipData prune_subchip_artifact(const std::shared_ptr<Chip> &subchip,
                                                 const std::vector<IdxType> &local_to_global,
                                                 const std::vector<QASMTrans::Gate> &local_gates,
                                                 const std::vector<IdxType> &logical_mapping,
                                                 const std::vector<IdxType> &measurement_mapping,
                                                 const std::string &subchip_name,
                                                 IdxType mode);

        void write_json_file(const nlohmann::json &data, const std::filesystem::path &output_path);
    } // namespace cli
} // namespace QASMTrans
