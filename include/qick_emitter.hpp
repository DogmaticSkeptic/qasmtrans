#pragma once

#include <string>

namespace QASMTrans
{
    namespace pulses
    {
        void emit_with_qick(const std::string &pulse_path,
                             const std::string &qick_config_path,
                             bool run_enabled,
                             bool summary_only,
                             const std::string &executable_path,
                             bool verbose);
    } // namespace pulses
} // namespace QASMTrans

