#include "../include/qick_emitter.hpp"

#include <pybind11/embed.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <mutex>
#include <vector>

namespace QASMTrans
{
    namespace pulses
    {
        namespace
        {
            namespace fs = std::filesystem;
            using GuardPtr = std::unique_ptr<pybind11::scoped_interpreter>;

            GuardPtr &interpreter_guard()
            {
                static GuardPtr guard;
                static std::mutex guard_mutex;
                std::lock_guard<std::mutex> lock(guard_mutex);
                if (!guard)
                {
                    guard = std::make_unique<pybind11::scoped_interpreter>();
                }
                return guard;
            }

            std::vector<fs::path> candidate_python_paths(const std::string &executable_path)
            {
                std::vector<fs::path> candidates;
                const fs::path exe_path = fs::absolute(fs::path(executable_path));
                const fs::path exe_dir = exe_path.parent_path();
                const fs::path current_python = fs::current_path() / "python";
                const fs::path sibling_python = exe_dir / "python";
                const fs::path parent_python = exe_dir.parent_path() / "python";
                const char *env_path = std::getenv("QASMTRANS_PYTHON_PATH");

                if (!current_python.empty())
                {
                    candidates.push_back(current_python);
                }
                if (!sibling_python.empty() && sibling_python != current_python)
                {
                    candidates.push_back(sibling_python);
                }
                if (!parent_python.empty() && parent_python != sibling_python)
                {
                    candidates.push_back(parent_python);
                }
                if (env_path != nullptr)
                {
                    fs::path env_candidate(env_path);
                    if (!env_candidate.empty())
                    {
                        candidates.push_back(fs::absolute(env_candidate));
                    }
                }
                return candidates;
            }

            void extend_sys_path(const std::vector<fs::path> &paths, bool verbose)
            {
                namespace py = pybind11;
                py::module sys = py::module::import("sys");
                py::list sys_path = sys.attr("path");
                for (const auto &candidate : paths)
                {
                    std::error_code ec;
                    if (!fs::exists(candidate, ec) || !fs::is_directory(candidate, ec))
                    {
                        continue;
                    }
                    std::string path_str = candidate.lexically_normal().string();
                    bool already_present = false;
                    for (const auto &entry : sys_path)
                    {
                        if (py::cast<std::string>(entry) == path_str)
                        {
                            already_present = true;
                            break;
                        }
                    }
                    if (!already_present)
                    {
                        sys_path.attr("insert")(0, path_str);
                        if (verbose)
                        {
                            std::cout << "[qick-emitter] Added Python path: " << path_str << std::endl;
                        }
                    }
                }
            }
        } // namespace

        void emit_with_qick(const std::string &pulse_path,
                            const std::string &qick_config_path,
                            bool run_enabled,
                            bool summary_only,
                            const std::string &executable_path,
                            bool verbose)
        {
            namespace py = pybind11;
            auto &guard = interpreter_guard();
            (void)guard;

            extend_sys_path(candidate_python_paths(executable_path), verbose);

            try
            {
                py::module runner = py::module::import("python.qick_pulse_runner");
                runner.attr("run_program")(pulse_path, qick_config_path, run_enabled, summary_only);
            }
            catch (const py::error_already_set &ex)
            {
                std::cerr << "Failed to emit QICK pulses: " << ex.what() << std::endl;
                throw;
            }
        }
    } // namespace pulses
} // namespace QASMTrans
