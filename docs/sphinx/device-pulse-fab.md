# Device pulse fabrication utilities

Utilities were added to help build and test pulse libraries, especially for Rigetti Ankaa-class devices. These live under `device_pulse_fab/` and complement the pulse dumping/QICK emission flow.

## Components
- `device_pulse_fab/build_rigetti_ankaa9q_configs.py`: generates Rigetti Ankaa 9-qubit configurations (coupling, gate definitions) for use as backend/pulse templates.
- `device_pulse_fab/merge_pulse.py`: merges pulse definitions from multiple sources/templates into a consolidated library JSON.
- `device_pulse_fab/simulate_ankaa9q.py`: simulates pulse behavior or schedules on Ankaa-9 layouts; useful for dry-runs before hardware.

## How it fits the pipeline
- These scripts produce or manipulate the **pulse templates** consumed by `dumpPulses` (`include/dump_pulses.hpp`) and by QICK emission (`src/qick_emitter.cpp`, `python/qasmtrans/qick_pulse_runner.py`).
- Generated configs can be passed via CLI `-p` (pulse template) or via `TranspileOptions.pulse_template_path` in Python.
- Use them to keep template JSONs aligned with device calibrations before running the transpiler or emitting pulses.

## Review pointers
- Check the scripts’ argparse options for input/output paths and device parameters.
- Expect JSON outputs aligned with the template schema (`pulse_definitions` with `id/gate/qubits/shape/waveform_type/width/amplitude/...`).
