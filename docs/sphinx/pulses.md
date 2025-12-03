# Pulses and QICK emission

QASMTrans can translate a routed circuit into a timed pulse schedule and either write it to disk or stream it to QICK hardware. The flow has two layers:

- C++ CLI: add `-p` to dump a pulse JSON next to the output QASM; add `-e` (and optionally `--emit-run`) to push that JSON to a QICK board.
- Python bindings: use `qasmtrans_core.transpile_qasm` to get the pulse JSON in memory and call `emit_qick` to stream or summarize it without touching the filesystem.

## Inputs

- **Pulse template** (`-p` / `TranspileOptions.pulse_template_path`): JSON describing calibrated pulses. Required to produce any pulse output. QASMTrans matches each logical gate+qubits to a `pulse_definitions` entry (see `python/demo_qick.py::build_minimal_pulse_template` for a tiny example).
- **Backend config** (`-m`/`-c` on CLI, `backend_config` in bindings): the same device JSON used for transpilation; timing info is embedded in the pulse JSON header.
- **QICK config** (`-e` on CLI, `emit_qick` argument): maps logical qubits to generator channels and carries board settings such as frequencies and readout channels. At minimum supply `{"qubit_gen_map": {"0": 0, "1": 1}}`; optional keys include `pulse_freq`, `readout_freq`, `ro_chs`, `readout_gen_ch`, `pulse_gain_scale`, `default_gain`, `default_phase_deg`, `start_offset`, etc. (see `python/qasmtrans/qick_pulse_runner.py`).

## CLI workflow

```bash
# Transpile, generate QASM + pulses; save pulses as <output>_pulses.json
./qasmtrans -i data/test_benchmark/bv10.qasm \
  -m ibmq -c data/devices/ibmq_toronto.json \
  -p data/pulse_templates/toronto.json

# Also emit to QICK using a board config; omit --emit-run for summary-only
./qasmtrans -i data/test_benchmark/bv10.qasm \
  -m ibmq -c data/devices/ibmq_toronto.json \
  -p data/pulse_templates/toronto.json \
  -e my_qick_config.json --emit-run
```

Behavior:

- Pulse JSON is written to `data/output_qasm_file/<stem>_pulses.json` (or near `-o` if you set a custom output). The JSON bundles `backend` metadata, `pulse_library` (only entries that were used), `schedule` (timed events with parameters), and `critical_path` summaries.
- `-e` embeds a Python interpreter, prepends likely `python/` paths (or `QASMTRANS_PYTHON_PATH`), loads `python/qasmtrans/qick_pulse_runner.py`, and calls `run_program(pulse_path, qick_config, run_enabled, summary_only)`.
- Without `--emit-run`, emission is summary-only (gate histogram, generators used, total duration). With `--emit-run`, QICK hardware is programmed; missing `qick` Python package or device errors are surfaced on stderr.

## Python bindings workflow

```python
import json
import qasmtrans as qt  # pybind11 module: qasmtrans_core

opts = qt.TranspileOptions()
opts.backend_config = "data/devices/ibmq_toronto.json"
opts.pulse_template_path = "data/pulse_templates/toronto.json"
opts.output_path = "data/output_py"

result = qt.transpile_qasm("OPENQASM 2.0; ...", opts)
print("QASM saved to:", result.output_qasm_path)
print("Pulses saved to:", result.pulse_path)

# In-memory pulse document
pulse_doc = result.pulse_doc or json.loads(result.pulse_schedule)
qick_cfg = {"qubit_gen_map": {"0": 0, "1": 1}, "pulse_freq": 6.0e9}

# Summary only (no hardware)
summary = qt.emit_qick(pulse_doc, qick_cfg, run_enabled=False, summary_only=True)
print(summary)

# Attempt a hardware run (requires qick on a QICK board)
# qt.emit_qick(pulse_doc, qick_cfg, run_enabled=True)
```

Notes:

- `result.pulse_doc` is populated when the pulse JSON parses cleanly; otherwise fall back to `result.pulse_schedule` (raw string).
- `qt.load_qick_config(path)` loads a QICK config JSON into a Python dict; useful when driving everything from Python.
- The same runner module is reused by both C++ and Python flows, so behavior matches the CLI flags.

## What the runner does

- Parses the QASMTrans pulse JSON and builds events per qubit:
  - Virtual Z pulses become per-qubit phase accumulation.
  - Non-virtual pulses map qubits to generator channels via `qubit_gen_map`; waveform amplitude comes from samples (arbitrary/flat-top) or `amplitude`.
  - Events are sorted by `start_time`; the runner reports a gate histogram and total scheduled duration.
- If `run_enabled` and `qick` is installed:
  - Uploads arbitrary/flat-top waveforms (scaled to `gain_max`), sets pulse registers (freq/gain/phase), schedules pulses, optional readout, and runs an `AveragerProgram`.
  - Returns averaged I/Q data when hardware executes; otherwise returns a summary dict.

## Minimal example references

- `python/demo_qick.py`: builds a tiny pulse template and backend, runs `transpile_qasm`, plots the pulse schedule, and calls `emit_qick` in summary mode.
- `python/qasmtrans/qick_pulse_runner.py`: full runner implementation, including optional readout wiring and waveform preparation.

