# Python package & bindings

The Python package (`qasmtrans`) wraps the pybind11 core (`qasmtrans_core`) and exposes the transpiler, pulse dumping, and QICK helpers.

## Package layout
- Core extension: built via CMake/pybind11 (`qasmtrans_core` .so).
- Shim module: `python/qasmtrans/__init__.py` tries to import `qasmtrans_core` from site-packages, bundled wheel layouts, or packaged resources.
- Demos/helpers: `python/demo.py`, `python/demo_qick.py`, plotting helper `python/qasmtrans/plot_pulses.py`, QICK runner `python/qasmtrans/qick_pulse_runner.py`.
- Packaging: `pyproject.toml`, `requirements*.txt` (dev vs minimal vs py38).

## Public API (from `src/python_binding.cpp`)

### TranspileOptions
Fields control transpilation, outputs, and pulse/QICK hooks:
- `mode`: basis set/mode (matches CLI `-m`).
- `backend_config`: path to device JSON (CLI `-c`).
- `pulse_template_path`: pulse template JSON; enables pulse dumping.
- `output_path`: base path for emitted QASM (CLI `-o`).
- `pulse_output_path`: override pulse JSON path (defaults to `<stem>_pulses.json` near `output_path`).
- `emit_run`: placeholder for future CLI parity (CLI uses `-e/--emit-run`).
- `limited_qubits`: limit mapping (CLI `-limited`).
- `disable_mapomatic`, `mapomatic_limit`: control mapomatic pass.
- `full_fidelity`: skip approximations.
- `allow_parameterized_merge`: allow parameterized merged gates when dumping pulses.
- `verbose`: verbosity/debug (matches CLI `-v`).

### TranspileResult
Outputs from `transpile_qasm`:
- `output_qasm`: QASM string.
- `output_qasm_path`: where QASM was written.
- `pulse_schedule`: pulse JSON as a string (if pulses generated).
- `pulse_doc`: parsed JSON (dict) when available.
- `pulse_path`: where the pulse JSON was written.
- `log`: captured stdout/diagnostics from the run.

### Functions
- `transpile_qasm(qasm_source, options=TranspileOptions())`: transpile from file path or raw QASM text; returns `TranspileResult`. Writes QASM (and pulses if template provided).
- `emit_qick(pulse_doc, qick_config, run_enabled=False, summary_only=False)`: emit pulses to QICK using in-memory documents (reuses `qasmtrans.qick_pulse_runner.run_program_obj`). Summary mode avoids hardware.
- `load_qick_config(path)`: load a QICK config JSON into a dict.

## Usage patterns

Minimal transpile + pulse dump:
```python
import json, qasmtrans as qt
opts = qt.TranspileOptions()
opts.backend_config = "data/devices/ibmq_toronto.json"
opts.pulse_template_path = "data/pulse_templates/toronto.json"
opts.output_path = "data/output_py"
res = qt.transpile_qasm("OPENQASM 2.0; ...", opts)
pulse_doc = res.pulse_doc or json.loads(res.pulse_schedule)
```

QICK summary (no hardware):
```python
qick_cfg = {"qubit_gen_map": {"0": 0, "1": 1}, "pulse_freq": 6.0e9}
summary = qt.emit_qick(pulse_doc, qick_cfg, run_enabled=False, summary_only=True)
print(summary)
```

Hardware attempt (requires `qick` on a QICK board):
```python
qt.emit_qick(pulse_doc, qick_cfg, run_enabled=True)  # may raise if qick or hardware unavailable
```

Plotting helper:
```python
import qasmtrans.plot_pulses as plot_pulses
analog, virtual, labels = plot_pulses.build_qubit_events(pulse_doc)
plot_pulses.plot_events(analog, virtual, labels, output_path="pulses.png", samples_per_us=200.0)
```

## Path resolution & packaging notes
- The shim first imports `qasmtrans_core` normally; if not found, it tries package-relative or bundled resource paths.
- CLI QICK emission also prepends `./python`, sibling/parent `python`, or `QASMTRANS_PYTHON_PATH` to `sys.path` (see `src/qick_emitter.cpp`).
- Wheels can bundle the extension inside `qasmtrans/`; the shim handles both layouts.
