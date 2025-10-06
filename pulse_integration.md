# Pulse Integration Notes

This document outlines the JSON artifacts produced and consumed when connecting QASMTrans pulse dumps to the QICK runner.

## Pulse Template (`*.pulses.json`)

Pulse templates describe the basis-library waveforms that QASMTrans can reference while compiling circuits. They are typically authored once per backend and passed to `qasmtrans` via `-p`.

```json
{
  "name": "ibmq_dummy",
  "version": "0.1",
  "num_qubits": 12,
  "basis_gates": ["id", "rz", "sx", "x", "cx", "reset"],
  "pulse_definitions": [
    {
      "id": "sx_q0",                  // optional; generated if omitted
      "gate": "sx",                    // gate name (lowercase)
      "qubits": [0],                    // participating qubits
      "shape": "gaussian",             // legacy label (optional)
      "waveform_type": "gaussian",     // explicit envelope intent
      "width": 5.7e-09,                 // seconds; used for scheduling
      "amplitude": 0.1,                 // envelope scale when samples are omitted
      "parameters": {"theta": 1.5708}, // optional constrained parameters
      "samples_i": [0.0, 0.02, 0.08],   // optional DAC I samples
      "samples_q": [0.0, 0.01, 0.03],   // optional DAC Q samples
      "virtual": false                  // optional; mark frame-change pulses when true
    },
    "... more entries ..."
  ]
}
```

Fields beyond `pulse_definitions` (e.g., coherence data) are ignored by the pulse exporter but can be retained for bookkeeping. Every gate/qubit combination referenced by the transpiled circuit must appear in the template; missing entries now raise an error during pulse dumping instead of generating fallbacks. Use `waveform_type` to declare how the pulse should be interpreted (`gaussian`, `flat_top`, `constant`, `arbitrary`, etc.). When `waveform_type` is `arbitrary`, `samples_i`/`samples_q` must be provided.

## Circuit Pulse Dump (`<stem>_pulses.json`)

Running `qasmtrans -p template.json` alongside the usual chip config produces a circuit-specific pulse dump next to the compiled QASM. The file has three top-level sections:

```json
{
  "backend": {
    "name": "ibmq_dummy",          // copied from backend/template
    "version": "0.1",
    "input_qasm": "bv10.qasm",
    "backend_config": "../data/devices/backend_config.json",
    "pulse_template": "../data/devices/dummy_ibmq12_pulses.json",
    "num_qubits": 5,
    "total_pulses": 42,
    "total_duration": 8.4e-06,
    "auto_generated_pulses": ["cx:0_1"] // only when fallbacks were needed
  },
  "pulse_library": [
    {
      "id": "sx_q0",
      "gate": "sx",
      "qubits": [0],
      "shape": "gaussian",
      "width": 5.7e-09,
      "amplitude": 0.1,
      "samples_i": [0.0, 0.02, ...],
      "samples_q": [0.0, 0.01, ...]
    },
    {
      "id": "rz_q0",
      "gate": "rz",
      "qubits": [0],
      "waveform_type": "virtual",
      "width": 0.0,
      "amplitude": 0.0,
      "virtual": true
    },
    {
      "id": "rx_q0_pi",
      "gate": "rx",
      "qubits": [0],
      "waveform_type": "arbitrary",
      "width": 4.0e-08,
      "amplitude": 0.34,
      "parameters": {"theta": 3.14159},
      "samples_i": [0.02, 0.03, ...],
      "samples_q": [-0.00, -0.01, ...]
    },
    "... unique pulses used in this circuit ..."
  ],
  "schedule": [
    {
      "index": 0,
      "gate": "sx",
      "qubits": [0],
      "pulse_id": "sx_q0",
      "start_time": 0.0,
      "duration": 5.7e-09,
      "parameters": {"theta": 1.5708} // optional
    },
    {
      "index": 1,
      "gate": "cx",
      "qubits": [0, 1],
      "pulse_id": "cx_q0_q1",
      "start_time": 5.7e-09,
      "duration": 3.2e-08
    }
  ]
}
```

- `pulse_library` collapses the template to only the pulses actually referenced by the circuit, preserving `samples_i`/`samples_q` arrays, `waveform_type`, and any constrained `parameters` (e.g., `theta` for Rigetti `rx`). Downstream tools can render exact DAC shapes or analytic envelopes (`gaussian`, `flat_top`, `constant`, `arbitrary`). Multi-qubit entries list every participating qubit so the runner can trigger the mapped generator channels simultaneously, and arbitrary/flat-top pulses ferry their envelope samples for upload to the QICK waveform memory.
- `schedule` is time-ordered and ensures single-qubit contention is serialized via start times chosen during export.
- Virtual frame-change operations (e.g., `rz`) appear with `"virtual": true`, zero duration, and no waveforms. The runner and plotting utilities skip them while the schedule retains their ordering and parameter metadata.

## QICK Configuration (`qick_config.json`)

The QICK runner consumes a hardware-specific configuration describing channel assignments and acquisition parameters. The exact schema mirrors what is expected by `TranspiledPulseProgram` and the underlying QICK API.

```json
{
  "qubit_gen_map": {"0": 6, "1": 7},    // pulse driver channel per qubit
  "generators": [6, 7],                   // optional; recomputed if omitted
  "ro_chs": [0],                          // readout ADC channels
  "readout_gen_ch": 8,                    // generator driving readout (optional)
  "pulse_freq": 250.0,                    // MHz unless otherwise noted by hardware
  "readout_freq": 250.0,
  "qubit_freqs": {"0": 240.0, "1": 255.0},
  "pulse_gain_scale": 30000,
  "default_gain": 28000,
  "gain_max": 32767,
  "start_offset": 0.0,                    // seconds offset applied before playback
  "default_phase_deg": 0.0,
  "init_synci": 200,
  "relax_delay": 1.5,                     // microseconds
  "adc_trig_offset": 100,
  "readout_length": 200,
  "enable_readout": true,
  "scope_pin": true
}
```

Only the keys used by `qick_pulse_runner.py` are required, so the document can be pared down to the essentials (e.g., `qubit_gen_map`, frequency/gain scales, readout settings).

## End-to-End Flow

1. **Author template**: Create or update the pulse template JSON describing the basis gate waveforms (`data/devices/dummy_ibmq12_pulses.json` provides a scaffold).
2. **Transpile**: Run QASMTrans with both the chip backend (`-c`) and pulse template (`-p`). This produces the routed QASM and `<stem>_pulses.json`, containing the pulse schedule tailored to the circuit and hardware timing data.
3. **Prepare QICK config**: Maintain a `qick_config.json` per board capturing generator, frequency, and acquisition settings.
4. **Play pulses**: Either invoke the standalone runner (`python3 python/qick_pulse_runner.py <stem>_pulses.json qick_config.json --summary`) or trigger emission directly during transpilation with `./QASMTrans -i circuit.qasm -c backend.json -p template.json -e qick_config.json [--emit-run]`. Without `--emit-run` the executable prints the summary and exits; adding the flag streams the pulses to the configured QICK hardware in the same invocation.

## Embedded QICK Emission (`-e`)

The C++ driver links against `pybind11::embed` so pulse playback can happen inside a single `QASMTrans` invocation:

```bash
./QASMTrans -i examples/ghz.qasm \
            -c data/devices/dummy_ibmq12.json \
            -p data/devices/dummy_ibmq12_pulses.json \
            -e qick_config.json           # emit summary

./QASMTrans ... -e qick_config.json --emit-run   # emit and stream to hardware
```

- `-e` requires `-p`; the transpiler must produce `<stem>_pulses.json` before emission.
- If `--emit-run` is omitted the run stays in “summary only” mode (counts, duration, generator map). With the flag, the embedded runner calls `TranspiledPulseProgram` exactly like the Python CLI.
- The pybind shim auto-adds `<repo>/python` to `sys.path`, and honours `QASMTRANS_PYTHON_PATH` when the helper scripts live elsewhere.
- Errors raised by the Python runner bubble back through stdout/stderr so the CLI exit code reflects hardware/programming failures.

### Frame-Change Support

During `prepare_events` the Python runner now accumulates RZ angles from virtual pulses and advances the generator phase for every subsequent waveform on the affected qubits. The integrated pathway shares the same code path, so hardware playback preserves QASM frame changes without manually expanding them into physical pulses.

For quick visualization without rerunning the transpiler, `python/test_circuit_pulses.json` provides a canned pulse dump containing constant, flat-top, Gaussian, and arbitrary waveforms (with I/Q samples) that you can feed to the plotter or runner. Two-qubit entries (e.g., `cx_q0_q1`) list both participating qubits and the shared `pulse_id`; during playback the runner resolves each qubit through `qubit_gen_map`, programming both generators and issuing the `pulse()` calls in the same body cycle so the channels fire in lockstep. Try it out with:

```bash
python3 python/plot_pulses.py python/test_circuit_pulses.json --title "Waveform Demo"
python3 python/qick_pulse_runner.py python/test_circuit_pulses.json qick_config.json --summary
```

The plot highlights per-qubit DAC traces, while the runner verifies that multi-qubit pulses fan out to every mapped generator simultaneously.

The integration hinges on matching the `pulse_id` names in `schedule` with entries in the `pulse_library`, and mapping the first qubit in each scheduled operation to a valid generator channel in the QICK config. With these artifacts in place, the pipeline carries a QASM input through transpilation, pulse construction, and hardware playback without additional manual transformations.

### Critical Path Metadata

Every pulse dump now records a `critical_path` object summarizing the longest serialized chain of circuit operations:

```json
"critical_path": {
  "total_duration": 8.40e-06,
  "gate_count": 42,
  "path_indices": [0, 1, 4, ...],
  "top_gates": [
    {"gate": "cx", "count": 12, "total_duration": 4.2e-06},
    {"gate": "rz", "count": 8, "total_duration": 0.0},
    {"gate": "sx", "count": 10, "total_duration": 2.1e-06}
  ]
}
```

`path_indices` reference entries in the `schedule` array; `top_gates` lists up to the ten gate types contributing the most cumulative duration along that path, together with their occurrence counts.
