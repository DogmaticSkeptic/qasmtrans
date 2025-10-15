# Rigetti Pulse Simulator Overview

This repository includes a lightweight Python driver (`device_pulse_fab/simulate_ankaa9q.py`) that replays pulse schedules emitted by QASMTrans against a calibrated Rigetti Ankaa-style device model. The simulator uses QuTiP throughout; every fidelity number reported in the CLI output or JSON summaries comes directly from QuTiP density-matrix evolutions—there are no analytical shortcuts or surrogate scoring functions.

## High-Level Flow

1. **Prepare calibrations**  
   `device_pulse_fab/build_rigetti_ankaa9q_configs.py` extracts the target device metadata (pulse envelopes, gate durations, target fidelities) and generates:
   - a device JSON (`rigetti_ankaaNq_device.json`) that records hardware topology, gate errors, and basis gates;  
   - a pulse template (`rigetti_ankaaNq_pulses.json`) containing sampled I/Q waveforms, areas, and the target rotation (`parameters.theta`) for each calibrated gate.

2. **Transpile the circuit**  
   `./build/QASMTrans -m rigetti -c <device> -p <pulses>` converts a logical QASM file into the calibrated basis, writes the routed QASM, and produces a pulse schedule (`*_pulses.json`) whose events reference the pulse library IDs.

3. **Replay the schedule with QuTiP**  
   `python device_pulse_fab/simulate_ankaa9q.py --pulse … --device … --qasm …` performs the following:
   - Canonicalises the pulse JSON (normalises gate labels, verifies samples).
   - Builds an “ideal progression” by applying each pulse’s *target* unitary (the `parameters.theta` stored in the pulse library) to a logical statevector, keeping track of virtual-Z frame updates.
   - Constructs the lab-frame Hamiltonian: single-qubit drives come from the rotated I/Q samples, and two-qubit exchange terms come from the `samples_i` array for the calibrated √iSWAP pulses.  
   - Calls `qutip.mesolve` with the assembled Hamiltonian list and the chosen initial density matrix (usually `|000…0⟩`; the sanity checks use |10⟩ for exchange gates) to obtain the full time evolution.  
   - Compares the simulated state to the ideal state after every event (if `-v 1`) and reports the overall fidelity plus any residual virtual-Z phase corrections applied at readout.

   The per-pulse “sanity check” table is also QuTiP-backed: each entry re-runs `simulate_pulse_schedule` on a single pulse in isolation and measures its fidelity against the target unitary.

4. **Reporting and tooling**  
   With `--output`, the driver saves a JSON report comprising the fidelity, number of pulses, and input file references—handy for batch sweeps. Optional CLI verbosity (`-v 1`) prints per-event fidelity, Bloch vectors, and inferred rotation axes to help debug calibration mismatches.

## Key Implementation Details

- **Virtual-Z handling**: Virtual frame updates never generate lab-frame Hamiltonian terms. The simulator stores per-qubit phase accumulators, applies them when rotating drive waveforms, and keeps the logical and physical frames aligned for fidelity comparisons.
- **Target vs. actual rotation**: Each pulse definition contains both the calibrated “actual” area (`calibration.theta_actual_rad`) and the intended logical rotation (`parameters.theta`). The simulator uses the latter when building the ideal state progression so we can observe how calibration errors accrue over a long schedule.
- **Initial states for sanity tests**: Two-qubit pulses are evaluated on the |10⟩ input so that √iSWAP fidelity is sensitive to the exchange leakage and phase errors; single-qubit pulses start from |0⟩.
- **No short-cuts**: Every fidelity number visible in the CLI output, sanity checks, and JSON reports comes from `qutip.fidelity` applied to the propagated density matrices. There is no analytical shortcutting or table look-up.

## Reproducing Results

```bash
# 1. Generate calibrations (example for 4 qubits)
conda run -n qasmtrans python device_pulse_fab/build_rigetti_ankaa9q_configs.py \
    --num-qubits 4 \
    --out-device data/devices/rigetti_ankaa4q_device.json \
    --out-pulses data/devices/rigetti_ankaa4q_pulses.json

# 2. Transpile a circuit (Hadamards, GHZ, brickwork, …)
./build/QASMTrans \
    -i data/test_benchmark/brickwork4d3.qasm \
    -c data/devices/rigetti_ankaa4q_device.json \
    -p data/devices/rigetti_ankaa4q_pulses.json \
    -m rigetti \
    -o data/output/brickwork4d3_rigetti.qasm

# 3. Replay the pulse schedule with QuTiP
conda run -n qasmtrans python device_pulse_fab/simulate_ankaa9q.py \
    --pulse data/output/brickwork4d3_rigetti_pulses.json \
    --device data/devices/rigetti_ankaa4q_device.json \
    --qasm data/output/brickwork4d3_rigetti.qasm \
    -v 1 \
    --output data/output/brickwork4d3_sim.json
```

Check the saved JSON report(s) or the terminal output to inspect the fidelity, number of pulses, and per-event diagnostics.
