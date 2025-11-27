# Benchmarking scripts

This directory mirrors the benchmarking helpers from `scripts/` but keeps all outputs local to avoid cluttering the repo root.

## Compilation benchmark (QASMTrans API)
Runs the Python API against a folder of QASM circuits, optionally comparing to Qiskit. Outputs a CSV and per-circuit logs here by default.

Example (skip Qiskit, 5s timeout, 3 retries):
```
python test/benchmarking/run_compilation_bench_api.py \
  --qasm_dir data/test_benchmark \
  --toronto_config data/devices/ibmq_toronto.json \
  --brisbane_config data/devices/ibm_brisbane.json \
  --mode ibmq \
  --max_time_ms 5000 \
  --retries 3 \
  --skip_qiskit
```
Outputs:
- `test/benchmarking/compilation_benchmarks.csv`
- Logs: `test/benchmarking/logs/`

## Adaptive timings (NWQ-Sim)
Python port of the adaptive timing shell script. Runs QASMTrans on AdaptiveData circuits (multiple copies per run), parses partition/total timings from stdout, and simulates subchips with NWQ-Sim for fidelity stats.

Example:
```
python test/benchmarking/run_adaptive_timings.py \
  --qasmtrans_bin build/QASMTrans \
  --nwqsim_exe ../NWQ-Sim/build/qasm/nwq_qasm
```
Default inputs target the AdaptiveData circuits living one level above the repo (`../AdaptiveData`). If you keep them elsewhere, point `--input_dir` at the correct folder. Outputs:
- CSV: `test/benchmarking/adaptive_timing_ibm_brisbane.csv`
- Transpiled outputs: `test/benchmarking/output_adaptive/`
Requires NWQ-Sim and the AdaptiveData QASM files. Adjust `--counts`/`--device_json`/`--mapomatic_limit` as needed.

Example override for a custom dataset:
```
python test/benchmarking/run_adaptive_timings.py \
  --input_dir ../AdaptiveData \
  --qasmtrans_bin build/QASMTrans \
  --nwqsim_exe ../NWQ-Sim/build/qasm/nwq_qasm
```
Outputs:
 - CSV: `test/benchmarking/adaptive_timing_ibm_brisbane.csv`
 - Transpiled outputs: `test/benchmarking/output_adaptive/`
Requires NWQ-Sim and the AdaptiveData QASM files. Adjust `--counts`/`--device_json`/`--mapomatic_limit` as needed.

## QASMTrans VQE/UCCSD metrics
Runs the QASMTrans Python API on generated VQE/UCCSD circuits and records wall/reported times plus simple metrics. Defaults to the generated circuits under `data/test_benchmark/vqe_uccsd_generated`.

Example (skip Qiskit timings):
```
python test/benchmarking/run_qasmtrans_vqe.py --skip_qiskit
```
Outputs:
- CSV: `test/benchmarking/vqe_uccsd_qasmtrans_metrics.csv`
Use `--device_json`/`--input_dir`/`--output_csv` to override paths; drop `--skip_qiskit` to collect Qiskit O1/O2/O3 timings too.

## Alg4 (Rigetti pulse merge) metrics
Collects fidelity/latency metrics for 4-qubit circuits (qft/shor/adder/w/bv/ghz) using the QASMTrans Python API plus the pulse merge/sim helpers.

Example:
```
python test/benchmarking/collect_alg4_metrics.py \
  --merge-limit 1 \
  --base-device data/devices/rigetti_ankaa4q_ideal_device.json \
  --base-pulses data/devices/rigetti_ankaa4q_ideal_pulses.json
```
Outputs stay under `test/benchmarking/` by default:
- CSV: `test/benchmarking/alg4_metrics.csv`
- Per-algorithm artefacts: `test/benchmarking/output_alg4/<algo>/`
Requires the baseline Rigetti device/pulse JSONs (generate with `device_pulse_fab/build_rigetti_ankaa9q_configs.py`), and uses `device_pulse_fab/merge_pulse.py` + `simulate_ankaa9q.py` under the hood.

## Mapomatic fidelity/timing suite
Runs QASMTrans (Python API) across Mapomatic modes (nomap/product/full/hybrid), simulates with NWQ-Sim to extract fidelity, and generates plots. Requires NWQ-Sim built at `../NWQ-Sim/build/qasm/nwq_qasm` (override with `--nwqsim_exe`).

Example (default circuits from `scripts/tmp`):
```
python test/benchmarking/run_mapomatic_suite.py \
  --nwqsim_exe ../NWQ-Sim/build/qasm/nwq_qasm
```
Outputs:
- `test/benchmarking/mapomatic_benchmarks.csv`
- Transpiled outputs: `test/benchmarking/output_qasm_file/`
- Plots: `test/benchmarking/mapomatic_fidelity.pdf` and `mapomatic_timing.pdf`

Use `--skip_plots` to skip PDF generation. Use `--modes nomap,product,full,hybrid` to select modes. Use `-i my.qasm` to override circuits.
