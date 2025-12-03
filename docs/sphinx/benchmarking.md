# Benchmarking suite

This branch adds a full benchmarking stack: new datasets, batch drivers, and plotting utilities. Use these pointers to review or extend runs.

## Datasets
- Expanded QASM corpora in `data/test_benchmark/`: adder variants, brickwork, GHZ, QAOA (logical and rxyz), QFT, BV, VQE/VQE-UCCSD families, QUGAN, and more. Larger generated sets live under `data/test_benchmark/vqe_uccsd_generated/`.
- Experimental mini-circuits under `data/test_benchmark/experiments/` for smoke testing and demos.
- Legacy `test.qasm` was removed; expect new stems in scripts and outputs.

## Batch runners (scripts/)
- `run_compilation_bench.py` / `run_compilation_bench.sh`: drive QASMTrans across the corpus; options for backend/mode, output dirs, and metrics.
- `run_compilation_bench_api.py`: Python API variant that uses `qasmtrans_core` bindings to collect results programmatically.
- `run_qasmtrans_vqe.py` / `run_qasmtrans_vqe.sh`: VQE-focused sweeps, including generated UCCSD inputs.
- `run_adaptive_timings.sh`: measure timing under different settings.
- `run_mapomatic_suite.py` / `run_mapomatic_bench.sh` / `compare_mapomatic_fidelity.sh`: evaluate the new mapomatic pass and plot fidelity/timing.
- `run_ansatz_merge_sweep.sh`, `run_compilation_bench_api.py`: utilities for merged-gate experiments.
- `collect_alg4_metrics.sh`: helper to gather metrics for Alg4 runs.

## Plotting utilities (scripts/)
- `plot_mapomatic_fidelity.py`, `plot_mapomatic_timing.py`: visualize mapomatic results.
- `plot_alg4_metrics.py`, `plot_vqe_timings.py`: plot algorithm- and VQE-specific metrics.
- `qiskit_transpile.py`: helper to compare with Qiskit transpilation baselines.

## Test mirrors (test/benchmarking/)
- Python harnesses mirror the scripts: `run_compilation_bench_api.py`, `run_mapomatic_suite.py`, `run_adaptive_timings.py`, `run_qasmtrans_vqe.py`, plus `README.md` describing usage.

## Outputs and conventions
- Scripts typically emit QASM/pulse outputs under `data/output_*` (configurable), plus JSON/CSV summaries for plotting.
- Mapomatic and VQE scripts accept backend configs (`-c`), modes (`-m`), and toggles for limits/merged gates; see script `--help` for arguments.

## Review guide
- For performance/latency logic changes, see `include/circuit_passes/mapomatic.hpp` and `include/circuit_passes/routing_mapping.hpp`.
- For data provenance, inspect individual QASM files and generated UCCSD sets under `data/test_benchmark/vqe_uccsd_generated/`.
