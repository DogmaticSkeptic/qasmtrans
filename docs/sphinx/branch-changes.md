# Branch changes vs master

This branch adds pulses/QICK support, new routing heuristics, a Python package, expanded benchmarks, and supporting tooling. Highlights and file anchors follow so reviewers can jump to the code.

## Pulses and QICK
- Pulse JSON generation: `include/dump_pulses.hpp`, `src/cli_support.cpp`, `src/qasmtrans.cpp` (`-p` flag derives `<stem>_pulses.json`).
- QICK emission: `include/qick_emitter.hpp`, `src/qick_emitter.cpp`, Python runner `python/qasmtrans/qick_pulse_runner.py` (CLI `-e/--emit-run`).
- Demos/utilities: `python/demo_qick.py`, plotting helper `python/qasmtrans/plot_pulses.py`, docs `pulse_integration.md`, `PULSE_SIMULATOR.md`.

## Routing/Mapomatic and IR changes
- New mapping pass and heuristics: `include/circuit_passes/mapomatic.hpp`, scripts `scripts/run_mapomatic_suite.py`, `scripts/compare_mapomatic_fidelity.sh`, plotting scripts.
- Routing/mapping refactors: `include/circuit_passes/routing_mapping.hpp`, `include/circuit_passes/transpiler.hpp`, `include/circuit_passes/decompose.hpp`.
- IR updates: `include/IR/{chip,circuit,gate}.hpp`, new utilities `include/util/chip_partition.hpp`, graph support via vendored lemon.
- Lemon graph library vendored: `third_party/lemon/...` plus `include/IR/graph.hpp` hooks.

## Python package and bindings
- Pybind module expanded: `src/python_binding.cpp` now exposes `transpile_qasm`, `emit_qick`, `load_qick_config`, `TranspileOptions/Result`.
- Packaging shim: `python/qasmtrans/__init__.py`, wheel scaffolding `pyproject.toml`, `requirements*.txt`.
- CLI/runner reuse: bindings share pulse dumping and QICK runner logic; demos in `python/demo.py` and `python/demo_qick.py`.

## Benchmarks and datasets
- Large benchmark expansion under `data/test_benchmark/` (brickwork, GHZ, QAOA, VQE, QUGAN, UCCSD, etc.).
- Benchmark scripts: `scripts/run_compilation_bench.py`, `scripts/run_compilation_bench_api.py`, `scripts/run_qasmtrans_vqe.py`, `scripts/run_adaptive_timings.sh`, and plotting scripts.
- Test harnesses mirror scripts in `test/benchmarking/`.

## Device pulse fabrication utilities
- Pulse synthesis and Rigetti helpers: `device_pulse_fab/build_rigetti_ankaa9q_configs.py`, `device_pulse_fab/merge_pulse.py`, `device_pulse_fab/simulate_ankaa9q.py`.

## Docs, notebooks, and CI
- Sphinx+Doxygen docs: `docs/` (including this page), MkDocs shim `mkdocs.yml`, GitHub Pages workflow `.github/workflows/docs.yml`.
- Notebooks for demos: `notebooks/demo.ipynb`, `notebooks/qick_demo.ipynb`.
- README and requirements updated to cover docs/pulses/Python bindings.

## Build/tooling changes
- `CMakeLists.txt` updated for pybind targets and doc tooling.
- New requirements files and `.gitignore` entries.
- Removed legacy `data/test_benchmark/test.qasm`; new configs and scripts added across `scripts/`.
