#!/usr/bin/env python3
"""Run a fixed-layout routing-pressure sweep for QASMTrans vs Qiskit."""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

REPO_ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = REPO_ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

import qasmtrans  # noqa: E402


def count_gate_types(qc: QuantumCircuit) -> Tuple[int, int]:
    one = two = 0
    for inst, qargs, _ in qc.data:
        if inst.name in {"barrier", "measure", "delay", "reset"}:
            continue
        if len(qargs) == 1:
            one += 1
        elif len(qargs) == 2:
            two += 1
    return one, two


def load_device(cfg_path: Path) -> Dict[str, object]:
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    coupling = cfg.get("cx_coupling") or cfg.get("coupling_map")
    if not coupling:
        raise ValueError(f"No coupling map in {cfg_path}")
    edges = []
    for entry in coupling:
        if isinstance(entry, str):
            left, right = entry.replace("-", "_").split("_")
            edges.append((int(left), int(right)))
        else:
            edges.append((int(entry[0]), int(entry[1])))
    basis = list(dict.fromkeys(cfg.get("basis_gates") or ["rz", "sx", "x", "cx"]))
    return {
        "name": cfg.get("name") or cfg_path.stem,
        "map": CouplingMap(edges),
        "basis": basis,
    }


def parse_qasmtrans_swap_count(log_text: str) -> int | None:
    match = re.search(r"STEP-2\. Routing stats: .*swap=(\d+)", log_text)
    if not match:
        return None
    return int(match.group(1))


def run_qasmtrans_case(
    qasm_path: Path,
    device_path: Path,
    layout: List[int],
    optimize_1q: bool,
    optimize_2q_cancel: bool,
    optimize_commute_2q: bool,
    optimize_2q_synth: bool,
) -> Dict[str, object]:
    opts = qasmtrans.TranspileOptions()
    opts.mode = "ibmq"
    opts.backend_config = str(device_path)
    opts.disable_mapomatic = True
    opts.optimize_1q = optimize_1q
    opts.optimize_2q_cancel = optimize_2q_cancel
    opts.optimize_commute_2q = optimize_commute_2q
    opts.optimize_2q_synth = optimize_2q_synth
    opts.initial_layout = layout
    opts.verbose = 1

    start = time.perf_counter()
    res = qasmtrans.transpile_qasm(str(qasm_path), opts)
    elapsed_ms = (time.perf_counter() - start) * 1e3
    qc = QuantumCircuit.from_qasm_str(res.output_qasm)
    one, two = count_gate_types(qc)
    return {
        "time_ms": elapsed_ms,
        "depth": qc.depth(),
        "one_qubit": one,
        "two_qubit": two,
        "swap_count": parse_qasmtrans_swap_count(res.log),
        "logical_to_physical": json.dumps(list(res.logical_to_physical)),
    }


def run_qiskit_case(
    qasm_path: Path,
    device: Dict[str, object],
    layout: List[int],
    optimization_level: int,
    seed: int,
) -> Dict[str, object]:
    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    start = time.perf_counter()
    tc = transpile(
        qc,
        coupling_map=device["map"],
        basis_gates=device["basis"],
        optimization_level=optimization_level,
        seed_transpiler=seed,
        initial_layout=layout,
        layout_method="trivial",
        routing_method="sabre",
    )
    elapsed_ms = (time.perf_counter() - start) * 1e3
    one, two = count_gate_types(tc)
    return {
        "time_ms": elapsed_ms,
        "depth": tc.depth(),
        "one_qubit": one,
        "two_qubit": two,
    }


def iter_manifest_rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as fh:
        yield from csv.DictReader(fh)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark fixed-layout routing-pressure sweep.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/routing_pressure_manifest.csv"),
        help="Manifest produced by generate_routing_pressure_qasm.py",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/routing_pressure_bench.csv"),
        help="Where to write the benchmark results CSV.",
    )
    parser.add_argument("--seed", type=int, default=1, help="Qiskit transpiler seed.")
    parser.add_argument("--qiskit-level", type=int, default=3, help="Qiskit optimization level.")
    parser.add_argument("--optimize-1q", action="store_true", help="Enable QASMTrans 1q optimization.")
    parser.add_argument("--optimize-2q-cancel", action="store_true", help="Enable QASMTrans 2q cancellation.")
    parser.add_argument("--optimize-commute-2q", action="store_true", help="Enable QASMTrans commute pass.")
    parser.add_argument("--optimize-2q-synth", action="store_true", help="Enable QASMTrans 2q synthesis.")
    args = parser.parse_args()

    rows_out = []
    for row in iter_manifest_rows(args.manifest):
        qasm_path = Path(row["circuit_file"])
        device_path = Path(row["device_file"])
        layout = json.loads(row["initial_layout"])
        device = load_device(device_path)

        print(
            f"[{row['case_id']}] n={row['num_qubits']} pressure={row['routing_pressure']} layout={row['layout_tag']}",
            flush=True,
        )

        qt = run_qasmtrans_case(
            qasm_path,
            device_path,
            layout,
            optimize_1q=args.optimize_1q,
            optimize_2q_cancel=args.optimize_2q_cancel,
            optimize_commute_2q=args.optimize_commute_2q,
            optimize_2q_synth=args.optimize_2q_synth,
        )
        qk = run_qiskit_case(
            qasm_path,
            device,
            layout,
            optimization_level=args.qiskit_level,
            seed=args.seed,
        )

        rows_out.append(
            {
                **row,
                "qasmtrans_time_ms": f"{qt['time_ms']:.6f}",
                "qasmtrans_depth": qt["depth"],
                "qasmtrans_1q": qt["one_qubit"],
                "qasmtrans_2q": qt["two_qubit"],
                "qasmtrans_swap_count": qt["swap_count"] if qt["swap_count"] is not None else "",
                "qasmtrans_mapping": qt["logical_to_physical"],
                "qiskit_time_ms": f"{qk['time_ms']:.6f}",
                "qiskit_depth": qk["depth"],
                "qiskit_1q": qk["one_qubit"],
                "qiskit_2q": qk["two_qubit"],
                "depth_ratio_qt_over_qk": f"{qt['depth'] / qk['depth']:.6f}",
                "twoq_ratio_qt_over_qk": f"{qt['two_qubit'] / qk['two_qubit']:.6f}",
                "time_ratio_qt_over_qk": f"{qt['time_ms'] / qk['time_ms']:.6f}",
            }
        )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "case_id",
                "num_qubits",
                "reps",
                "circuit_file",
                "device_file",
                "layout_tag",
                "level_index",
                "routing_pressure",
                "max_hop",
                "initial_layout",
                "qasmtrans_time_ms",
                "qasmtrans_depth",
                "qasmtrans_1q",
                "qasmtrans_2q",
                "qasmtrans_swap_count",
                "qasmtrans_mapping",
                "qiskit_time_ms",
                "qiskit_depth",
                "qiskit_1q",
                "qiskit_2q",
                "depth_ratio_qt_over_qk",
                "twoq_ratio_qt_over_qk",
                "time_ratio_qt_over_qk",
            ],
        )
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote benchmark results to {args.output_csv}")


if __name__ == "__main__":
    main()
