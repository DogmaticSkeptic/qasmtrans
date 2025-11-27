#!/usr/bin/env python3
"""
Run QASMTrans (Python API) on generated VQE/UCCSD circuits and collect timings/metrics.
Outputs stay in this benchmarking directory by default.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from pathlib import Path
from typing import Dict, Tuple

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

import qasmtrans


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]


def load_device(cfg_path: Path) -> Dict[str, object]:
    with cfg_path.open("r", encoding="utf-8") as fh:
        cfg = json.load(fh)

    coupling = cfg.get("cx_coupling") or cfg.get("coupling_map")
    if not coupling:
        raise ValueError(f"No coupling_map/cx_coupling found in {cfg_path}")

    edges = []
    for entry in coupling:
        if isinstance(entry, str):
            parts = entry.replace("-", "_").split("_")
            if len(parts) != 2:
                continue
            u, v = map(int, parts)
        else:
            if len(entry) != 2:
                continue
            u, v = entry
        edges.append((u, v))

    cmap = CouplingMap(edges)
    basis = cfg.get("basis_gates") or ["rz", "sx", "x", "cx"]
    basis_gates = list(dict.fromkeys(basis))

    return {
        "map": cmap,
        "basis": basis_gates,
        "name": cfg.get("name") or cfg_path.stem,
    }


def run_qasmtrans(circuit_qasm: str, device_json: Path) -> Tuple[float, float, Tuple[str, str, str]]:
    opts = qasmtrans.TranspileOptions()
    opts.mode = "ibmq"
    opts.backend_config = str(device_json)
    opts.disable_mapomatic = True
    opts.verbose = 1

    start = time.perf_counter()
    res = qasmtrans.transpile_qasm(circuit_qasm, opts)
    elapsed_ms = (time.perf_counter() - start) * 1e3

    reported_ms = None
    one_qubit = two_qubit = depth = ""
    if isinstance(res.log, str):
        match = re.search(r"total QASMTrans time:\s*([-+]?\d*\.?\d+)", res.log)
        if match:
            reported_ms = float(match.group(1))

        metrics_match = re.search(
            r"\[metrics\]\s*one_qubit_gates=(\d+)\s+two_qubit_gates=(\d+)\s+depth=(\d+)",
            res.log,
        )
        if metrics_match:
            one_qubit, two_qubit, depth = metrics_match.groups()

    return elapsed_ms, reported_ms if reported_ms is not None else elapsed_ms, (one_qubit, two_qubit, depth)


def run_qiskit_times(qasm_path: Path, device: Dict[str, object]) -> Tuple[str, str, str]:
    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    times = []
    for level in (1, 2, 3):
        start = time.perf_counter()
        transpile(qc, coupling_map=device["map"], basis_gates=device["basis"], optimization_level=level)
        elapsed_ms = (time.perf_counter() - start) * 1e3
        times.append(f"{elapsed_ms:.3f}")
    return tuple(times)  # type: ignore


def main() -> None:
    parser = argparse.ArgumentParser(description="QASMTrans VQE/UCCSD benchmark via Python API.")
    parser.add_argument(
        "--input_dir",
        default=str(REPO_ROOT / "data/test_benchmark/vqe_uccsd_generated"),
        help="Directory containing generated VQE/UCCSD QASM files.",
    )
    parser.add_argument(
        "--device_json",
        default=str(REPO_ROOT / "data/devices/ibm_brisbane.json"),
        help="Device JSON.",
    )
    parser.add_argument(
        "--output_csv",
        default=str(SCRIPT_DIR / "vqe_uccsd_qasmtrans_metrics.csv"),
        help="Destination CSV (default: test/benchmarking/vqe_uccsd_qasmtrans_metrics.csv).",
    )
    parser.add_argument(
        "--skip_qiskit",
        action="store_true",
        help="Skip Qiskit timings.",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    device_json = Path(args.device_json).resolve()
    device = load_device(device_json)

    circuits = sorted(input_dir.glob("*.qasm"))
    print(f"Found {len(circuits)} circuits under {input_dir}")
    if not circuits:
        return

    out_csv = Path(args.output_csv).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "circuit",
                "qasmtrans_reported_ms",
                "qasmtrans_elapsed_ms",
                "qiskit_o1_time_ms",
                "qiskit_o2_time_ms",
                "qiskit_o3_time_ms",
                "one_qubit",
                "two_qubit",
                "depth",
            ]
        )

        for circuit_path in circuits:
            name = circuit_path.name
            print(f"=== {name} ===")
            qasm_text = circuit_path.read_text(encoding="utf-8")
            elapsed_ms, reported_ms, metrics = run_qasmtrans(qasm_text, device_json)
            one_qubit, two_qubit, depth = metrics

            if args.skip_qiskit:
                qiskit_times = ("", "", "")
            else:
                qiskit_times = run_qiskit_times(circuit_path, device)

            writer.writerow(
                [
                    name,
                    f"{reported_ms:.3f}",
                    f"{elapsed_ms:.3f}",
                    qiskit_times[0],
                    qiskit_times[1],
                    qiskit_times[2],
                    one_qubit,
                    two_qubit,
                    depth,
                ]
            )

    print(f"Saved metrics to {out_csv}")


if __name__ == "__main__":
    main()
