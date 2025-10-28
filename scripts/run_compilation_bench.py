#!/usr/bin/env python3
"""
Benchmark runner that collects Qiskit (O1/O2/O3) and QASMTrans metrics for all
OpenQASM circuits under data/test_benchmark, retrying QASMTrans runs that exceed
the time budget. Results are written to a single CSV.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap


def load_device(cfg_path: Path) -> Dict[str, object]:
    with cfg_path.open("r") as fh:
        cfg = json.load(fh)

    coupling = cfg.get("cx_coupling") or cfg.get("coupling_map")
    if not coupling:
        raise ValueError(f"No coupling_map/cx_coupling found in {cfg_path}")

    edges = []
    for entry in coupling:
        if isinstance(entry, str):
            parts = entry.replace("-", "_").split("_")
            if len(parts) != 2:
                raise ValueError(f"Unexpected coupling entry {entry!r} in {cfg_path}")
            u, v = map(int, parts)
        else:
            if len(entry) != 2:
                raise ValueError(f"Unexpected coupling entry {entry!r} in {cfg_path}")
            u, v = entry
        edges.append((u, v))

    cmap = CouplingMap(edges)
    basis = cfg.get("basis_gates") or ["rz", "sx", "x", "cx"]
    basis_gates = list(dict.fromkeys(basis))

    return {
        "map": cmap,
        "basis": basis_gates,
        "name": cfg.get("name") or cfg_path.stem,
        "num_qubits": int(cfg.get("num_qubits", cmap.size())),
    }


def transpile_with_levels(
    qc: QuantumCircuit,
    device: Dict[str, object],
    levels: Sequence[int],
    seed: int,
) -> Dict[int, Dict[str, float]]:
    results: Dict[int, Dict[str, float]] = {}
    for lvl in levels:
        start = time.perf_counter()
        tc = transpile(
            qc,
            coupling_map=device["map"],
            basis_gates=device["basis"],
            optimization_level=lvl,
            seed_transpiler=seed,
        )
        elapsed_ms = (time.perf_counter() - start) * 1e3
        one, two = count_gate_types(tc)
        results[lvl] = {
            "time_ms": elapsed_ms,
            "depth": tc.depth(),
            "one_qubit": one,
            "two_qubit": two,
        }
    return results


def count_gate_types(qc: QuantumCircuit) -> tuple[int, int]:
    one = two = 0
    for inst, qargs, _ in qc.data:
        if inst.name in {"barrier", "measure", "delay"}:
            continue
        qubit_count = len(qargs)
        if qubit_count == 1:
            one += 1
        elif qubit_count == 2:
            two += 1
    return one, two


METRICS_RE = re.compile(
    r"\[metrics\]\s*one_qubit_gates=(\d+)\s+two_qubit_gates=(\d+)\s+depth=(\d+)"
)
FAST_TOTAL_RE = re.compile(r"total QASMTrans time:\s*([-+]?\d*\.?\d+)\s*ms")
SLOW_TOTAL_RE = re.compile(r"\[timing\]\s*total_ms=(\d+)")


def run_qasmtrans(
    exe: Path,
    circuit: Path,
    device_cfg: Path,
    backend_name: str,
    output_dir: Path,
    max_time_ms: int,
    retries: int,
) -> Dict[str, Optional[float]]:
    timeout_sec = max_time_ms / 1000.0
    for attempt in range(1, retries + 1):
        tmp_qasm = output_dir / f"{circuit.stem}_attempt{attempt}.qasm"
        cmd = [
            str(exe),
            "-i",
            str(circuit),
            "-m",
            "ibmq",
            "-c",
            str(device_cfg),
            "-o",
            str(tmp_qasm),
            "-v",
            "1",
            "--disable_mapomatic",
        ]
        print(f"    QASMTrans attempt {attempt}: {' '.join(cmd)}", flush=True)
        start = time.perf_counter()
        duration_ms = None
        try:
            result = subprocess.run(
                cmd,
                text=True,
                capture_output=True,
                check=False,
                timeout=timeout_sec,
            )
            output = result.stdout or ""
            duration_ms = (time.perf_counter() - start) * 1e3
        except subprocess.TimeoutExpired:
            print(
                f"      Terminated after {max_time_ms} ms timeout.",
                flush=True,
            )
            continue

        metrics_match = METRICS_RE.search(output)
        total_match = FAST_TOTAL_RE.search(output) or SLOW_TOTAL_RE.search(output)

        if not metrics_match or not total_match or result.returncode != 0:
            print("      QASMTrans failed or missing metrics.", flush=True)
            if result.stderr:
                print(result.stderr.strip(), flush=True)
            if attempt == retries:
                return {
                    "time_ms": None,
                    "one_qubit": None,
                    "two_qubit": None,
                    "depth": None,
                    "backend": backend_name,
                }
            continue

        total_ms = float(total_match.group(1))

        one_qubit, two_qubit, depth = map(int, metrics_match.groups())
        print(
            f"      Success: {total_ms:.1f} ms (elapsed {duration_ms:.1f} ms), 1q={one_qubit}, 2q={two_qubit}, depth={depth}",
            flush=True,
        )
        return {
            "time_ms": total_ms,
            "one_qubit": one_qubit,
            "two_qubit": two_qubit,
            "depth": depth,
            "backend": backend_name,
        }

    return {
        "time_ms": None,
        "one_qubit": None,
        "two_qubit": None,
        "depth": None,
        "backend": backend_name,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Qiskit + QASMTrans benchmarks.")
    parser.add_argument(
        "--qasm_dir",
        default="data/test_benchmark",
        help="Directory containing *.qasm circuits to benchmark.",
    )
    parser.add_argument(
        "--toronto_config",
        default="data/devices/ibmq_toronto.json",
        help="Device JSON for <=27 qubits.",
    )
    parser.add_argument(
        "--brisbane_config",
        default="data/devices/ibm_brisbane.json",
        help="Device JSON for >27 qubits.",
    )
    parser.add_argument(
        "--qasmtrans_bin",
        default="build/QASMTrans",
        help="Path to QASMTrans executable.",
    )
    parser.add_argument(
        "--max_qubits",
        type=int,
        default=150,
        help="Skip circuits with more than this many qubits.",
    )
    parser.add_argument(
        "--max_time_ms",
        type=int,
        default=1000,
        help="Maximum acceptable QASMTrans time before retry/skip.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=2,
        help="Number of QASMTrans attempts before giving up.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Seed passed to Qiskit's transpiler.",
    )
    parser.add_argument(
        "--skip_qiskit",
        action="store_true",
        help="Skip Qiskit transpilation and only collect QASMTrans metrics.",
    )
    parser.add_argument(
        "--output_csv",
        default="data/compilation_benchmarks.csv",
        help="Destination CSV file.",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    qasm_dir = (repo_root / args.qasm_dir).resolve()
    if not qasm_dir.is_dir():
        raise FileNotFoundError(f"QASM directory not found: {qasm_dir}")

    toronto_cfg = (repo_root / args.toronto_config).resolve()
    brisbane_cfg = (repo_root / args.brisbane_config).resolve()
    qasmtrans_bin = (repo_root / args.qasmtrans_bin).resolve()
    if not qasmtrans_bin.is_file():
        raise FileNotFoundError(f"QASMTrans binary not found: {qasmtrans_bin}")

    devices = {
        "toronto": load_device(toronto_cfg),
        "brisbane": load_device(brisbane_cfg),
    }

    qasm_files = sorted(qasm_dir.glob("*.qasm"))
    print(f"Discovered {len(qasm_files)} circuit(s) under {qasm_dir}.", flush=True)
    if not qasm_files:
        return

    rows: List[Dict[str, object]] = []

    tmp_parent = repo_root / "tmp"
    tmp_parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="qasmtrans_", dir=tmp_parent) as tmp:
        tmp_dir = Path(tmp)
        total = len(qasm_files)
        for idx, qasm_path in enumerate(qasm_files, start=1):
            qc = QuantumCircuit.from_qasm_file(str(qasm_path))
            qubits = qc.num_qubits
            if args.max_qubits and qubits > args.max_qubits:
                print(
                    f"[{idx}/{total}] Skipping {qasm_path.name}: {qubits} qubits (> {args.max_qubits}).",
                    flush=True,
                )
                continue

            device_key = "toronto" if qubits <= 27 else "brisbane"
            device = devices[device_key]
            print(
                f"[{idx}/{total}] {qasm_path.name}: {qubits} qubits → {device['name']} ({device['num_qubits']} phys)",
                flush=True,
            )

            if args.skip_qiskit:
                print("  Qiskit transpilation skipped (--skip_qiskit).", flush=True)
                qiskit_metrics = {
                    1: {"time_ms": None, "depth": None, "one_qubit": None, "two_qubit": None},
                    2: {"time_ms": None, "depth": None, "one_qubit": None, "two_qubit": None},
                    3: {"time_ms": None, "depth": None, "one_qubit": None, "two_qubit": None},
                }
            else:
                qiskit_metrics = transpile_with_levels(
                    qc, device, levels=[1, 2, 3], seed=args.seed
                )
                for level in [1, 2, 3]:
                    res = qiskit_metrics[level]
                    print(
                        f"  Qiskit level {level}: {res['time_ms']:.3f} ms, depth={res['depth']}, "
                        f"1q={res['one_qubit']}, 2q={res['two_qubit']}",
                        flush=True,
                    )

            device_cfg = toronto_cfg if device_key == "toronto" else brisbane_cfg
            qasmtrans_metrics = run_qasmtrans(
                qasmtrans_bin,
                qasm_path,
                device_cfg,
                device["name"],
                tmp_dir,
                max_time_ms=args.max_time_ms,
                retries=args.retries,
            )

            qt_time = qasmtrans_metrics["time_ms"]
            ratio = ""
            o1_time = qiskit_metrics[1]["time_ms"]
            if qt_time is not None and qt_time > 0 and o1_time:
                ratio = f"{o1_time / qt_time:.2f}"

            rows.append(
                {
                    "name": qasm_path.stem,
                    "circuit_file": str(qasm_path),
                    "logical_qubits": qubits,
                    "qiskit_o1_time_ms": f"{qiskit_metrics[1]['time_ms']:.3f}"
                    if qiskit_metrics[1]["time_ms"] is not None
                    else "",
                    "qiskit_o2_time_ms": f"{qiskit_metrics[2]['time_ms']:.3f}"
                    if qiskit_metrics[2]["time_ms"] is not None
                    else "",
                    "qiskit_o3_time_ms": f"{qiskit_metrics[3]['time_ms']:.3f}"
                    if qiskit_metrics[3]["time_ms"] is not None
                    else "",
                    "qmap_time_ms": "",
                    "qasmtrans_time_ms": f"{qt_time:.6f}" if qt_time is not None else "",
                    "ratio_o1_over_qt": ratio,
                    "qiskit_o1_single_qubit": qiskit_metrics[1]["one_qubit"] or "",
                    "qiskit_o2_single_qubit": qiskit_metrics[2]["one_qubit"] or "",
                    "qiskit_o3_single_qubit": qiskit_metrics[3]["one_qubit"] or "",
                    "qasmtrans_single_qubit": qasmtrans_metrics["one_qubit"]
                    if qasmtrans_metrics["one_qubit"] is not None
                    else "",
                    "qiskit_o1_two_qubit": qiskit_metrics[1]["two_qubit"] or "",
                    "qiskit_o2_two_qubit": qiskit_metrics[2]["two_qubit"] or "",
                    "qiskit_o3_two_qubit": qiskit_metrics[3]["two_qubit"] or "",
                    "qasmtrans_two_qubit": qasmtrans_metrics["two_qubit"]
                    if qasmtrans_metrics["two_qubit"] is not None
                    else "",
                    "qiskit_o1_depth": qiskit_metrics[1]["depth"] or "",
                    "qiskit_o2_depth": qiskit_metrics[2]["depth"] or "",
                    "qiskit_o3_depth": qiskit_metrics[3]["depth"] or "",
                    "qasmtrans_depth": qasmtrans_metrics["depth"]
                    if qasmtrans_metrics["depth"] is not None
                    else "",
                }
            )

    fieldnames = [
        "name",
        "circuit_file",
        "logical_qubits",
        "qiskit_o1_time_ms",
        "qiskit_o2_time_ms",
        "qiskit_o3_time_ms",
        "qmap_time_ms",
        "qasmtrans_time_ms",
        "ratio_o1_over_qt",
        "qiskit_o1_single_qubit",
        "qiskit_o2_single_qubit",
        "qiskit_o3_single_qubit",
        "qasmtrans_single_qubit",
        "qiskit_o1_two_qubit",
        "qiskit_o2_two_qubit",
        "qiskit_o3_two_qubit",
        "qasmtrans_two_qubit",
        "qiskit_o1_depth",
        "qiskit_o2_depth",
        "qiskit_o3_depth",
        "qasmtrans_depth",
    ]

    output_path = (repo_root / args.output_csv).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote combined results to {output_path}", flush=True)


if __name__ == "__main__":
    main()
