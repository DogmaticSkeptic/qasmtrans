#!/usr/bin/env python3
"""Run the large-device Brisbane routing-pressure benchmark."""
from __future__ import annotations

import argparse
import csv
import json
import re
import signal
import subprocess
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
    edges = []
    for entry in cfg["cx_coupling"]:
        src, dst = map(int, str(entry).split("_"))
        edges.append((src, dst))
    return {
        "name": cfg.get("name") or cfg_path.stem,
        "map": CouplingMap(edges),
        "basis": list(dict.fromkeys(cfg.get("basis_gates") or ["rz", "sx", "x", "cx"])),
    }


def parse_swap_count(log_text: str) -> int | None:
    match = re.search(r"STEP-2\. Routing stats: .*swap=(\d+)", log_text)
    return int(match.group(1)) if match else None


def run_qasmtrans_case_inner(
    qasm_path: Path,
    device_path: Path,
    layout: List[int],
    optimize_1q: bool,
    optimize_2q_cancel: bool,
    optimize_commute_2q: bool,
) -> Dict[str, object]:
    import qasmtrans

    opts = qasmtrans.TranspileOptions()
    opts.mode = "ibmq"
    opts.backend_config = str(device_path)
    opts.disable_mapomatic = True
    opts.initial_layout = layout
    opts.optimize_1q = optimize_1q
    opts.optimize_2q_cancel = optimize_2q_cancel
    opts.optimize_commute_2q = optimize_commute_2q
    opts.verbose = 1

    start = time.perf_counter()
    res = qasmtrans.transpile_qasm(str(qasm_path), opts)
    elapsed_ms = (time.perf_counter() - start) * 1e3
    qc = QuantumCircuit.from_qasm_str(res.output_qasm)
    one, two = count_gate_types(qc)
    return {
        "status": "ok",
        "time_ms": elapsed_ms,
        "depth": qc.depth(),
        "one_qubit": one,
        "two_qubit": two,
        "swap_count": parse_swap_count(res.log),
        "error": "",
    }


def run_qasmtrans_case(
    qasm_path: Path,
    device_path: Path,
    layout: List[int],
    optimize_1q: bool,
    optimize_2q_cancel: bool,
    optimize_commute_2q: bool,
) -> Dict[str, object]:
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--child-qasmtrans",
        "--child-qasm-path",
        str(qasm_path),
        "--child-device-path",
        str(device_path),
        "--child-layout-json",
        json.dumps(layout),
    ]
    if optimize_1q:
        cmd.append("--optimize-1q")
    if optimize_2q_cancel:
        cmd.append("--optimize-2q-cancel")
    if optimize_commute_2q:
        cmd.append("--optimize-commute-2q")
    proc = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        return {
            "status": "error",
            "time_ms": None,
            "depth": None,
            "one_qubit": None,
            "two_qubit": None,
            "swap_count": None,
            "error": (proc.stderr or proc.stdout).strip(),
        }
    lines = [line for line in proc.stdout.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError("QASMTrans child process produced no JSON output")
    return json.loads(lines[-1])


def run_qiskit_case(
    qasm_path: Path,
    device: Dict[str, object],
    layout: List[int],
    optimization_level: int,
    seed: int,
    timeout_sec: int,
) -> Dict[str, object]:
    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
    class QiskitTimeout(RuntimeError):
        pass

    def handler(signum, frame):
        raise QiskitTimeout("Qiskit transpile timed out")

    old_handler = signal.signal(signal.SIGALRM, handler)
    start = time.perf_counter()
    try:
        signal.alarm(max(1, timeout_sec))
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
    except QiskitTimeout:
        elapsed_ms = (time.perf_counter() - start) * 1e3
        return {
            "time_ms": elapsed_ms,
            "depth": None,
            "one_qubit": None,
            "two_qubit": None,
            "status": "timeout",
        }
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
    elapsed_ms = (time.perf_counter() - start) * 1e3
    one, two = count_gate_types(tc)
    return {
        "time_ms": elapsed_ms,
        "depth": tc.depth(),
        "one_qubit": one,
        "two_qubit": two,
        "status": "ok",
    }


def iter_manifest_rows(path: Path) -> Iterable[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as fh:
        yield from csv.DictReader(fh)


def load_completed_case_ids(path: Path) -> set[str]:
    if not path.exists():
        return set()
    with path.open("r", newline="", encoding="utf-8") as fh:
        return {row["case_id"] for row in csv.DictReader(fh)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the Brisbane >100q routing-pressure benchmark.")
    parser.add_argument("--manifest", type=Path, default=Path("data/brisbane_routing_pressure_manifest.csv"))
    parser.add_argument("--output-csv", type=Path, default=Path("data/brisbane_routing_pressure_bench.csv"))
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--qiskit-level", type=int, default=3)
    parser.add_argument("--qiskit-timeout-sec", type=int, default=120)
    parser.add_argument("--optimize-1q", action="store_true")
    parser.add_argument("--optimize-2q-cancel", action="store_true")
    parser.add_argument("--optimize-commute-2q", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--child-qasmtrans", action="store_true")
    parser.add_argument("--child-qasm-path", type=Path)
    parser.add_argument("--child-device-path", type=Path)
    parser.add_argument("--child-layout-json")
    args = parser.parse_args()

    if args.child_qasmtrans:
        result = run_qasmtrans_case_inner(
            args.child_qasm_path,
            args.child_device_path,
            json.loads(args.child_layout_json),
            optimize_1q=args.optimize_1q,
            optimize_2q_cancel=args.optimize_2q_cancel,
            optimize_commute_2q=args.optimize_commute_2q,
        )
        print(json.dumps(result))
        return

    fieldnames = [
        "case_id",
        "circuit_variant",
        "num_qubits",
        "input_depth",
        "reps",
        "circuit_file",
        "device_file",
        "pressure_bucket",
        "target_pressure",
        "routing_pressure",
        "layout_tag",
        "initial_layout",
        "qasmtrans_status",
        "qasmtrans_time_ms",
        "qasmtrans_depth",
        "qasmtrans_1q",
        "qasmtrans_2q",
        "qasmtrans_swap_count",
        "qasmtrans_error",
        "qiskit_status",
        "qiskit_time_ms",
        "qiskit_depth",
        "qiskit_1q",
        "qiskit_2q",
        "depth_ratio_qt_over_qk",
        "twoq_ratio_qt_over_qk",
        "time_ratio_qt_over_qk",
    ]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    completed_case_ids = load_completed_case_ids(args.output_csv) if args.resume else set()
    mode = "a" if args.resume and args.output_csv.exists() else "w"
    with args.output_csv.open(mode, newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        if mode == "w":
            writer.writeheader()
            fh.flush()
        for row in iter_manifest_rows(args.manifest):
            if row["case_id"] in completed_case_ids:
                continue
            layout = json.loads(row["initial_layout"])
            qasm_path = Path(row["circuit_file"])
            device_path = Path(row["device_file"])
            device = load_device(device_path)
            print(
                f"[{row['case_id']}] variant={row['circuit_variant']} pressure={row['routing_pressure']} bucket={row['pressure_bucket']}",
                flush=True,
            )
            qt = run_qasmtrans_case(
                qasm_path,
                device_path,
                layout,
                optimize_1q=args.optimize_1q,
                optimize_2q_cancel=args.optimize_2q_cancel,
                optimize_commute_2q=args.optimize_commute_2q,
            )
            qk = run_qiskit_case(
                qasm_path,
                device,
                layout,
                optimization_level=args.qiskit_level,
                seed=args.seed,
                timeout_sec=args.qiskit_timeout_sec,
            )
            writer.writerow(
                {
                    **row,
                    "qasmtrans_status": qt["status"],
                    "qasmtrans_time_ms": f"{qt['time_ms']:.6f}" if qt["time_ms"] is not None else "",
                    "qasmtrans_depth": qt["depth"] if qt["depth"] is not None else "",
                    "qasmtrans_1q": qt["one_qubit"] if qt["one_qubit"] is not None else "",
                    "qasmtrans_2q": qt["two_qubit"] if qt["two_qubit"] is not None else "",
                    "qasmtrans_swap_count": qt["swap_count"] if qt["swap_count"] is not None else "",
                    "qasmtrans_error": qt["error"],
                    "qiskit_status": qk["status"],
                    "qiskit_time_ms": f"{qk['time_ms']:.6f}",
                    "qiskit_depth": qk["depth"] if qk["depth"] is not None else "",
                    "qiskit_1q": qk["one_qubit"] if qk["one_qubit"] is not None else "",
                    "qiskit_2q": qk["two_qubit"] if qk["two_qubit"] is not None else "",
                    "depth_ratio_qt_over_qk": f"{qt['depth'] / qk['depth']:.6f}" if qt["depth"] is not None and qk["depth"] else "",
                    "twoq_ratio_qt_over_qk": f"{qt['two_qubit'] / qk['two_qubit']:.6f}" if qt["two_qubit"] is not None and qk["two_qubit"] else "",
                    "time_ratio_qt_over_qk": f"{qt['time_ms'] / qk['time_ms']:.6f}" if qt["time_ms"] is not None and qk["time_ms"] else "",
                }
            )
            fh.flush()
    print(f"Wrote benchmark results to {args.output_csv}")


if __name__ == "__main__":
    main()
