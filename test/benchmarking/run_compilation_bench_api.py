#!/usr/bin/env python3
"""Benchmark Qiskit O1/O2/O3 vs QASMTrans using the Python API."""
from __future__ import annotations

import argparse
import csv
import json
import re
import tempfile
import time
import subprocess
import sys
import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from qiskit import QuantumCircuit, transpile
from qiskit.transpiler import CouplingMap

import qasmtrans

FIELDNAMES = [
    "name",
    "circuit_file",
    "logical_qubits",
    "qiskit_o1_time_ms",
    "qiskit_o2_time_ms",
    "qiskit_o3_time_ms",
    "qmap_time_ms",
    "qasmtrans_time_ms",
    "qasmtrans_reported_ms",
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


def count_gates_from_qasm(qasm_text: str) -> Tuple[int, int, Optional[int]]:
    one = two = 0
    depth = None
    lines = [ln.strip() for ln in qasm_text.splitlines() if ln.strip() and not ln.startswith("//")]
    qubit_depth: Dict[int, int] = {}
    for ln in lines:
        if ln.startswith(("OPENQASM", "include", "qreg", "creg")):
            continue
        if ln.startswith("measure"):
            continue
        parts = ln.rstrip(";").split()
        if len(parts) < 2:
            continue
        qubit_tokens = parts[1].split(",")
        qubits = []
        for tok in qubit_tokens:
            if "[" in tok and "]" in tok:
                try:
                    idx = int(tok[tok.find("[") + 1 : tok.find("]")])
                    qubits.append(idx)
                except ValueError:
                    continue
        if len(qubits) == 1:
            one += 1
        elif len(qubits) == 2:
            two += 1
        if qubits:
            max_depth = max(qubit_depth.get(q, 0) for q in qubits)
            new_depth = max_depth + 1
            for q in qubits:
                qubit_depth[q] = new_depth
    if qubit_depth:
        depth = max(qubit_depth.values())
    return one, two, depth


def count_qubits_from_qasm(qasm_text: str) -> int:
    qubits = 0
    for ln in qasm_text.splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("//"):
            continue
        if ln.startswith("qreg"):
            # example: qreg q[16];
            start = ln.find("[")
            end = ln.find("]")
            if start != -1 and end != -1 and end > start + 1:
                try:
                    qubits += int(ln[start + 1 : end])
                except ValueError:
                    continue
    return qubits


def run_qasmtrans_api(
    qasm_text: str,
    backend_cfg: Path,
    mode: str,
    max_time_ms: int,
    retries: int,
    disable_mapomatic: bool,
    output_dir: Path,
    log_file: Optional[Path] = None,
) -> Dict[str, Optional[float]]:
    log_lines: List[str] = []
    helper_code = """
import json, sys, time, os
import qasmtrans

qasm_text = sys.stdin.read()
opts = qasmtrans.TranspileOptions()
opts.backend_config = os.environ["QASMTRANS_BACKEND"]
opts.output_path = os.environ["QASMTRANS_OUTDIR"]
opts.mode = os.environ.get("QASMTRANS_MODE", "ibmq")
opts.disable_mapomatic = os.environ.get("QASMTRANS_DISABLE_MAPOMATIC", "0") == "1"
opts.verbose = int(os.environ.get("QASMTRANS_VERBOSE", "1"))
start = time.perf_counter()
res = qasmtrans.transpile_qasm(qasm_text, opts)
elapsed_ms = (time.perf_counter() - start) * 1e3
out = {
    "time_ms": elapsed_ms,
    "log": res.log,
}
print(json.dumps(out))
"""

    timeout_sec = max_time_ms / 1000.0
    for attempt in range(1, retries + 1):
        print(f"    QASMTrans attempt {attempt} (timeout {max_time_ms} ms)...", flush=True)
        log_lines.append(f"attempt {attempt}, timeout={max_time_ms} ms")
        try:
            proc = subprocess.run(
                [sys.executable, "-c", helper_code],
                input=qasm_text,
                text=True,
                capture_output=True,
                timeout=timeout_sec,
                env={
                    **os.environ,
                    "QASMTRANS_BACKEND": str(backend_cfg),
                    "QASMTRANS_OUTDIR": str(output_dir),
                    "QASMTRANS_MODE": mode,
                    "QASMTRANS_DISABLE_MAPOMATIC": "1" if disable_mapomatic else "0",
                    "QASMTRANS_VERBOSE": "1",
                },
            )
        except subprocess.TimeoutExpired:
            print("      QASMTrans timed out.", flush=True)
            log_lines.append("  timed out")
            if attempt == retries:
                if log_file:
                    log_file.write_text("\n".join(log_lines), encoding="utf-8")
                return {
                    "time_ms": None,
                    "reported_ms": None,
                    "one_qubit": None,
                    "two_qubit": None,
                    "depth": None,
                }
            continue

        if proc.returncode != 0:
            print(
                f"      QASMTrans failed (returncode {proc.returncode}). stdout:\n{proc.stdout}\n"
                f"stderr:\n{proc.stderr}",
                flush=True,
            )
            log_lines.append(f"  failed rc={proc.returncode}")
            if proc.stdout:
                log_lines.append("  stdout:\n" + proc.stdout)
            if proc.stderr:
                log_lines.append("  stderr:\n" + proc.stderr)
            if attempt == retries:
                if log_file:
                    log_file.write_text("\n".join(log_lines), encoding="utf-8")
                return {
                    "time_ms": None,
                    "reported_ms": None,
                    "one_qubit": None,
                    "two_qubit": None,
                    "depth": None,
                }
            continue

        try:
            stdout_text = proc.stdout or ""
            lines = [ln for ln in stdout_text.splitlines() if ln.strip()]
            json_line = lines[-1] if lines else ""
            payload = json.loads(json_line)
            elapsed_ms = float(payload.get("time_ms", 0.0))
            log = payload.get("log", "")
            # Prepend any stdout content before the JSON line into the log
            if len(lines) > 1:
                prefix = "\n".join(lines[:-1])
                log = (prefix + "\n" + str(log)) if log else prefix
        except Exception as ex:
            print(f"      QASMTrans output parse failed: {ex}. stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}", flush=True)
            log_lines.append(f"  parse failed: {ex}")
            if proc.stdout:
                log_lines.append("  stdout:\n" + proc.stdout)
            if proc.stderr:
                log_lines.append("  stderr:\n" + proc.stderr)
            if attempt == retries:
                if log_file:
                    log_file.write_text("\n".join(log_lines), encoding="utf-8")
                return {
                    "time_ms": None,
                    "reported_ms": None,
                    "one_qubit": None,
                    "two_qubit": None,
                    "depth": None,
                }
            continue

        one = two = depth = None
        try:
            from qasmtrans_core import count_gate_types as _dummy  # type: ignore
        except Exception:
            pass
        # Gate counting skipped to avoid large payloads; focus on timings.
        reported_ms = None
        if isinstance(log, str) and log:
            std_match = re.search(r"total QASMTrans time:\s*([-+]?\d*\.?\d+)", log)
            alt_match = re.search(r"\[timing\]\s*total_ms\s*=\s*([-+]?\d*\.?\d+)", log)
            match = std_match or alt_match
            if match:
                try:
                    reported_ms = float(match.group(1))
                except Exception:
                    reported_ms = None

        print(
            f"      QASMTrans success: wall={elapsed_ms:.3f} ms, "
            f"reported={'' if reported_ms is None else f'{reported_ms:.3f} ms'}, "
            f"1q={one}, 2q={two}, depth={depth}",
            flush=True,
        )
        log_lines.append(
            f"  success wall_ms={elapsed_ms:.3f}, reported_ms={reported_ms}, 1q={one}, 2q={two}, depth={depth}"
        )
        if isinstance(log, str) and log:
            log_lines.append("  qasmtrans log:\n" + log)
        if log_file:
            log_file.write_text("\n".join(log_lines), encoding="utf-8")
        return {
            "time_ms": elapsed_ms,
            "reported_ms": reported_ms,
            "one_qubit": one,
            "two_qubit": two,
            "depth": depth,
        }

    if log_file:
        log_file.write_text("\n".join(log_lines), encoding="utf-8")
    return {
        "time_ms": None,
        "reported_ms": None,
        "one_qubit": None,
        "two_qubit": None,
        "depth": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Qiskit + QASMTrans benchmarks via Python API.")
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
        "--mode",
        default="ibmq",
        help="QASMTrans mode (ibmq/ionq/quantinuum/rigetti/quafu).",
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
        "--disable_mapomatic",
        action="store_true",
        help="Disable Mapomatic in QASMTrans.",
    )
    parser.add_argument(
        "--output_csv",
        default=None,
        help="Destination CSV file (default: test/benchmarking/compilation_benchmarks.csv).",
    )
    parser.add_argument(
        "--log_dir",
        default=None,
        help="If set, write per-circuit QASMTrans logs to this directory.",
    )
    parser.add_argument(
        "--skip_names",
        default="",
        help="Comma-separated circuit basenames to skip (without .qasm).",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    repo_root = Path.cwd()
    qasm_dir = (repo_root / args.qasm_dir).resolve()
    if not qasm_dir.is_dir():
        raise FileNotFoundError(f"QASM directory not found: {qasm_dir}")

    toronto_cfg = (repo_root / args.toronto_config).resolve()
    brisbane_cfg = (repo_root / args.brisbane_config).resolve()

    devices = {
        "toronto": load_device(toronto_cfg),
        "brisbane": load_device(brisbane_cfg),
    }

    skip_set = {name.strip() for name in args.skip_names.split(",") if name.strip()}

    log_dir: Optional[Path] = None
    if args.log_dir:
        log_dir = (repo_root / args.log_dir).resolve()
    else:
        log_dir = script_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    if args.output_csv:
        output_path = (repo_root / args.output_csv).resolve()
    else:
        output_path = (script_dir / "compilation_benchmarks.csv").resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows_written = 0
    out_fh = output_path.open("w", newline="", encoding="utf-8")
    writer = csv.DictWriter(out_fh, fieldnames=FIELDNAMES)
    writer.writeheader()

    qasm_files = sorted(qasm_dir.glob("*.qasm"))
    print(f"Discovered {len(qasm_files)} circuit(s) under {qasm_dir}.", flush=True)
    if not qasm_files:
        out_fh.close()
        return

    total = len(qasm_files)

    with tempfile.TemporaryDirectory(prefix="qasmtrans_bench_") as tmpdir:
        tmp_path = Path(tmpdir)
        for idx, qasm_path in enumerate(qasm_files, start=1):
            if qasm_path.stem in skip_set:
                print(f"[{idx}/{total}] Skipping {qasm_path.name}: in skip list.", flush=True)
                continue

            qasm_text = qasm_path.read_text(encoding="utf-8")

            qc: Optional[QuantumCircuit] = None
            qubits: Optional[int] = None
            if args.skip_qiskit:
                qubits = count_qubits_from_qasm(qasm_text)
            else:
                try:
                    qc = QuantumCircuit.from_qasm_file(str(qasm_path))
                    qubits = qc.num_qubits
                except Exception as exc:
                    print(
                        f"[{idx}/{total}] Skipping {qasm_path.name}: failed to parse ({exc}).",
                        flush=True,
                    )
                    continue

            qubits = qubits or 0
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
                    if res["time_ms"] is None:
                        time_str = ""
                    else:
                        time_str = f"{res['time_ms']:.3f} ms"
                    print(
                        f"  Qiskit level {level}: {time_str}, "
                        f"depth={res['depth']}, 1q={res['one_qubit']}, 2q={res['two_qubit']}",
                        flush=True,
                    )

            device_cfg = toronto_cfg if device_key == "toronto" else brisbane_cfg
            circuit_log_path = None
            if log_dir:
                circuit_log_path = log_dir / f"{idx:03d}_{qasm_path.stem}.log"
            qt_metrics = run_qasmtrans_api(
                qasm_text,
                device_cfg,
                args.mode,
                max_time_ms=args.max_time_ms,
                retries=args.retries,
                disable_mapomatic=args.disable_mapomatic,
                output_dir=tmp_path,
                log_file=circuit_log_path,
            )

            qt_time = qt_metrics["time_ms"]
            ratio = ""
            o1_time = qiskit_metrics[1]["time_ms"]
            if qt_time is not None and qt_time > 0 and o1_time:
                ratio = f"{o1_time / qt_time:.2f}"

            row = {
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
                "qasmtrans_reported_ms": f"{qt_metrics['reported_ms']:.6f}"
                if qt_metrics.get("reported_ms") is not None
                else "",
                "ratio_o1_over_qt": ratio,
                "qiskit_o1_single_qubit": qiskit_metrics[1]["one_qubit"] or "",
                "qiskit_o2_single_qubit": qiskit_metrics[2]["one_qubit"] or "",
                "qiskit_o3_single_qubit": qiskit_metrics[3]["one_qubit"] or "",
                "qasmtrans_single_qubit": qt_metrics["one_qubit"] or "",
                "qiskit_o1_two_qubit": qiskit_metrics[1]["two_qubit"] or "",
                "qiskit_o2_two_qubit": qiskit_metrics[2]["two_qubit"] or "",
                "qiskit_o3_two_qubit": qiskit_metrics[3]["two_qubit"] or "",
                "qasmtrans_two_qubit": qt_metrics["two_qubit"] or "",
                "qiskit_o1_depth": qiskit_metrics[1]["depth"] or "",
                "qiskit_o2_depth": qiskit_metrics[2]["depth"] or "",
                "qiskit_o3_depth": qiskit_metrics[3]["depth"] or "",
                "qasmtrans_depth": qt_metrics["depth"] or "",
            }
            writer.writerow(row)
            out_fh.flush()
            rows_written += 1

    out_fh.close()
    print(f"\nWrote {rows_written} rows to {output_path}")


if __name__ == "__main__":
    main()
