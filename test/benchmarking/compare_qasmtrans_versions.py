#!/usr/bin/env python3
"""Benchmark old vs new QASMTrans on the NWQBench-style corpus."""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional

from qiskit import QuantumCircuit

QREG_RE = re.compile(r"^\s*qreg\s+([A-Za-z_]\w*)\[(\d+)\]\s*;")
REF_RE = re.compile(r"([A-Za-z_]\w*)\[(\d+)\]")


def count_gate_types(qc: QuantumCircuit) -> tuple[int, int]:
    one = two = 0
    for instruction in qc.data:
        inst = instruction.operation
        qargs = instruction.qubits
        if inst.name in {"barrier", "measure", "delay"}:
            continue
        qubit_count = len(qargs)
        if qubit_count == 1:
            one += 1
        elif qubit_count == 2:
            two += 1
    return one, two


def choose_device_config(repo_root: Path, qubits: int) -> tuple[str, Path]:
    if qubits <= 27:
        return "ibmq_toronto", (repo_root / "data/devices/ibmq_toronto.json").resolve()
    if qubits <= 127:
        return "ibm_brisbane", (repo_root / "data/devices/ibm_brisbane.json").resolve()
    raise ValueError(f"Unsupported qubit count {qubits}: exceeds 127-qubit comparison target")


def parse_emitted_qasm_metrics(qasm_path: Path) -> Dict[str, int]:
    qreg_offsets: Dict[str, int] = {}
    total_qubits = 0
    one = two = depth = 0
    depths: List[int] = []
    with qasm_path.open("r") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line or line.startswith("//"):
                continue
            qreg_match = QREG_RE.match(line)
            if qreg_match:
                name = qreg_match.group(1)
                size = int(qreg_match.group(2))
                qreg_offsets[name] = total_qubits
                total_qubits += size
                depths.extend([0] * size)
                continue
            op_name = line.split(None, 1)[0].lower() if line.split(None, 1) else ""
            if op_name in {"openqasm", "include", "creg", "gate", "opaque", "measure", "barrier", "delay"}:
                continue

            qubits: List[int] = []
            for reg_name, index_text in REF_RE.findall(line):
                if reg_name not in qreg_offsets:
                    continue
                qubits.append(qreg_offsets[reg_name] + int(index_text))
            if not qubits:
                continue
            if len(qubits) == 1:
                one += 1
            elif len(qubits) == 2:
                two += 1
            current_depth = max(depths[q] for q in qubits) + 1
            for q in qubits:
                depths[q] = current_depth
            if current_depth > depth:
                depth = current_depth
    return {"one_qubit": one, "two_qubit": two, "depth": depth}


def run_binary(
    exe: Path,
    circuit: Path,
    device_cfg: Path,
    output_qasm: Path,
    timeout_ms: int,
    retries: int,
    extra_args: List[str],
) -> Dict[str, Optional[float]]:
    timeout_sec = timeout_ms / 1000.0
    last_status = "failed"
    for attempt in range(1, retries + 1):
        if output_qasm.exists():
            output_qasm.unlink()
        cmd = [
            str(exe),
            "-i",
            str(circuit),
            "-m",
            "ibmq",
            "-c",
            str(device_cfg),
            "-o",
            str(output_qasm),
            "-v",
            "1",
        ] + list(extra_args)
        start = time.perf_counter()
        try:
            result = subprocess.run(
                cmd,
                text=True,
                capture_output=True,
                timeout=timeout_sec,
                check=False,
            )
        except subprocess.TimeoutExpired:
            last_status = "timeout"
            continue
        duration_ms = (time.perf_counter() - start) * 1e3
        if result.returncode != 0:
            last_status = f"exit_{result.returncode}"
            continue
        if not output_qasm.exists():
            last_status = "missing_output"
            continue
        try:
            metrics = parse_emitted_qasm_metrics(output_qasm)
        except Exception:
            last_status = "parse_error"
            continue
        return {
            "status": "ok",
            "attempts": attempt,
            "time_ms": duration_ms,
            "one_qubit": metrics["one_qubit"],
            "two_qubit": metrics["two_qubit"],
            "depth": metrics["depth"],
        }
    return {
        "status": last_status,
        "attempts": retries,
        "time_ms": None,
        "one_qubit": None,
        "two_qubit": None,
        "depth": None,
    }


def format_optional_float(value: Optional[float], digits: int = 3) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare old vs new QASMTrans outputs and runtime.")
    parser.add_argument("--qasm_dir", default="data/test_benchmark", help="Directory containing *.qasm circuits.")
    parser.add_argument("--new_bin", default="build/QASMTrans", help="Path to the new QASMTrans binary.")
    parser.add_argument(
        "--old_bin",
        default="/home/hoyt924.linux/qasmtransoriginal/qasmtrans/build/QASMTrans",
        help="Path to the old QASMTrans binary.",
    )
    parser.add_argument(
        "--max_qubits",
        type=int,
        default=127,
        help="Skip circuits above this many qubits.",
    )
    parser.add_argument(
        "--timeout_ms",
        type=int,
        default=5000,
        help="Per-attempt timeout for each QASMTrans binary.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=2,
        help="Number of attempts per binary after timeout/failure.",
    )
    parser.add_argument(
        "--output_csv",
        default="data/qasmtrans_old_vs_new_nwqbench.csv",
        help="Destination CSV.",
    )
    parser.add_argument(
        "--small_device_cfg",
        default="",
        help="Optional backend JSON override for circuits with <=27 qubits.",
    )
    parser.add_argument(
        "--large_device_cfg",
        default="",
        help="Optional backend JSON override for circuits with 28-127 qubits.",
    )
    parser.add_argument(
        "--new-arg",
        action="append",
        default=[],
        help="Extra argument to pass to the new QASMTrans binary. Repeat for multiple flags.",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    qasm_dir = (repo_root / args.qasm_dir).resolve()
    new_bin = (repo_root / args.new_bin).resolve()
    old_bin = Path(args.old_bin).resolve()
    output_csv = (repo_root / args.output_csv).resolve()
    small_device_override = Path(args.small_device_cfg).resolve() if args.small_device_cfg else None
    large_device_override = Path(args.large_device_cfg).resolve() if args.large_device_cfg else None
    if not qasm_dir.is_dir():
        raise FileNotFoundError(f"Missing QASM directory: {qasm_dir}")
    if not new_bin.is_file():
        raise FileNotFoundError(f"Missing new binary: {new_bin}")
    if not old_bin.is_file():
        raise FileNotFoundError(f"Missing old binary: {old_bin}")
    if small_device_override is not None and not small_device_override.is_file():
        raise FileNotFoundError(f"Missing small-device override: {small_device_override}")
    if large_device_override is not None and not large_device_override.is_file():
        raise FileNotFoundError(f"Missing large-device override: {large_device_override}")

    rows: List[Dict[str, object]] = []
    qasm_files = sorted(qasm_dir.glob("*.qasm"))
    tmp_parent = repo_root / "tmp"
    tmp_parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="qasmtrans_compare_", dir=tmp_parent) as tmp:
        tmp_dir = Path(tmp)
        total = len(qasm_files)
        for idx, qasm_path in enumerate(qasm_files, start=1):
            try:
                source_qc = QuantumCircuit.from_qasm_file(str(qasm_path))
            except Exception as exc:
                print(f"[{idx}/{total}] Skipping {qasm_path.name}: parse failed ({exc})", flush=True)
                rows.append(
                    {
                        "name": qasm_path.stem,
                        "circuit_file": str(qasm_path),
                        "logical_qubits": "",
                        "device_name": "",
                        "status": "source_parse_failed",
                        "old_status": "",
                        "new_status": "",
                    }
                )
                continue

            qubits = source_qc.num_qubits
            if qubits > args.max_qubits:
                print(f"[{idx}/{total}] Skipping {qasm_path.name}: {qubits} qubits (> {args.max_qubits})", flush=True)
                rows.append(
                    {
                        "name": qasm_path.stem,
                        "circuit_file": str(qasm_path),
                        "logical_qubits": qubits,
                        "device_name": "",
                        "status": "skipped_too_wide",
                        "old_status": "",
                        "new_status": "",
                    }
                )
                continue

            device_name, device_cfg = choose_device_config(repo_root, qubits)
            if qubits <= 27 and small_device_override is not None:
                device_name = small_device_override.stem
                device_cfg = small_device_override
            elif qubits > 27 and large_device_override is not None:
                device_name = large_device_override.stem
                device_cfg = large_device_override
            print(f"[{idx}/{total}] {qasm_path.name}: {qubits} qubits -> {device_name}", flush=True)
            old_output = tmp_dir / f"{qasm_path.stem}_old.qasm"
            new_output = tmp_dir / f"{qasm_path.stem}_new.qasm"
            old_result = run_binary(
                old_bin,
                qasm_path,
                device_cfg,
                old_output,
                timeout_ms=args.timeout_ms,
                retries=args.retries,
                extra_args=[],
            )
            new_result = run_binary(
                new_bin,
                qasm_path,
                device_cfg,
                new_output,
                timeout_ms=args.timeout_ms,
                retries=args.retries,
                extra_args=args.new_arg,
            )

            status = "completed" if old_result["status"] == "ok" and new_result["status"] == "ok" else "partial_failure"
            old_time = old_result["time_ms"]
            new_time = new_result["time_ms"]
            speedup = ""
            if old_time is not None and new_time is not None and new_time > 0:
                speedup = f"{old_time / new_time:.2f}"
            rows.append(
                {
                    "name": qasm_path.stem,
                    "circuit_file": str(qasm_path),
                    "logical_qubits": qubits,
                    "device_name": device_name,
                    "status": status,
                    "old_status": old_result["status"],
                    "new_status": new_result["status"],
                    "old_attempts": old_result["attempts"],
                    "new_attempts": new_result["attempts"],
                    "old_time_ms": format_optional_float(old_time, 6),
                    "new_time_ms": format_optional_float(new_time, 6),
                    "old_over_new_time_ratio": speedup,
                    "old_1q": old_result["one_qubit"] if old_result["one_qubit"] is not None else "",
                    "new_1q": new_result["one_qubit"] if new_result["one_qubit"] is not None else "",
                    "delta_1q_new_minus_old": (
                        new_result["one_qubit"] - old_result["one_qubit"]
                        if old_result["one_qubit"] is not None and new_result["one_qubit"] is not None
                        else ""
                    ),
                    "old_2q": old_result["two_qubit"] if old_result["two_qubit"] is not None else "",
                    "new_2q": new_result["two_qubit"] if new_result["two_qubit"] is not None else "",
                    "delta_2q_new_minus_old": (
                        new_result["two_qubit"] - old_result["two_qubit"]
                        if old_result["two_qubit"] is not None and new_result["two_qubit"] is not None
                        else ""
                    ),
                    "old_depth": old_result["depth"] if old_result["depth"] is not None else "",
                    "new_depth": new_result["depth"] if new_result["depth"] is not None else "",
                    "delta_depth_new_minus_old": (
                        new_result["depth"] - old_result["depth"]
                        if old_result["depth"] is not None and new_result["depth"] is not None
                        else ""
                    ),
                }
            )
            print(
                "  old="
                f"{old_result['status']} time={format_optional_float(old_time)}ms "
                f"1q={old_result['one_qubit']} 2q={old_result['two_qubit']} depth={old_result['depth']}; "
                "new="
                f"{new_result['status']} time={format_optional_float(new_time)}ms "
                f"1q={new_result['one_qubit']} 2q={new_result['two_qubit']} depth={new_result['depth']}",
                flush=True,
            )

    fieldnames = [
        "name",
        "circuit_file",
        "logical_qubits",
        "device_name",
        "status",
        "old_status",
        "new_status",
        "old_attempts",
        "new_attempts",
        "old_time_ms",
        "new_time_ms",
        "old_over_new_time_ratio",
        "old_1q",
        "new_1q",
        "delta_1q_new_minus_old",
        "old_2q",
        "new_2q",
        "delta_2q_new_minus_old",
        "old_depth",
        "new_depth",
        "delta_depth_new_minus_old",
    ]
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote results to {output_csv}", flush=True)


if __name__ == "__main__":
    main()
