#!/usr/bin/env python3

"""Benchmark QASMTrans on generated VQE-UCCSD proxy circuits."""

from __future__ import annotations

import argparse
import re
import subprocess
import time
from pathlib import Path

import pandas as pd
from qiskit import QuantumCircuit


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


def run_qasmtrans(
    exe: Path,
    circuit: Path,
    device_json: Path,
    timeout_sec: float,
    attempts: int,
) -> dict[str, object]:
    tmp_dir = circuit.parent / "qasmtrans_outputs"
    tmp_dir.mkdir(exist_ok=True)

    for attempt in range(1, attempts + 1):
        out_qasm = tmp_dir / f"{circuit.stem}_attempt{attempt}.qasm"
        cmd = [
            str(exe),
            "-i",
            str(circuit),
            "-m",
            "ibmq",
            "-c",
            str(device_json),
            "-o",
            str(out_qasm),
            "-v",
            "1",
            "--disable_mapomatic",
        ]
        print(f"Running: {' '.join(cmd)}")
        start = time.perf_counter()
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout_sec,
                check=False,
            )
            elapsed_ms = (time.perf_counter() - start) * 1e3
        except subprocess.TimeoutExpired:
            print(f"  timeout after {timeout_sec}s")
            continue

        stdout = result.stdout or ""

        time_line = next(
            (line for line in stdout.splitlines() if line.startswith(" total QASMTrans time")),
            None,
        )
        if time_line:
            match = re.search(r"([-+]?\d*\.?\d+)\s*ms", time_line)
            reported_ms = float(match.group(1)) if match else float(elapsed_ms)
        else:
            reported_ms = float(elapsed_ms)

        if not out_qasm.exists():
            print("  output QASM missing, retrying")
            continue
        try:
            qc = QuantumCircuit.from_qasm_file(out_qasm)
            one_qubit, two_qubit = count_gate_types(qc)
            depth = qc.depth()
        except Exception as exc:
            print(f"  failed to parse QASM output: {exc}")
            continue

        return {
            "circuit": circuit.name,
            "reported_time_ms": reported_ms,
            "elapsed_time_ms": elapsed_ms,
            "one_qubit": one_qubit,
            "two_qubit": two_qubit,
            "depth": depth,
        }

    return {
        "circuit": circuit.name,
        "reported_time_ms": None,
        "elapsed_time_ms": None,
        "one_qubit": None,
        "two_qubit": None,
        "depth": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_dir",
        type=Path,
        default=Path("data/test_benchmark/vqe_uccsd_generated"),
        help="Directory containing generated VQE QASM files.",
    )
    parser.add_argument(
        "--device_json",
        type=Path,
        default=Path("data/devices/ibm_brisbane.json"),
        help="Device JSON for QASMTrans backend.",
    )
    parser.add_argument(
        "--qasmtrans_bin",
        type=Path,
        default=Path("build/QASMTrans"),
        help="Path to QASMTrans executable.",
    )
    parser.add_argument(
        "--timeout_ms",
        type=int,
        default=1000,
        help="Timeout per attempt in milliseconds.",
    )
    parser.add_argument(
        "--attempts",
        type=int,
        default=3,
        help="Retry count per circuit.",
    )
    parser.add_argument(
        "--output_csv",
        type=Path,
        default=Path("data/vqe_uccsd_qasmtrans_metrics.csv"),
        help="CSV destination for results.",
    )
    args = parser.parse_args()

    circuits = sorted(args.input_dir.glob("*.qasm"))
    if not circuits:
        raise FileNotFoundError(f"No QASM files found in {args.input_dir}")

    rows = []
    for circuit in circuits:
        print(f"=== {circuit.name} ===")
        row = run_qasmtrans(
            args.qasmtrans_bin,
            circuit,
            args.device_json,
            timeout_sec=args.timeout_ms / 1000.0,
            attempts=args.attempts,
        )
        rows.append(row)

    df = pd.DataFrame(rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_csv, index=False)
    print(f"Saved metrics to {args.output_csv}")


if __name__ == "__main__":
    main()
