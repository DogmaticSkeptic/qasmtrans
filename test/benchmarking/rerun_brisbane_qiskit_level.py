#!/usr/bin/env python3
"""Recompute the Qiskit side of a Brisbane routing-pressure benchmark CSV."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from run_brisbane_pressure_bench import load_device, run_qiskit_case


def main() -> None:
    parser = argparse.ArgumentParser(description="Reuse QASMTrans rows and rerun Qiskit at a different optimization level.")
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--qiskit-level", type=int, required=True)
    parser.add_argument("--qiskit-timeout-sec", type=int, default=120)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    with args.input_csv.open("r", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
        fieldnames = list(rows[0].keys())

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            print(
                f"[{row['case_id']}] pressure={row['target_pressure']} qiskit_level={args.qiskit_level}",
                flush=True,
            )
            device = load_device(Path(row["device_file"]))
            qk = run_qiskit_case(
                Path(row["circuit_file"]),
                device,
                layout=json.loads(row["initial_layout"]),
                optimization_level=args.qiskit_level,
                seed=args.seed,
                timeout_sec=args.qiskit_timeout_sec,
            )
            row["qiskit_status"] = qk["status"]
            row["qiskit_time_ms"] = f"{qk['time_ms']:.6f}"
            row["qiskit_depth"] = qk["depth"] if qk["depth"] is not None else ""
            row["qiskit_1q"] = qk["one_qubit"] if qk["one_qubit"] is not None else ""
            row["qiskit_2q"] = qk["two_qubit"] if qk["two_qubit"] is not None else ""

            qt_depth = int(row["qasmtrans_depth"]) if row.get("qasmtrans_depth") else None
            qt_twoq = int(row["qasmtrans_2q"]) if row.get("qasmtrans_2q") else None
            qt_time = float(row["qasmtrans_time_ms"]) if row.get("qasmtrans_time_ms") else None

            row["depth_ratio_qt_over_qk"] = f"{qt_depth / qk['depth']:.6f}" if qt_depth is not None and qk["depth"] else ""
            row["twoq_ratio_qt_over_qk"] = f"{qt_twoq / qk['two_qubit']:.6f}" if qt_twoq is not None and qk["two_qubit"] else ""
            row["time_ratio_qt_over_qk"] = f"{qt_time / qk['time_ms']:.6f}" if qt_time is not None and qk["time_ms"] else ""
            writer.writerow(row)
            fh.flush()

    print(f"Wrote updated benchmark results to {args.output_csv}")


if __name__ == "__main__":
    main()
