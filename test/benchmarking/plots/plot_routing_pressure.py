#!/usr/bin/env python3
"""Plot fixed-layout routing-pressure benchmark results."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def load_rows(path: Path):
    with path.open("r", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        row["num_qubits"] = int(row["num_qubits"])
        row["routing_pressure"] = float(row["routing_pressure"])
        row["depth_ratio_qt_over_qk"] = float(row["depth_ratio_qt_over_qk"])
        row["twoq_ratio_qt_over_qk"] = float(row["twoq_ratio_qt_over_qk"])
        row["time_ratio_qt_over_qk"] = float(row["time_ratio_qt_over_qk"])
        row["qasmtrans_swap_count"] = int(row["qasmtrans_swap_count"] or 0)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot routing-pressure sweep metrics.")
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("data/routing_pressure_bench.csv"),
        help="CSV produced by run_routing_pressure_bench.py",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/routing_pressure_summary.png"),
        help="Where to save the summary plot.",
    )
    parser.add_argument(
        "--absolute-output",
        type=Path,
        default=Path("data/routing_pressure_absolute.png"),
        help="Where to save the absolute-deterioration plot.",
    )
    args = parser.parse_args()

    rows = load_rows(args.input_csv)
    series = defaultdict(list)
    for row in rows:
        series[row["num_qubits"]].append(row)
    for items in series.values():
        items.sort(key=lambda row: row["routing_pressure"])

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    ax_depth, ax_twoq, ax_time, ax_swap = axes.flat

    for num_qubits, items in sorted(series.items()):
        xs = [row["routing_pressure"] for row in items]
        ax_depth.plot(xs, [row["depth_ratio_qt_over_qk"] for row in items], marker="o", label=f"n={num_qubits}")
        ax_twoq.plot(xs, [row["twoq_ratio_qt_over_qk"] for row in items], marker="o", label=f"n={num_qubits}")
        ax_time.plot(xs, [row["time_ratio_qt_over_qk"] for row in items], marker="o", label=f"n={num_qubits}")
        ax_swap.plot(xs, [row["qasmtrans_swap_count"] for row in items], marker="o", label=f"n={num_qubits}")

    for ax in (ax_depth, ax_twoq, ax_time):
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
        ax.set_xlabel("Routing pressure (avg excess hops)")
        ax.grid(alpha=0.3)

    ax_swap.set_xlabel("Routing pressure (avg excess hops)")
    ax_swap.grid(alpha=0.3)

    ax_depth.set_title("Depth Ratio: QASMTrans / Qiskit")
    ax_depth.set_ylabel("Ratio")
    ax_twoq.set_title("2Q Gate Ratio: QASMTrans / Qiskit")
    ax_twoq.set_ylabel("Ratio")
    ax_time.set_title("Transpile Time Ratio: QASMTrans / Qiskit")
    ax_time.set_ylabel("Ratio")
    ax_swap.set_title("QASMTrans Routing SWAP Count")
    ax_swap.set_ylabel("SWAPs")

    handles, labels = ax_depth.get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=len(handles), frameon=False)

    fig.suptitle("Routing-Pressure Sweep on Fixed-Layout VQE-Style Circuits", fontsize=14)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)
    print(f"Wrote plot to {args.output}")

    absolute_fig, absolute_axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    ax_qt_depth, ax_qk_depth, ax_qt_twoq, ax_qk_twoq = absolute_axes.flat
    for num_qubits, items in sorted(series.items()):
        base = items[0]
        base_qt_depth = float(base["qasmtrans_depth"])
        base_qk_depth = float(base["qiskit_depth"])
        base_qt_twoq = float(base["qasmtrans_2q"])
        base_qk_twoq = float(base["qiskit_2q"])
        xs = [row["routing_pressure"] for row in items]
        ax_qt_depth.plot(
            xs,
            [float(row["qasmtrans_depth"]) / base_qt_depth for row in items],
            marker="o",
            label=f"n={num_qubits}",
        )
        ax_qk_depth.plot(
            xs,
            [float(row["qiskit_depth"]) / base_qk_depth for row in items],
            marker="o",
            label=f"n={num_qubits}",
        )
        ax_qt_twoq.plot(
            xs,
            [float(row["qasmtrans_2q"]) / base_qt_twoq for row in items],
            marker="o",
            label=f"n={num_qubits}",
        )
        ax_qk_twoq.plot(
            xs,
            [float(row["qiskit_2q"]) / base_qk_twoq for row in items],
            marker="o",
            label=f"n={num_qubits}",
        )

    for ax in (ax_qt_depth, ax_qk_depth, ax_qt_twoq, ax_qk_twoq):
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
        ax.set_xlabel("Routing pressure (avg excess hops)")
        ax.set_ylabel("Inflation vs pressure=0")
        ax.grid(alpha=0.3)

    ax_qt_depth.set_title("QASMTrans Depth Inflation")
    ax_qk_depth.set_title("Qiskit Depth Inflation")
    ax_qt_twoq.set_title("QASMTrans 2Q Inflation")
    ax_qk_twoq.set_title("Qiskit 2Q Inflation")

    handles, labels = ax_qt_depth.get_legend_handles_labels()
    if handles:
        absolute_fig.legend(handles, labels, loc="upper center", ncol=len(handles), frameon=False)

    absolute_fig.suptitle("Absolute Deterioration vs Routing Pressure", fontsize=14)
    args.absolute_output.parent.mkdir(parents=True, exist_ok=True)
    absolute_fig.savefig(args.absolute_output, dpi=200)
    print(f"Wrote plot to {args.absolute_output}")


if __name__ == "__main__":
    main()
