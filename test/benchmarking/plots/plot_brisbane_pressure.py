#!/usr/bin/env python3
"""Plot the >100-qubit Brisbane routing-pressure benchmark with error bars."""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, pstdev

import matplotlib.pyplot as plt


def load_rows(path: Path):
    with path.open("r", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    filtered = []
    for row in rows:
        if not row.get("depth_ratio_qt_over_qk"):
            continue
        row["pressure_bucket"] = int(row["pressure_bucket"])
        row["target_pressure"] = float(row["target_pressure"]) if row.get("target_pressure") else float(row["routing_pressure"])
        row["routing_pressure"] = float(row["routing_pressure"])
        row["depth_ratio_qt_over_qk"] = float(row["depth_ratio_qt_over_qk"])
        row["twoq_ratio_qt_over_qk"] = float(row["twoq_ratio_qt_over_qk"])
        row["time_ratio_qt_over_qk"] = float(row["time_ratio_qt_over_qk"])
        row["qasmtrans_swap_count"] = int(row["qasmtrans_swap_count"] or 0)
        row["qasmtrans_depth"] = int(row["qasmtrans_depth"])
        row["qiskit_depth"] = int(row["qiskit_depth"])
        filtered.append(row)
    return filtered


def bucket_stats(rows, key):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["target_pressure"]].append(row)
    stats = []
    for target_pressure in sorted(grouped):
        items = grouped[target_pressure]
        ys = [row[key] for row in items]
        stats.append(
            {
                "target_pressure": target_pressure,
                "x": target_pressure,
                "y": mean(ys),
                "yerr": pstdev(ys) if len(ys) > 1 else 0.0,
                "count": len(items),
                "items": items,
            }
        )
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Brisbane routing-pressure benchmark with error bars.")
    parser.add_argument("--input-csv", type=Path, default=Path("data/brisbane_routing_pressure_bench.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/brisbane_routing_pressure_errorbars.png"))
    parser.add_argument("--title-suffix", default="", help="Optional suffix appended to the figure title.")
    args = parser.parse_args()

    rows = load_rows(args.input_csv)
    depth_stats = bucket_stats(rows, "depth_ratio_qt_over_qk")
    twoq_stats = bucket_stats(rows, "twoq_ratio_qt_over_qk")
    time_stats = bucket_stats(rows, "time_ratio_qt_over_qk")
    swap_stats = bucket_stats(rows, "qasmtrans_swap_count")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    ax_depth, ax_twoq, ax_time, ax_swaps = axes.flat

    def draw(ax, stats, title, ylabel, baseline=None):
        if baseline is not None:
            ax.axhline(baseline, color="black", linestyle="--", linewidth=1)
        xs = [item["x"] for item in stats]
        ys = [item["y"] for item in stats]
        yerr = [item["yerr"] for item in stats]
        ax.errorbar(xs, ys, yerr=yerr, marker="o", capsize=4, linewidth=2)
        ax.set_title(title)
        ax.set_xlabel("Routing pressure (avg excess hops)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)

    for stat in depth_stats:
        jitter_scale = 0.03
        xs = [stat["x"] + jitter_scale * ((idx - (len(stat["items"]) - 1) / 2.0) / max(1.0, len(stat["items"]) / 2.0)) for idx, _ in enumerate(stat["items"])]
        ys = [item["depth_ratio_qt_over_qk"] for item in stat["items"]]
        ax_depth.scatter(xs, ys, alpha=0.25, color="tab:blue")
    for stat in twoq_stats:
        xs = [stat["x"] + 0.03 * ((idx - (len(stat["items"]) - 1) / 2.0) / max(1.0, len(stat["items"]) / 2.0)) for idx, _ in enumerate(stat["items"])]
        ys = [item["twoq_ratio_qt_over_qk"] for item in stat["items"]]
        ax_twoq.scatter(xs, ys, alpha=0.25, color="tab:orange")
    for stat in time_stats:
        xs = [stat["x"] + 0.03 * ((idx - (len(stat["items"]) - 1) / 2.0) / max(1.0, len(stat["items"]) / 2.0)) for idx, _ in enumerate(stat["items"])]
        ys = [item["time_ratio_qt_over_qk"] for item in stat["items"]]
        ax_time.scatter(xs, ys, alpha=0.25, color="tab:green")
    for stat in swap_stats:
        xs = [stat["x"] + 0.03 * ((idx - (len(stat["items"]) - 1) / 2.0) / max(1.0, len(stat["items"]) / 2.0)) for idx, _ in enumerate(stat["items"])]
        ys = [item["qasmtrans_swap_count"] for item in stat["items"]]
        ax_swaps.scatter(xs, ys, alpha=0.25, color="tab:red")

    def draw_errorbar(ax, stats, key, title, ylabel, baseline=None):
        if baseline is not None:
            ax.axhline(baseline, color="black", linestyle="--", linewidth=1)
        xs = [item["x"] for item in stats]
        ys = [item["y"] for item in stats]
        yerr = [item["yerr"] for item in stats]
        ax.errorbar(xs, ys, yerr=yerr, marker="o", capsize=4, linewidth=2, color="black")
        ax.set_title(title)
        ax.set_xlabel("Routing pressure (avg excess hops)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)

    draw_errorbar(ax_depth, depth_stats, "depth_ratio_qt_over_qk", "Depth Ratio: QASMTrans / Qiskit", "Ratio", baseline=1.0)
    draw_errorbar(ax_twoq, twoq_stats, "twoq_ratio_qt_over_qk", "2Q Gate Ratio: QASMTrans / Qiskit", "Ratio", baseline=1.0)
    draw_errorbar(ax_time, time_stats, "time_ratio_qt_over_qk", "Transpile Time Ratio: QASMTrans / Qiskit", "Ratio", baseline=1.0)
    draw_errorbar(ax_swaps, swap_stats, "qasmtrans_swap_count", "QASMTrans SWAP Count", "SWAPs")

    if rows:
        n = rows[0]["num_qubits"]
        input_depth = rows[0]["input_depth"]
        suffix = f" {args.title_suffix}" if args.title_suffix else ""
        fig.suptitle(
            f"Brisbane-Derived Routing-Pressure Benchmark{suffix} (n={n}, input depth={input_depth})",
            fontsize=14,
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)
    print(f"Wrote plot to {args.output}")


if __name__ == "__main__":
    main()
