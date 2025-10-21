#!/usr/bin/env python3
"""Plot fidelity, latency, and pulse-count improvements from alg4_metrics.csv."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


COLOR_BASELINE = "#0072B2"  # Blue
COLOR_MERGED = "#009E73"  # Green
COLOR_CHANGE = "#E69F00"  # Orange


def read_metrics(csv_path: Path):
    rows = []
    with csv_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["baseline_fidelity"] = float(row["baseline_fidelity"])
            row["merged_fidelity"] = float(row["merged_fidelity"])
            row["baseline_latency_s"] = float(row["baseline_latency_s"])
            row["merged_latency_s"] = float(row["merged_latency_s"])
            row["baseline_pulses"] = int(row["baseline_pulses"])
            row["merged_pulses"] = int(row["merged_pulses"])
            rows.append(row)
    return rows


def annotate_bars(ax, bars, fmt="{:.3f}"):
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            fmt.format(height),
            ha="center",
            va="bottom",
            fontsize=8,
        )


def annotate_percent(ax, bars):
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{height:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )


def plot_fidelity(records, output_dir: Path):
    algorithms = [row["algorithm"] for row in records]
    baseline = [row["baseline_fidelity"] for row in records]
    merged = [row["merged_fidelity"] for row in records]

    positions = range(len(algorithms))
    width = 0.35

    fig, ax = plt.subplots(figsize=(5.0, 5.0), constrained_layout=True)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Fidelity")
    ax.set_title("Baseline vs. Merged Fidelity")

    bars_base = ax.bar(
        [p - width / 2 for p in positions],
        baseline,
        width,
        label="Baseline",
        color=COLOR_BASELINE,
    )
    bars_merged = ax.bar(
        [p + width / 2 for p in positions],
        merged,
        width,
        label="Merged",
        color=COLOR_MERGED,
    )
    annotate_bars(ax, bars_base)
    annotate_bars(ax, bars_merged)

    ax.set_xticks(list(positions))
    ax.set_xticklabels(algorithms)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    out_path = output_dir / "alg4_fidelity_comparison.pdf"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote {out_path}")


def plot_percent_change(values_base, values_merged, algorithms, ylabel, title, filename):
    reductions = []
    for base, merged in zip(values_base, values_merged):
        if base == 0:
            reductions.append(0.0)
        else:
            reductions.append((base - merged) / base * 100.0)

    positions = range(len(algorithms))

    fig, ax = plt.subplots(figsize=(5.0, 5.0), constrained_layout=True)
    bars = ax.bar(
        positions,
        reductions,
        color=COLOR_CHANGE,
    )
    annotate_percent(ax, bars)

    ax.set_xticks(list(positions))
    ax.set_xticklabels(algorithms)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.yaxis.set_major_formatter(PercentFormatter(decimals=0))
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    out_path = filename
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    repo_root = Path(__file__).resolve().parents[1]
    csv_path = repo_root / "data" / "output" / "alg4" / "alg4_metrics.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Metrics CSV not found at {csv_path}")

    output_dir = csv_path.parent
    records = read_metrics(csv_path)
    algorithms = [row["algorithm"] for row in records]

    plot_fidelity(records, output_dir)

    baseline_latency = [row["baseline_latency_s"] for row in records]
    merged_latency = [row["merged_latency_s"] for row in records]
    plot_percent_change(
        baseline_latency,
        merged_latency,
        algorithms,
        "Latency Reduction (%)",
        "Critical Path Reduction",
        output_dir / "alg4_latency_reduction.pdf",
    )

    baseline_pulses = [row["baseline_pulses"] for row in records]
    merged_pulses = [row["merged_pulses"] for row in records]
    plot_percent_change(
        baseline_pulses,
        merged_pulses,
        algorithms,
        "Pulse Reduction (%)",
        "Pulse Count Reduction",
        output_dir / "alg4_pulse_reduction.pdf",
    )


if __name__ == "__main__":
    main()
