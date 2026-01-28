#!/usr/bin/env python3
"""Plot Mapomatic timing results from data/mapomatic_benchmarks.csv."""

from __future__ import annotations

import csv
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Dict, List

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise SystemExit(
        "matplotlib is required to create the timing plot. "
        "Install it with `pip install matplotlib`."
    ) from exc


REPO_ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = REPO_ROOT / "data" / "mapomatic_benchmarks.csv"
PDF_PATH = REPO_ROOT / "data" / "mapomatic_timing.pdf"

MODE_ORDER = ["product"]
MODE_LABELS = OrderedDict(
    [
        ("product", "Critical Path"),
    ]
)
MODE_COLORS = {
    "product": "#1f77b4",
}
CIRCUIT_LABELS = {
    "vqe8": "VQE (8 qubits)",
    "qpe9": "QPE (9 qubits)",
    "sat7": "SAT (7 qubits)",
    "shor7": "Shor (7 qubits)",
}


def load_timing_data(csv_path: Path) -> Dict[str, Dict[str, float]]:
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing benchmark CSV at {csv_path}")

    data: Dict[str, Dict[str, float]] = defaultdict(dict)
    with csv_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            circuit = row["circuit"]
            mode = row["mode"]
            time_ms = float(row["mapomatic_ms"]) if row["mapomatic_ms"] else 0.0
            data[circuit][mode] = time_ms
    return data


def compute_percentages(data: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    percentages: Dict[str, Dict[str, float]] = defaultdict(dict)
    for circuit, timings in data.items():
        full_time = timings.get("full", 0.0)
        if full_time <= 0.0:
            continue
        for mode in MODE_ORDER:
            mode_time = timings.get(mode, 0.0)
            if mode_time > 0.0:
                percentages[circuit][mode] = (mode_time / full_time) * 100.0
            else:
                percentages[circuit][mode] = 0.0
    return percentages


def plot_percentages(percentages: Dict[str, Dict[str, float]], pdf_path: Path) -> None:
    circuits: List[str] = sorted(percentages.keys())
    x = range(len(circuits))
    width = 0.25

    fig, ax = plt.subplots(figsize=(6.5, 4.0))

    for idx, (mode, label) in enumerate(MODE_LABELS.items()):
        offsets = [pos + (idx - (len(MODE_LABELS) - 1) / 2) * width for pos in x]
        values = [percentages[circuit].get(mode, 0.0) for circuit in circuits]
        bars = ax.bar(
            offsets,
            values,
            width=width,
            color=MODE_COLORS.get(mode, None),
        )
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + 2,
                f"{height:.0f}%",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax.set_ylabel("Critical Path Time vs Full Fidelity", fontsize=11)
    ax.set_xticks(list(x))
    ax.set_xticklabels(
        [CIRCUIT_LABELS.get(circuit, circuit.upper()) for circuit in circuits],
        fontsize=11,
    )
    max_val = max((max(vals.values()) for vals in percentages.values()), default=0.0)
    ax.set_ylim(0.0, max(110.0, max_val * 1.1))
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.6)

    fig.tight_layout()
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    timing_data = load_timing_data(CSV_PATH)
    percentages = compute_percentages(timing_data)
    if not percentages:
        raise SystemExit("No valid timing data found in CSV.")
    plot_percentages(percentages, PDF_PATH)
    print(f"Wrote timing comparison plot to {PDF_PATH}")


if __name__ == "__main__":
    main()
