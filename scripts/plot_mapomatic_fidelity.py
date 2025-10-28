#!/usr/bin/env python3
"""
Generate a bar plot comparing circuit fidelities across Mapomatic modes.

The script reads data/mapomatic_benchmarks.csv (produced by
scripts/run_mapomatic_bench.sh) and writes a publication-quality PDF to
data/mapomatic_fidelity.pdf.
"""

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
    raise SystemExit("matplotlib is required to generate the fidelity plot. "
                     "Install it with `pip install matplotlib`.") from exc


REPO_ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = REPO_ROOT / "data" / "mapomatic_benchmarks.csv"
PDF_PATH = REPO_ROOT / "data" / "mapomatic_fidelity.pdf"

MODE_ORDER = ["nomap", "product", "full"]
MODE_LABELS = OrderedDict(
    [
        ("nomap", "Baseline"),
        ("product", "Critical Path"),
        ("full", "Full Circuit"),
    ]
)
MODE_COLORS = {
    "nomap": "#6c757d",       # muted gray for baseline
    "product": "#1f77b4",     # blue
    "full": "#2ca02c",        # green
}
CIRCUIT_LABELS = {
    "vqe8": "VQE (8 qubits)",
    "qpe9": "QPE (9 qubits)",
    "sat7": "SAT (7 qubits)",
    "shor7": "Shor (7 qubits)",
}


def load_fidelity_data(csv_path: Path) -> Dict[str, Dict[str, float]]:
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing benchmark CSV at {csv_path}")

    data: Dict[str, Dict[str, float]] = defaultdict(dict)
    with csv_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            circuit = row["circuit"]
            mode = row["mode"]
            fidelity = float(row["fidelity"])
            data[circuit][mode] = fidelity
    return data


def plot_fidelities(data: Dict[str, Dict[str, float]], pdf_path: Path) -> None:
    circuits: List[str] = sorted(data.keys())
    num_modes = len(MODE_LABELS)
    x = range(len(circuits))
    width = 0.18

    fig, ax = plt.subplots(figsize=(7.5, 4.2))

    for idx, (mode, label) in enumerate(MODE_LABELS.items()):
        offsets = [pos + (idx - (num_modes - 1) / 2) * width for pos in x]
        fidelities = [data[circuit].get(mode, 0.0) for circuit in circuits]
        bars = ax.bar(
            offsets,
            fidelities,
            width=width,
            label=label,
            color=MODE_COLORS.get(mode, None),
        )
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + 0.01,
                f"{height:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax.set_ylabel("Simulated Fidelity", fontsize=11)
    ax.set_xticks(list(x))
    ax.set_xticklabels(
        [CIRCUIT_LABELS.get(circuit, circuit.upper()) for circuit in circuits],
        fontsize=11,
    )
    ax.set_ylim(0.0, 1.05)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.6)

    fig.tight_layout()
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    data = load_fidelity_data(CSV_PATH)
    plot_fidelities(data, PDF_PATH)
    print(f"Wrote fidelity comparison plot to {PDF_PATH}")


if __name__ == "__main__":
    main()
