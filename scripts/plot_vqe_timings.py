#!/usr/bin/env python3
"""Plot VQE-UCCSD transpilation times: QASMTrans vs. Qiskit O3."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import csv
import matplotlib.pyplot as plt


VQE_PATTERN = re.compile(r"vqe_uccsd_n(?P<qubits>\d+)\.qasm", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a QASMTrans vs. Qiskit O3 timing plot for VQE-UCCSD circuits."
    )
    parser.add_argument(
        "--input_csv",
        default="data/vqe_uccsd_qasmtrans_metrics.csv",
        help="CSV file produced by run_qasmtrans_vqe.sh.",
    )
    parser.add_argument(
        "--output_pdf",
        default="figures/vqe_uccsd_qasmtrans_vs_qiskit_o3.pdf",
        help="Destination PDF path for the plot.",
    )
    parser.add_argument(
        "--min_qubits",
        type=int,
        default=10,
        help="Minimum logical qubit count to include (inclusive).",
    )
    parser.add_argument(
        "--max_qubits",
        type=int,
        default=100,
        help="Maximum logical qubit count to include (inclusive).",
    )
    return parser.parse_args()


def extract_qubits(circuit_name: str) -> int | None:
    match = VQE_PATTERN.search(circuit_name)
    return int(match.group("qubits")) if match else None


def load_and_filter(csv_path: Path, min_qubits: int, max_qubits: int) -> list[dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    with csv_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            qubits = extract_qubits(raw.get("circuit", ""))
            if qubits is None:
                continue
            if qubits < min_qubits or qubits > max_qubits:
                continue

            try:
                qt_time = float(raw["qasmtrans_reported_ms"])
                qiskit_time = float(raw["qiskit_o3_time_ms"])
            except (KeyError, TypeError, ValueError):
                continue

            rows[qubits] = {
                "qubits": qubits,
                "qasmtrans_time": qt_time,
                "qiskit_o3_time": qiskit_time,
            }

    data = [rows[key] for key in sorted(rows)]
    if not data:
        raise ValueError("No VQE-UCCSD rows found in the requested qubit range.")
    return data


def make_plot(data: list[dict[str, float]], output_pdf: Path) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(6, 6))

    qubits = [row["qubits"] for row in data]
    qasmtrans_times = [row["qasmtrans_time"] for row in data]
    qiskit_times = [row["qiskit_o3_time"] for row in data]

    ax.plot(
        qubits,
        qasmtrans_times,
        marker="o",
        linestyle="-",
        linewidth=2,
        markersize=6,
        color="#1f77b4",
        label="QASMTrans",
    )
    ax.plot(
        qubits,
        qiskit_times,
        marker="s",
        linestyle="-",
        linewidth=2,
        markersize=6,
        color="#d62728",
        label="Qiskit (O3)",
    )

    label_kwargs = {"fontsize": 18}
    tick_fontsize = 14

    ax.set_xlabel("Qubit #", **label_kwargs)
    ax.set_ylabel("Transpilation time (ms)", **label_kwargs)
    ax.set_xlim(min(qubits) - 2, max(qubits) + 2)
    ax.legend(frameon=True, fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=tick_fontsize)
    ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.7)

    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_pdf, format="pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    csv_path = Path(args.input_csv)
    if not csv_path.is_file():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")

    data = load_and_filter(csv_path, args.min_qubits, args.max_qubits)

    output_pdf = Path(args.output_pdf)
    make_plot(data, output_pdf)
    print(f"Saved plot to {output_pdf}")


if __name__ == "__main__":
    main()
