#!/usr/bin/env python3

"""
Utility to generate synthetic VQE-UCCSD style ansatz circuits for a range of
qubit counts and dump them as OpenQASM files.  The circuits are built using a
TwoLocal ansatz (Ry rotations with CZ entanglement) as a light‑weight proxy so
we can benchmark transpilation behaviour across different sizes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import TwoLocal
from qiskit.qasm2 import dumps as qasm2_dumps


def build_vqe_uccsd_ansatz(num_qubits: int, reps: int) -> QuantumCircuit:
    """Construct a simple layered ansatz mimicking VQE-UCCSD structure."""
    ansatz = TwoLocal(
        num_qubits=num_qubits,
        rotation_blocks="ry",
        entanglement_blocks="cz",
        entanglement="linear",
        reps=reps,
        insert_barriers=True,
    )
    qc = QuantumCircuit(num_qubits, name=f"vqe_uccsd_n{num_qubits}")
    qc.compose(ansatz, inplace=True)
    return qc


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate TwoLocal VQE-UCCSD proxy circuits as QASM."
    )
    parser.add_argument(
        "--min_qubits",
        type=int,
        default=10,
        help="Minimum number of qubits (inclusive). Default: 10",
    )
    parser.add_argument(
        "--max_qubits",
        type=int,
        default=120,
        help="Maximum number of qubits (inclusive). Default: 120",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=10,
        help="Qubit interval between generated circuits. Default: 10",
    )
    parser.add_argument(
        "--reps",
        type=int,
        default=6,
        help="Number of TwoLocal repetition layers. Default: 6",
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        default=Path("data/test_benchmark/vqe_uccsd_generated"),
        help="Directory to place generated QASM files. Default: data/test_benchmark/vqe_uccsd_generated",
    )
    args = parser.parse_args()

    if args.min_qubits < 1 or args.max_qubits < args.min_qubits:
        raise ValueError("Invalid min/max qubit range.")
    if args.step < 1:
        raise ValueError("Step must be >= 1.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for n in range(args.min_qubits, args.max_qubits + 1, args.step):
        qc = build_vqe_uccsd_ansatz(n, reps=args.reps)
        expanded = qc.decompose().decompose()  # expand TwoLocal instructions
        if expanded.parameters:
            param_map = {param: 0.0 for param in expanded.parameters}
            expanded = expanded.assign_parameters(param_map)
        qasm_path = args.output_dir / f"vqe_uccsd_n{n}.qasm"

        qasm_text = qasm2_dumps(expanded)
        qasm_path.write_text(qasm_text, encoding="utf-8")
        print(f"Wrote {qasm_path}")


if __name__ == "__main__":
    main()
