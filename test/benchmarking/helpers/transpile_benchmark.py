#!/usr/bin/env python3
from __future__ import annotations

import time
from pathlib import Path

from qiskit import QuantumCircuit, transpile


def count_1q_2q(circ: QuantumCircuit) -> tuple[int, int]:
    one_q = 0
    two_q = 0
    for inst, qargs, _cargs in circ.data:
        # Skip barriers to keep counts comparable
        if inst.name == "barrier":
            continue
        nq = inst.num_qubits
        if nq == 1:
            one_q += 1
        elif nq == 2:
            two_q += 1
    return one_q, two_q


def main() -> None:
    base = Path(__file__).resolve().parents[1] / "data" / "test_benchmark"
    qasm_files = sorted(p for p in base.rglob("*.qasm") if p.is_file())
    if not qasm_files:
        print(f"No .qasm files found under {base}")
        return

    print("file,opt_level,time_ms,depth,one_q,two_q")

    for qasm_path in qasm_files:
        try:
            circ = QuantumCircuit.from_qasm_file(str(qasm_path))
        except Exception as exc:  # noqa: BLE001
            print(f"{qasm_path.name},ERR,0,0,0,0  # failed to load: {exc}")
            continue

        for level in (1, 2, 3):
            start = time.perf_counter()
            try:
                tcirc = transpile(circ, optimization_level=level)
            except Exception as exc:  # noqa: BLE001
                elapsed_ms = (time.perf_counter() - start) * 1000.0
                print(f"{qasm_path.name},{level},{elapsed_ms:.3f},0,0,0  # transpile failed: {exc}")
                continue
            elapsed_ms = (time.perf_counter() - start) * 1000.0
            depth = tcirc.depth()
            one_q, two_q = count_1q_2q(tcirc)
            print(f"{qasm_path.name},{level},{elapsed_ms:.3f},{depth},{one_q},{two_q}")


if __name__ == "__main__":
    main()
