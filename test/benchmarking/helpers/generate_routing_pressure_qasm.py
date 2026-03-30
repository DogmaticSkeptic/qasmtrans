#!/usr/bin/env python3
"""Generate topology-aware VQE-style circuits plus layout-stress manifests."""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import Counter
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from qiskit import QuantumCircuit
from qiskit.circuit.library import TwoLocal
from qiskit.qasm2 import dumps as qasm2_dumps


def build_vqe_proxy(num_qubits: int, reps: int) -> QuantumCircuit:
    ansatz = TwoLocal(
        num_qubits=num_qubits,
        rotation_blocks="ry",
        entanglement_blocks="cz",
        entanglement="linear",
        reps=reps,
        insert_barriers=False,
    )
    qc = QuantumCircuit(num_qubits, name=f"vqe_layoutstress_n{num_qubits}")
    qc.compose(ansatz, inplace=True)
    expanded = qc.decompose().decompose()
    if expanded.parameters:
        # Keep non-zero angles so post-routing cleanups do not trivially erase the ansatz.
        param_map = {
            param: 0.113 * (idx + 1)
            for idx, param in enumerate(sorted(expanded.parameters, key=lambda p: p.name))
        }
        expanded = expanded.assign_parameters(param_map)
    return expanded


def count_weighted_interactions(qc: QuantumCircuit) -> Counter[Tuple[int, int]]:
    weights: Counter[Tuple[int, int]] = Counter()
    for inst, qargs, _ in qc.data:
        if inst.name in {"barrier", "measure", "delay", "reset"}:
            continue
        if len(qargs) != 2:
            continue
        q0 = qc.find_bit(qargs[0]).index
        q1 = qc.find_bit(qargs[1]).index
        edge = (q0, q1) if q0 < q1 else (q1, q0)
        weights[edge] += 1
    return weights


def layout_pressure(weights: Counter[Tuple[int, int]], layout: Sequence[int]) -> Tuple[float, int]:
    total_weight = sum(weights.values())
    if total_weight == 0:
        return 0.0, 1
    weighted_excess = 0.0
    max_hop = 1
    for (u, v), w in weights.items():
        hop = abs(layout[u] - layout[v])
        max_hop = max(max_hop, hop)
        weighted_excess += w * max(0, hop - 1)
    return weighted_excess / total_weight, max_hop


def random_layouts(num_qubits: int, seed: int, samples: int) -> Iterable[Tuple[str, List[int]]]:
    rng = random.Random(seed)
    identity = list(range(num_qubits))
    yield ("identity", identity.copy())
    yield ("reverse", list(reversed(identity)))
    weave = []
    lo = 0
    hi = num_qubits - 1
    while lo <= hi:
        weave.append(lo)
        if lo != hi:
            weave.append(hi)
        lo += 1
        hi -= 1
    yield ("weave", weave)

    for idx in range(samples):
        layout = identity.copy()
        swap_budget = max(1, math.ceil((idx + 1) * num_qubits / 3))
        for _ in range(swap_budget):
            i = rng.randrange(num_qubits)
            j = rng.randrange(num_qubits)
            layout[i], layout[j] = layout[j], layout[i]
        yield (f"random_{idx:03d}", layout)


def select_pressure_levels(
    weights: Counter[Tuple[int, int]],
    num_qubits: int,
    levels: int,
    seed: int,
    samples: int,
) -> List[Tuple[str, List[int], float, int]]:
    candidates = []
    seen = set()
    for tag, layout in random_layouts(num_qubits, seed=seed, samples=samples):
        key = tuple(layout)
        if key in seen:
            continue
        seen.add(key)
        pressure, max_hop = layout_pressure(weights, layout)
        candidates.append((tag, layout, pressure, max_hop))

    candidates.sort(key=lambda item: (item[2], item[3], item[0]))
    if not candidates:
        raise RuntimeError("No candidate layouts generated")

    chosen = []
    used = set()
    target_count = min(levels, len(candidates))
    for level in range(target_count):
        idx = round(level * (len(candidates) - 1) / max(1, target_count - 1))
        while idx < len(candidates) and tuple(candidates[idx][1]) in used:
            idx += 1
        if idx >= len(candidates):
            idx = len(candidates) - 1
            while idx >= 0 and tuple(candidates[idx][1]) in used:
                idx -= 1
        if idx < 0:
            break
        chosen.append(candidates[idx])
        used.add(tuple(candidates[idx][1]))
    chosen.sort(key=lambda item: (item[2], item[3], item[0]))
    return chosen


def write_path_device(path: Path, num_qubits: int) -> None:
    couplings = []
    for src in range(num_qubits - 1):
        couplings.append(f"{src}_{src + 1}")
        couplings.append(f"{src + 1}_{src}")
    payload = {
        "name": f"line_n{num_qubits}_bidirectional",
        "num_qubits": num_qubits,
        "basis_gates": ["rz", "sx", "x", "cx"],
        "cx_coupling": couplings,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate routing-pressure VQE proxy inputs.")
    parser.add_argument("--qubits", nargs="+", type=int, default=[8, 12, 16], help="Qubit counts to generate.")
    parser.add_argument("--reps", type=int, default=4, help="TwoLocal repetition count.")
    parser.add_argument("--levels", type=int, default=8, help="Number of routing-pressure levels per circuit size.")
    parser.add_argument("--random-samples", type=int, default=256, help="Number of random layout candidates per size.")
    parser.add_argument("--seed", type=int, default=7, help="Seed for reproducible layout generation.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/test_benchmark/routing_pressure_vqe"),
        help="Directory for generated QASM files.",
    )
    parser.add_argument(
        "--device-dir",
        type=Path,
        default=Path("data/devices/routing_pressure"),
        help="Directory for generated path-device JSON files.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/routing_pressure_manifest.csv"),
        help="Manifest CSV describing each fixed-layout benchmark case.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.device_dir.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for num_qubits in args.qubits:
        qc = build_vqe_proxy(num_qubits, args.reps)
        qasm_path = args.output_dir / f"vqe_layoutstress_n{num_qubits}.qasm"
        qasm_path.write_text(qasm2_dumps(qc), encoding="utf-8")

        device_path = args.device_dir / f"line_n{num_qubits}_bidirectional.json"
        write_path_device(device_path, num_qubits)

        weights = count_weighted_interactions(qc)
        selected = select_pressure_levels(
            weights,
            num_qubits=num_qubits,
            levels=args.levels,
            seed=args.seed + num_qubits,
            samples=args.random_samples,
        )
        for level_index, (tag, layout, pressure, max_hop) in enumerate(selected):
            rows.append(
                {
                    "case_id": f"n{num_qubits}_level{level_index:02d}",
                    "num_qubits": num_qubits,
                    "reps": args.reps,
                    "circuit_file": str(qasm_path.resolve()),
                    "device_file": str(device_path.resolve()),
                    "layout_tag": tag,
                    "level_index": level_index,
                    "routing_pressure": f"{pressure:.6f}",
                    "max_hop": max_hop,
                    "initial_layout": json.dumps(layout),
                }
            )

    with args.manifest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "case_id",
                "num_qubits",
                "reps",
                "circuit_file",
                "device_file",
                "layout_tag",
                "level_index",
                "routing_pressure",
                "max_hop",
                "initial_layout",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote manifest to {args.manifest}")


if __name__ == "__main__":
    main()
