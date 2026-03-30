#!/usr/bin/env python3
"""Generate a >100-qubit Brisbane-derived routing-pressure benchmark."""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import Counter, deque
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from qiskit import QuantumCircuit
from qiskit.qasm2 import dumps as qasm2_dumps


DirectedEdge = Tuple[int, int]
UndirectedEdge = Tuple[int, int]


def load_directed_device(path: Path) -> Tuple[int, List[DirectedEdge], Dict[int, set[int]]]:
    cfg = json.loads(path.read_text(encoding="utf-8"))
    num_qubits = int(cfg["num_qubits"])
    directed_edges: List[DirectedEdge] = []
    undirected_adj: Dict[int, set[int]] = {idx: set() for idx in range(num_qubits)}
    for entry in cfg["cx_coupling"]:
        src, dst = map(int, str(entry).split("_"))
        directed_edges.append((src, dst))
        undirected_adj[src].add(dst)
        undirected_adj[dst].add(src)
    return num_qubits, directed_edges, undirected_adj


def choose_connected_subset(adjacency: Dict[int, set[int]], subset_size: int) -> List[int]:
    if subset_size > len(adjacency):
        raise ValueError("subset_size exceeds device size")
    start = max(adjacency, key=lambda node: (len(adjacency[node]), -node))
    order: List[int] = []
    seen = {start}
    queue = deque([start])
    while queue and len(order) < subset_size:
        node = queue.popleft()
        order.append(node)
        for nxt in sorted(adjacency[node], key=lambda nb: (-len(adjacency[nb]), nb)):
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
        if len(queue) == 0 and len(order) < subset_size:
            remaining = [node for node in adjacency if node not in seen]
            if remaining:
                extra = max(remaining, key=lambda node: (len(adjacency[node]), -node))
                seen.add(extra)
                queue.append(extra)
    if len(order) != subset_size:
        raise RuntimeError(f"Unable to find connected subset of size {subset_size}")
    return order


def make_subdevice(
    directed_edges: Sequence[DirectedEdge],
    selected_nodes: Sequence[int],
    output_path: Path,
) -> Tuple[Dict[int, int], Dict[int, set[int]], List[DirectedEdge], List[UndirectedEdge]]:
    old_to_new = {old: new for new, old in enumerate(selected_nodes)}
    local_directed: List[DirectedEdge] = []
    local_undirected_set = set()
    local_adj: Dict[int, set[int]] = {idx: set() for idx in range(len(selected_nodes))}
    for src, dst in directed_edges:
        if src not in old_to_new or dst not in old_to_new:
            continue
        new_src = old_to_new[src]
        new_dst = old_to_new[dst]
        local_directed.append((new_src, new_dst))
        a, b = sorted((new_src, new_dst))
        local_undirected_set.add((a, b))
        local_adj[new_src].add(new_dst)
        local_adj[new_dst].add(new_src)

    payload = {
        "name": f"ibm_brisbane_sub{len(selected_nodes)}",
        "num_qubits": len(selected_nodes),
        "basis_gates": ["rz", "sx", "x", "cx"],
        "cx_coupling": [f"{src}_{dst}" for src, dst in local_directed],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    local_undirected = sorted(local_undirected_set)
    return old_to_new, local_adj, local_directed, local_undirected


def greedy_edge_coloring(edges: Sequence[UndirectedEdge]) -> List[List[UndirectedEdge]]:
    colors: List[List[UndirectedEdge]] = []
    used_by_vertex: Dict[int, set[int]] = {}
    for u, v in sorted(edges):
        used = used_by_vertex.get(u, set()) | used_by_vertex.get(v, set())
        color = 0
        while color in used:
            color += 1
        while color >= len(colors):
            colors.append([])
        colors[color].append((u, v))
        used_by_vertex.setdefault(u, set()).add(color)
        used_by_vertex.setdefault(v, set()).add(color)
    return colors


def select_oriented_layers(
    undirected_layers: Sequence[Sequence[UndirectedEdge]],
    directed_edges: Sequence[DirectedEdge],
) -> List[List[DirectedEdge]]:
    directed_lookup = set(directed_edges)
    layers: List[List[DirectedEdge]] = []
    for layer in undirected_layers:
        directed_layer: List[DirectedEdge] = []
        for u, v in layer:
            if (u, v) in directed_lookup:
                directed_layer.append((u, v))
            elif (v, u) in directed_lookup:
                directed_layer.append((v, u))
            else:
                raise RuntimeError(f"Missing directed orientation for edge {(u, v)}")
        layers.append(directed_layer)
    return layers


def build_circuit(num_qubits: int, layers: Sequence[Sequence[DirectedEdge]], reps: int, variant_seed: int) -> QuantumCircuit:
    rng = random.Random(variant_seed)
    qc = QuantumCircuit(num_qubits, name=f"brisbane_topology_n{num_qubits}_v{variant_seed}")
    for rep in range(reps):
        for qubit in range(num_qubits):
            theta = 0.031 * ((variant_seed + 1) * (rep + 1) * (qubit + 1))
            qc.ry(theta, qubit)
        layer_order = list(range(len(layers)))
        rng.shuffle(layer_order)
        for layer_idx in layer_order:
            layer = list(layers[layer_idx])
            rng.shuffle(layer)
            for ctrl, tgt in layer:
                qc.cx(ctrl, tgt)
        for qubit in range(num_qubits):
            theta = 0.017 * ((variant_seed + 3) * (rep + 2) * (qubit + 1))
            qc.rz(theta, qubit)
    return qc


def interaction_weights(layers: Sequence[Sequence[DirectedEdge]], reps: int) -> Counter[UndirectedEdge]:
    weights: Counter[UndirectedEdge] = Counter()
    for _ in range(reps):
        for layer in layers:
            for ctrl, tgt in layer:
                edge = (ctrl, tgt) if ctrl < tgt else (tgt, ctrl)
                weights[edge] += 1
    return weights


def shortest_paths(adjacency: Dict[int, set[int]]) -> List[List[int]]:
    num_nodes = len(adjacency)
    distances = [[10**9] * num_nodes for _ in range(num_nodes)]
    for src in range(num_nodes):
        distances[src][src] = 0
        queue = deque([src])
        while queue:
            node = queue.popleft()
            for nxt in adjacency[node]:
                if distances[src][nxt] > distances[src][node] + 1:
                    distances[src][nxt] = distances[src][node] + 1
                    queue.append(nxt)
    return distances


def routing_pressure(weights: Counter[UndirectedEdge], layout: Sequence[int], distances: Sequence[Sequence[int]]) -> float:
    total_weight = sum(weights.values())
    weighted = 0.0
    for (u, v), weight in weights.items():
        dist = distances[layout[u]][layout[v]]
        weighted += weight * max(0, dist - 1)
    return weighted / max(1, total_weight)


def build_weighted_neighbors(weights: Counter[UndirectedEdge], num_qubits: int) -> List[Dict[int, int]]:
    neighbors: List[Dict[int, int]] = [dict() for _ in range(num_qubits)]
    for (u, v), weight in weights.items():
        neighbors[u][v] = weight
        neighbors[v][u] = weight
    return neighbors


def excess_distance_matrix(distances: Sequence[Sequence[int]]) -> List[List[int]]:
    return [[max(0, dist - 1) for dist in row] for row in distances]


def pressure_numerator(weights: Counter[UndirectedEdge], layout: Sequence[int], excess: Sequence[Sequence[int]]) -> int:
    total = 0
    for (u, v), weight in weights.items():
        total += weight * excess[layout[u]][layout[v]]
    return total


def layout_hamming(layout: Sequence[int]) -> int:
    return sum(1 for logical, physical in enumerate(layout) if logical != physical)


def swap_delta_numerator(
    layout: Sequence[int],
    logical_a: int,
    logical_b: int,
    neighbors: Sequence[Dict[int, int]],
    excess: Sequence[Sequence[int]],
) -> int:
    phys_a = layout[logical_a]
    phys_b = layout[logical_b]
    if phys_a == phys_b:
        return 0
    delta = 0
    for other, weight in neighbors[logical_a].items():
        if other == logical_b:
            continue
        phys_other = layout[other]
        delta += weight * (excess[phys_b][phys_other] - excess[phys_a][phys_other])
    for other, weight in neighbors[logical_b].items():
        if other == logical_a:
            continue
        phys_other = layout[other]
        delta += weight * (excess[phys_a][phys_other] - excess[phys_b][phys_other])
    return delta


def deterministic_candidate_pool(
    num_qubits: int,
    weights: Counter[UndirectedEdge],
    distances: Sequence[Sequence[int]],
    random_samples: int,
    seed: int,
    max_pressure: float | None = None,
) -> List[Dict[str, object]]:
    rng = random.Random(seed)
    identity = list(range(num_qubits))
    candidates = []
    seen = {tuple(identity)}
    swap_budgets = sorted({1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64})
    per_budget = max(4, random_samples // max(1, len(swap_budgets)))
    sample_idx = 0
    for budget in swap_budgets:
        for _ in range(per_budget):
            layout = identity.copy()
            for __ in range(budget):
                i = rng.randrange(num_qubits)
                j = rng.randrange(num_qubits)
                layout[i], layout[j] = layout[j], layout[i]
            key = tuple(layout)
            if key in seen:
                continue
            seen.add(key)
            pressure = routing_pressure(weights, layout, distances)
            candidates.append(
                {
                    "layout": layout,
                    "pressure": pressure,
                    "tag": f"swap{budget:02d}_{sample_idx:04d}",
                }
            )
            sample_idx += 1
    while sample_idx < random_samples:
        layout = identity.copy()
        rng.shuffle(layout)
        key = tuple(layout)
        if key in seen:
            continue
        seen.add(key)
        pressure = routing_pressure(weights, layout, distances)
        candidates.append({"layout": layout, "pressure": pressure, "tag": f"random_{sample_idx:04d}"})
        sample_idx += 1
    if max_pressure is not None:
        candidates = [item for item in candidates if item["pressure"] <= max_pressure]
    candidates.sort(key=lambda item: item["pressure"])
    return [{"layout": identity, "pressure": 0.0, "tag": "identity"}] + candidates


def refine_layout_toward_target(
    target_pressure: float,
    seed_layout: Sequence[int],
    weights: Counter[UndirectedEdge],
    neighbors: Sequence[Dict[int, int]],
    excess: Sequence[Sequence[int]],
    total_weight: int,
    max_steps: int = 96,
) -> Dict[str, object]:
    layout = list(seed_layout)
    numerator = pressure_numerator(weights, layout, excess)
    hamming = layout_hamming(layout)
    current_pressure = numerator / max(1, total_weight)
    current_key = (abs(current_pressure - target_pressure), 0 if current_pressure <= target_pressure else 1, hamming)

    for _ in range(max_steps):
        best = None
        for logical_a in range(len(layout) - 1):
            for logical_b in range(logical_a + 1, len(layout)):
                delta = swap_delta_numerator(layout, logical_a, logical_b, neighbors, excess)
                new_numerator = numerator + delta
                new_pressure = new_numerator / max(1, total_weight)
                old_hd = int(layout[logical_a] != logical_a) + int(layout[logical_b] != logical_b)
                new_hd = int(layout[logical_b] != logical_a) + int(layout[logical_a] != logical_b)
                candidate_hamming = hamming - old_hd + new_hd
                candidate_key = (
                    abs(new_pressure - target_pressure),
                    0 if new_pressure <= target_pressure else 1,
                    candidate_hamming,
                    logical_a,
                    logical_b,
                )
                if candidate_key < current_key and (best is None or candidate_key < best[0]):
                    best = (candidate_key, logical_a, logical_b, new_numerator, new_pressure, candidate_hamming)
        if best is None:
            break
        _, logical_a, logical_b, numerator, current_pressure, hamming = best
        layout[logical_a], layout[logical_b] = layout[logical_b], layout[logical_a]
        current_key = (abs(current_pressure - target_pressure), 0 if current_pressure <= target_pressure else 1, hamming)

    return {
        "layout": layout,
        "pressure": current_pressure,
        "pressure_error": current_pressure - target_pressure,
        "hamming": hamming,
    }


def select_target_layouts(
    num_qubits: int,
    weights: Counter[UndirectedEdge],
    distances: Sequence[Sequence[int]],
    targets: Sequence[float],
    layouts_per_target: int,
    random_samples: int,
    seed: int,
    max_pressure: float | None = None,
    search_seeds: int = 10,
) -> List[Dict[str, object]]:
    total_weight = sum(weights.values())
    neighbors = build_weighted_neighbors(weights, num_qubits)
    excess = excess_distance_matrix(distances)
    pool = deterministic_candidate_pool(
        num_qubits,
        weights=weights,
        distances=distances,
        random_samples=random_samples,
        seed=seed,
        max_pressure=max_pressure,
    )

    selected: List[Dict[str, object]] = []
    prior_layouts: List[Dict[str, object]] = [pool[0]]
    for target_index, target_pressure in enumerate(targets):
        seed_candidates: List[Dict[str, object]] = []
        seen = set()
        for item in prior_layouts:
            key = tuple(item["layout"])
            if key not in seen:
                seen.add(key)
                seed_candidates.append(item)
        for item in sorted(
            pool,
            key=lambda cand: (
                abs(float(cand["pressure"]) - target_pressure),
                layout_hamming(cand["layout"]),
                str(cand["tag"]),
            ),
        ):
            key = tuple(item["layout"])
            if key in seen:
                continue
            seen.add(key)
            seed_candidates.append(item)
            if len(seed_candidates) >= max(layouts_per_target * 4, search_seeds):
                break

        refined: List[Dict[str, object]] = []
        refined_seen = set()
        for seed_index, item in enumerate(seed_candidates):
            refined_item = refine_layout_toward_target(
                target_pressure,
                seed_layout=item["layout"],
                weights=weights,
                neighbors=neighbors,
                excess=excess,
                total_weight=total_weight,
            )
            key = tuple(refined_item["layout"])
            if key in refined_seen:
                continue
            refined_seen.add(key)
            refined_item["tag"] = f"target{target_pressure:.1f}_seed{seed_index:02d}"
            refined.append(refined_item)

        refined.sort(
            key=lambda item: (
                abs(float(item["pressure"]) - target_pressure),
                abs(float(item["pressure_error"])),
                int(item["hamming"]),
                str(item["tag"]),
            )
        )

        picked: List[Dict[str, object]] = []
        for layout_index, item in enumerate(refined[:layouts_per_target]):
            picked_item = dict(item)
            picked_item["pressure_bucket"] = target_index
            picked_item["target_pressure"] = target_pressure
            picked_item["layout_index"] = layout_index
            picked.append(picked_item)
        while len(picked) < layouts_per_target and picked:
            dup = dict(picked[0])
            dup["layout_index"] = len(picked)
            dup["tag"] = f"{dup['tag']}_dup{len(picked):02d}"
            picked.append(dup)
        if not picked:
            raise RuntimeError(f"Unable to find layout for target pressure {target_pressure}")
        selected.extend(picked)
        prior_layouts = picked
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a large-device routing-pressure benchmark manifest.")
    parser.add_argument("--device", type=Path, default=Path("data/devices/ibm_brisbane.json"))
    parser.add_argument("--num-qubits", type=int, default=110)
    parser.add_argument("--reps", type=int, default=20)
    parser.add_argument("--circuit-variants", type=int, default=5)
    parser.add_argument("--layouts-per-target", type=int, default=2)
    parser.add_argument("--pressure-step", type=float, default=0.2)
    parser.add_argument("--max-pressure", type=float, default=2.0)
    parser.add_argument("--random-layout-samples", type=int, default=2000)
    parser.add_argument("--search-seeds", type=int, default=12)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument(
        "--device-out",
        type=Path,
        default=Path("data/devices/routing_pressure/ibm_brisbane_sub110.json"),
    )
    parser.add_argument(
        "--circuit-dir",
        type=Path,
        default=Path("data/test_benchmark/routing_pressure_brisbane"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/brisbane_routing_pressure_manifest.csv"),
    )
    args = parser.parse_args()

    _, directed_edges, adjacency = load_directed_device(args.device)
    selected_nodes = choose_connected_subset(adjacency, args.num_qubits)
    _, local_adj, local_directed, local_undirected = make_subdevice(directed_edges, selected_nodes, args.device_out)
    layers = select_oriented_layers(greedy_edge_coloring(local_undirected), local_directed)
    weights = interaction_weights(layers, args.reps)
    distances = shortest_paths(local_adj)
    num_targets = int(round(args.max_pressure / args.pressure_step))
    targets = [round(idx * args.pressure_step, 10) for idx in range(num_targets + 1)]
    layouts = select_target_layouts(
        args.num_qubits,
        weights=weights,
        distances=distances,
        targets=targets,
        layouts_per_target=args.layouts_per_target,
        random_samples=args.random_layout_samples,
        seed=args.seed,
        max_pressure=args.max_pressure,
        search_seeds=args.search_seeds,
    )

    args.circuit_dir.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for circuit_variant in range(args.circuit_variants):
        qc = build_circuit(args.num_qubits, layers, args.reps, variant_seed=args.seed + circuit_variant)
        circuit_path = args.circuit_dir / f"brisbane_topology_n{args.num_qubits}_v{circuit_variant:02d}.qasm"
        circuit_path.write_text(qasm2_dumps(qc), encoding="utf-8")
        input_depth = qc.depth()
        for sample_idx, layout_info in enumerate(layouts):
            rows.append(
                {
                    "case_id": f"v{circuit_variant:02d}_t{layout_info['pressure_bucket']:02d}_l{layout_info['layout_index']:02d}",
                    "circuit_variant": circuit_variant,
                    "num_qubits": args.num_qubits,
                    "input_depth": input_depth,
                    "reps": args.reps,
                    "circuit_file": str(circuit_path.resolve()),
                    "device_file": str(args.device_out.resolve()),
                    "pressure_bucket": int(layout_info["pressure_bucket"]),
                    "target_pressure": f"{float(layout_info['target_pressure']):.6f}",
                    "routing_pressure": f"{float(layout_info['pressure']):.6f}",
                    "layout_tag": str(layout_info["tag"]),
                    "initial_layout": json.dumps(layout_info["layout"]),
                }
            )

    with args.manifest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "case_id",
                "circuit_variant",
                "num_qubits",
                "input_depth",
                "reps",
                "circuit_file",
                "device_file",
                "pressure_bucket",
                "target_pressure",
                "routing_pressure",
                "layout_tag",
                "initial_layout",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote manifest to {args.manifest}")
    print(f"Wrote subdevice to {args.device_out}")
    print(f"Generated {len(rows)} benchmark cases")


if __name__ == "__main__":
    main()
