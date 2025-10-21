#!/usr/bin/env python3
import argparse
import copy
import json
import math
import re
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from qutip_qtrl import pulseoptim
import qutip_qtrl.pulsegen as pulsegen

try:
    from device_pulse_fab.simulate_ankaa9q import simulate_pulse_schedule, phase_unitary
except ModuleNotFoundError:
    sys.path.append(str(Path(__file__).resolve().parent))
    from simulate_ankaa9q import simulate_pulse_schedule, phase_unitary


def evaluate_expression(expr) -> float:
    if expr is None:
        return 0.0
    if isinstance(expr, (int, float, np.generic)):
        return float(expr)
    if isinstance(expr, str):
        expr = expr.strip()
        if not expr:
            return 0.0
        safe_env = {"pi": math.pi}
        return float(eval(expr, {"__builtins__": None}, safe_env))
    raise TypeError(f"unsupported numeric expression type: {type(expr)}")


def normalize_gate_label(gate: str) -> Tuple[str, float | None]:
    gate_l = gate.lower()
    if gate_l == "x":
        return "rx", math.pi
    if gate_l == "sx":
        return "rx", 0.5 * math.pi
    return gate_l, None


def canonicalize_pulse_library(pulse_doc: dict) -> List[dict]:
    lib = pulse_doc.get("pulse_library")
    if not isinstance(lib, list):
        lib = pulse_doc.get("pulse_definitions")
    if not isinstance(lib, list):
        return []
    for entry in lib:
        gate_norm, default_theta = normalize_gate_label(str(entry.get("gate", "")))
        entry["gate"] = gate_norm
        params = entry.setdefault("parameters", {})
        if default_theta is not None:
            params.setdefault("theta", default_theta)
        if "width" in entry:
            entry["width"] = evaluate_expression(entry["width"])
        if "amplitude" in entry:
            entry["amplitude"] = evaluate_expression(entry["amplitude"])
        if "duration" in entry:
            entry["duration"] = evaluate_expression(entry["duration"])
        if "frequency" in entry:
            entry["frequency"] = evaluate_expression(entry["frequency"])
        for key in ("samples", "samples_i", "samples_q"):
            arr = entry.get(key)
            if isinstance(arr, Sequence) and not isinstance(arr, (bytes, str)):
                entry[key] = [evaluate_expression(v) for v in arr]
    return lib


def extract_width(entry: dict) -> float:
    width = entry.get("width")
    if width is None or width == 0:
        params = entry.get("parameters") or {}
        dur_ns = params.get("duration_ns")
        if dur_ns is not None:
            width = float(dur_ns) * 1e-9
    width = evaluate_expression(width)
    return float(width)


def extract_theta(entry: dict) -> float:
    params = entry.get("parameters") or {}
    theta = evaluate_expression(params.get("theta"))
    if abs(theta) > 1e-12:
        return float(theta)
    amp = evaluate_expression(entry.get("amplitude"))
    width = extract_width(entry)
    if width > 0.0 and amp != 0.0:
        return float(amp * width)
    return float(theta)


def extract_samples(entry: dict, key: str) -> np.ndarray:
    samples = entry.get(key)
    if samples is None and key == "samples_i":
        samples = entry.get("samples")
    if samples is None:
        return np.array([], dtype=float)
    if isinstance(samples, np.ndarray):
        return samples.astype(float)
    try:
        return np.asarray(samples, dtype=float)
    except Exception:
        return np.array([], dtype=float)


def resample_to_length(samples: np.ndarray, width: float, length: int) -> np.ndarray:
    length = max(1, int(length))
    if samples.size == 0:
        return np.zeros(length, dtype=float)
    if samples.size == 1:
        return np.full(length, float(samples[0]), dtype=float)
    if samples.size == length:
        return samples.astype(float)
    source = np.linspace(0.0, width, samples.size, dtype=float)
    target = np.linspace(0.0, width, length, dtype=float)
    return np.interp(target, source, samples).astype(float)


@dataclass
class PulseSpec:
    gate: str
    theta: float
    qubits: Tuple[int, ...]
    entry: dict
    width: float
    local_qubits: Optional[Tuple[int, ...]] = None


@dataclass
class ControlBasis:
    ctrls: List[qt.Qobj]
    labels: List[str]
    ix_index: Dict[int, int]
    iy_index: Dict[int, int]
    coupling_index: Dict[Tuple[int, int], int]


@dataclass
class MergeCandidate:
    gate: str
    qubits: Tuple[int, ...]
    score: float
    instance_id: Optional[int] = None


@dataclass
class OptimizationSettings:
    compression_ratio: Optional[float]
    compression_shift_ns: float
    evo_time_ns: Optional[float]
    amp_bound: Optional[float]
    seed_from_library: bool
    min_merged_fidelity: float
    max_iter: int = 400
    max_wall_time: int = 120


@dataclass
class OptimizationResult:
    candidate: MergeCandidate
    seq: List[PulseSpec]
    unique_qubits: List[int]
    logical_to_physical: Dict[int, int]
    basis: ControlBasis
    control_waveforms: Dict[str, np.ndarray]
    evo_time: float
    dt_effective: float
    native_time: float
    native_fidelity: float
    optimized_fidelity: float
    requested_ratio: float
    actual_ratio: float
    final_fid_err: float
    n_ts: int
    plot_path: Optional[Path]


CONTROL_SINGLE_RE = re.compile(r"([IQ])(\d+)$")
CONTROL_COUPLING_RE = re.compile(r"J(\d+)(\d+)$")


def parse_control_label(label: str) -> Tuple[str, Tuple[int, ...]]:
    match_single = CONTROL_SINGLE_RE.fullmatch(label)
    if match_single:
        axis = match_single.group(1)
        qubit = int(match_single.group(2))
        return axis, (qubit,)
    match_coupling = CONTROL_COUPLING_RE.fullmatch(label)
    if match_coupling:
        q0 = int(match_coupling.group(1))
        q1 = int(match_coupling.group(2))
        return "J", (q0, q1)
    raise ValueError(f"unrecognized control label: {label}")


def collect_single_qubit_specs(pulse_lib: Sequence[dict], gate_names: Tuple[str, ...] = ("rx",)) -> Dict[Tuple[int, ...], PulseSpec]:
    specs: Dict[Tuple[int, ...], PulseSpec] = {}
    for entry in pulse_lib:
        gate = str(entry.get("gate", "")).lower()
        if gate not in gate_names:
            continue
        qubits_raw = entry.get("qubits") or []
        if len(qubits_raw) != 1:
            continue
        if entry.get("virtual"):
            continue
        qubit = int(qubits_raw[0])
        key = (qubit,)
        if key in specs:
            continue
        theta = extract_theta(entry)
        width = extract_width(entry)
        specs[key] = PulseSpec(gate, theta, (qubit,), entry, width)
    return specs


def collect_two_qubit_specs(pulse_lib: Sequence[dict], gate_name: str = "iswap") -> Dict[Tuple[int, int], PulseSpec]:
    specs: Dict[Tuple[int, int], PulseSpec] = {}
    gate_name = gate_name.lower()
    for entry in pulse_lib:
        gate = str(entry.get("gate", "")).lower()
        if gate != gate_name:
            continue
        qubits_raw = entry.get("qubits") or []
        if len(qubits_raw) != 2:
            continue
        if entry.get("virtual"):
            continue
        physical = tuple(int(q) for q in qubits_raw)
        key = tuple(sorted(physical))
        if key in specs:
            continue
        theta = extract_theta(entry)
        width = extract_width(entry)
        specs[key] = PulseSpec(gate, theta, physical, entry, width)
    return specs


def extract_native_gates(device_doc: dict) -> set[str]:
    gates: set[str] = set()
    basis = device_doc.get("basis_gates") or []
    for gate in basis:
        if isinstance(gate, str):
            gates.add(gate.lower())
    metadata = device_doc.get("metadata")
    if isinstance(metadata, dict):
        meta_basis = metadata.get("basis_gates") or []
        for gate in meta_basis:
            if isinstance(gate, str):
                gates.add(gate.lower())
    return gates


def aggregate_schedule_candidates(schedule: Sequence[dict]) -> Dict[Tuple[str, Tuple[int, ...]], float]:
    totals: Dict[Tuple[str, Tuple[int, ...]], float] = {}
    for event in schedule:
        gate_raw = str(event.get("gate", "")).lower()
        gate_norm, _ = normalize_gate_label(gate_raw)
        qubits_raw = event.get("qubits") or []
        if not qubits_raw:
            continue
        qubits = tuple(int(q) for q in qubits_raw)
        if len(qubits) > 2:
            continue
        duration = float(event.get("duration", 0.0) or 0.0)
        if duration <= 0.0:
            continue
        if len(qubits) == 2:
            qubit_key = tuple(sorted(qubits))
        else:
            qubit_key = qubits
        if gate_norm not in ("rx", "iswap"):
            continue
        totals[(gate_norm, qubit_key)] = totals.get((gate_norm, qubit_key), 0.0) + duration
    return totals


def extract_top_gate_entries(pulse_doc: dict) -> List[dict]:
    top_list = pulse_doc.get("top_gates")
    if not isinstance(top_list, list):
        top_list = (pulse_doc.get("critical_path") or {}).get("top_gates")
    if not isinstance(top_list, list):
        return []
    entries: List[dict] = []
    for item in top_list:
        if isinstance(item, dict) and item.get("gate"):
            entries.append(item)
    return entries


def select_merge_candidates(
    rx_specs: Dict[Tuple[int, ...], PulseSpec],
    iswap_specs: Dict[Tuple[int, int], PulseSpec],
    merge_limit: Optional[int],
    native_gates: Optional[set[str]] = None,
    schedule_totals: Optional[Dict[Tuple[str, Tuple[int, ...]], float]] = None,
    instance_scores: Optional[List[Tuple[str, Tuple[int, ...], int, float]]] = None,
    top_gate_entries: Optional[Sequence[dict]] = None,
    logical_instances: Optional[Dict[int, dict]] = None,
) -> List[MergeCandidate]:
    candidates: List[MergeCandidate] = []
    used_instances: set[int] = set()

    if top_gate_entries and logical_instances:
        gate_to_instances: Dict[str, List[Tuple[int, dict]]] = {}
        for inst_id, info in logical_instances.items():
            label = info.get("label")
            if not label:
                continue
            if native_gates and label in native_gates:
                continue
            gate_to_instances.setdefault(label, []).append((inst_id, info))
        for label, entries in gate_to_instances.items():
            entries.sort(
                key=lambda item: (
                    -float(item[1].get("duration") or 0.0),
                    float(item[1].get("first_index") or 0.0),
                )
            )
        sorted_top = sorted(
            (
                entry
                for entry in top_gate_entries
                if isinstance(entry, dict) and entry.get("gate")
            ),
            key=lambda entry: float(entry.get("total_duration") or 0.0),
            reverse=True,
        )
        for entry in sorted_top:
            label = str(entry.get("gate", "")).lower()
            if not label or (native_gates and label in native_gates):
                continue
            available = gate_to_instances.get(label)
            if not available:
                continue
            count = int(entry.get("count") or 1)
            added = 0
            for inst_id, info in available:
                if inst_id in used_instances:
                    continue
                qubit_order = info.get("ordered_qubits")
                if not qubit_order:
                    qubit_order = sorted(info.get("qubits") or [])
                qubits_tuple = tuple(int(q) for q in qubit_order)
                if not qubits_tuple:
                    continue
                score = float(entry.get("total_duration") or info.get("duration") or 0.0)
                candidates.append(
                    MergeCandidate(label, qubits_tuple, score, instance_id=inst_id)
                )
                used_instances.add(inst_id)
                added += 1
                if merge_limit and merge_limit > 0 and len(candidates) >= merge_limit:
                    break
                if added >= count:
                    break
            if merge_limit and merge_limit > 0 and len(candidates) >= merge_limit:
                break

    if merge_limit and merge_limit > 0 and len(candidates) >= merge_limit:
        return candidates[: merge_limit]

    fallback: List[MergeCandidate] = []
    if instance_scores:
        for gate, qubits, instance_id, score in instance_scores:
            if native_gates and gate in native_gates:
                continue
            if instance_id in used_instances:
                continue
            fallback.append(MergeCandidate(gate, qubits, score, instance_id=instance_id))
    elif schedule_totals:
        for (gate, qubits), score in schedule_totals.items():
            if native_gates and gate in native_gates:
                continue
            if gate == "rx" and qubits in rx_specs:
                fallback.append(MergeCandidate("rx", qubits, score))
            elif gate == "iswap":
                key = tuple(sorted(qubits))
                if key in iswap_specs:
                    fallback.append(MergeCandidate("iswap", key, score))
    else:
        for qubits, spec in rx_specs.items():
            fallback.append(MergeCandidate("rx", qubits, spec.width))
        for qubits, spec in iswap_specs.items():
            fallback.append(MergeCandidate("iswap", qubits, spec.width))

    if not instance_scores:
        fallback.sort(key=lambda c: c.score, reverse=True)
    for cand in fallback:
        candidates.append(cand)
        if merge_limit and merge_limit > 0 and len(candidates) >= merge_limit:
            break
    return candidates


def build_sequence_for_candidate(
    candidate: MergeCandidate,
    rx_specs: Dict[Tuple[int, ...], PulseSpec],
    iswap_specs: Dict[Tuple[int, int], PulseSpec],
    logical_instances: Optional[Dict[int, dict]] = None,
    library_by_id: Optional[Dict[str, dict]] = None,
) -> Optional[Tuple[List[PulseSpec], List[int], Dict[int, int]]]:
    if candidate.instance_id is not None and logical_instances and library_by_id:
        instance = logical_instances.get(candidate.instance_id)
        result = build_sequence_from_instance(candidate, instance, library_by_id)
        if result is not None:
            return result
    if candidate.gate == "cx":
        if len(candidate.qubits) != 2:
            return None
        ctrl, target = candidate.qubits
        pair = tuple(sorted(candidate.qubits))
        spec_iswap = iswap_specs.get(pair)
        spec_rx_ctrl = rx_specs.get((ctrl,))
        spec_rx_target = rx_specs.get((target,))
        if spec_iswap is None or spec_rx_ctrl is None or spec_rx_target is None:
            return None
        seq: List[PulseSpec] = []
        unique_qubits = [ctrl, target]
        logical_to_physical = {0: ctrl, 1: target}

        local_map = {ctrl: 0, target: 1}

        def localized_qubits(spec: PulseSpec) -> Tuple[int, ...]:
            return tuple(local_map.get(q, 0) for q in spec.qubits)

        seq.append(replace(spec_iswap, local_qubits=localized_qubits(spec_iswap)))
        seq.append(replace(spec_rx_ctrl, theta=-0.5 * math.pi, local_qubits=localized_qubits(spec_rx_ctrl)))
        seq.append(replace(spec_iswap, local_qubits=localized_qubits(spec_iswap)))
        seq.append(replace(spec_rx_target, theta=-0.5 * math.pi, local_qubits=localized_qubits(spec_rx_target)))

        return seq, unique_qubits, logical_to_physical

    if candidate.gate == "rx":
        spec = rx_specs.get(candidate.qubits)
        if spec is None:
            return None
        physical_qubit = spec.qubits[0]
        seq = [replace(spec, local_qubits=(0,))]
        unique_qubits = [physical_qubit]
        logical_to_physical = {0: physical_qubit}
        return seq, unique_qubits, logical_to_physical

    if candidate.gate == "iswap":
        key = tuple(sorted(candidate.qubits))
        spec = iswap_specs.get(key)
        if spec is None:
            return None
        q0, q1 = spec.qubits
        seq: List[PulseSpec] = []
        spec_rx0 = rx_specs.get((q0,))
        spec_rx1 = rx_specs.get((q1,))
        if spec_rx0 is not None:
            seq.append(replace(spec_rx0, local_qubits=(0,)))
        seq.append(replace(spec, local_qubits=(0, 1)))
        if spec_rx1 is not None:
            seq.append(replace(spec_rx1, local_qubits=(1,)))
        unique_qubits = [q0, q1]
        logical_to_physical = {0: q0, 1: q1}
        return seq, unique_qubits, logical_to_physical

    return None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as h:
        return json.load(h)


def single_xy_unitary(theta: float, phase: float) -> qt.Qobj:
    return (-1j * theta / 2.0 * (math.cos(phase) * qt.sigmax() + math.sin(phase) * qt.sigmay())).expm()


def two_qubit_xy_unitary(theta: float) -> qt.Qobj:
    X = qt.sigmax()
    Y = qt.sigmay()
    H = -1j * theta / 2.0 * (qt.tensor(X, X) + qt.tensor(Y, Y))
    return H.expm()


def build_sequence_from_instance(
    candidate: MergeCandidate,
    instance: Optional[dict],
    library_by_id: Dict[str, dict],
) -> Optional[Tuple[List[PulseSpec], List[int], Dict[int, int]]]:
    if not instance:
        return None
    events: List[dict] = instance.get("events") or []
    if not events:
        return None
    unique_qubits: List[int] = list(candidate.qubits)
    phys_to_local: Dict[int, int] = {q: idx for idx, q in enumerate(unique_qubits)}
    seq: List[PulseSpec] = []
    for event in events:
        qubits_raw = event.get("qubits") or []
        if not qubits_raw:
            continue
        qubits = tuple(int(q) for q in qubits_raw)
        for q in qubits:
            if q not in phys_to_local:
                phys_to_local[q] = len(unique_qubits)
                unique_qubits.append(q)
        local_qubits = tuple(phys_to_local[q] for q in qubits)
        gate = str(event.get("gate", "")).lower()
        duration = float(event.get("duration", 0.0) or 0.0)
        params = event.get("parameters") or {}
        theta = evaluate_expression(params.get("theta")) if params else None
        if event.get("virtual"):
            continue
        pulse_id = event.get("pulse_id")
        if not pulse_id:
            continue
        library_entry = library_by_id.get(pulse_id)
        if library_entry is None:
            continue
        gate_from_library = str(library_entry.get("gate", "")).lower()
        if gate_from_library not in ("rx", "iswap"):
            continue
        entry_copy = copy.deepcopy(library_entry)
        entry_copy["qubits"] = list(qubits)
        if duration > 0.0:
            entry_copy["width"] = duration
        if params:
            entry_copy.setdefault("parameters", {}).update(params)
        if theta is None:
            theta = extract_theta(entry_copy)
        width = float(entry_copy.get("width", duration))
        seq.append(PulseSpec(gate, theta or 0.0, qubits, entry_copy, width, local_qubits=local_qubits))
    if not seq:
        return None
    logical_to_physical = {idx: q for idx, q in enumerate(unique_qubits)}
    return seq, unique_qubits, logical_to_physical


def embed_single(U: qt.Qobj, target: int, nq: int) -> qt.Qobj:
    ops = [qt.qeye(2)] * nq
    ops[target] = U
    return qt.tensor(ops)


def embed_two(U: qt.Qobj, q0: int, q1: int, nq: int) -> qt.Qobj:
    q0, q1 = sorted((q0, q1))
    dim = 2 ** nq
    data = np.zeros((dim, dim), dtype=complex)
    for c in range(dim):
        bc = [(c >> (nq - 1 - k)) & 1 for k in range(nq)]
        for r in range(dim):
            br = [(r >> (nq - 1 - k)) & 1 for k in range(nq)]
            ok = True
            for k in range(nq):
                if k == q0 or k == q1:
                    continue
                if br[k] != bc[k]:
                    ok = False
                    break
            if not ok:
                continue
            sc = (bc[q0] << 1) | bc[q1]
            sr = (br[q0] << 1) | br[q1]
            data[r, c] = U.full()[sr, sc]
    return qt.Qobj(data, dims=[[2] * nq, [2] * nq])


def local_qubits(spec: PulseSpec) -> Tuple[int, ...]:
    if spec.local_qubits is not None:
        return spec.local_qubits
    return spec.qubits


def ideal_from_sequence(seq: List[PulseSpec], nq: int) -> qt.Qobj:
    U = qt.qeye([2] * nq)
    for spec in seq:
        lq = local_qubits(spec)
        if spec.gate == "rx" and len(lq) == 1:
            U = embed_single(single_xy_unitary(spec.theta, 0.0), lq[0], nq) * U
        elif spec.gate == "iswap" and len(lq) == 2:
            U = embed_two(two_qubit_xy_unitary(spec.theta), lq[0], lq[1], nq) * U
        else:
            raise ValueError("unsupported gate")
    return U


def unitary_fidelity(U_targ: qt.Qobj, U_test: qt.Qobj) -> float:
    overlap = (U_targ.dag() * U_test).tr()
    return float(abs(overlap) / U_targ.shape[0])


def pauli_ops(nq: int):
    e0 = qt.basis(2, 0)
    e1 = qt.basis(2, 1)
    s01 = e0 * e1.dag()
    s10 = s01.dag()
    X = s01 + s10
    Y = -1j * s01 + 1j * s10
    def lift_one(op, i):
        ops = [qt.qeye(2)] * nq
        ops[i] = op
        return qt.tensor(ops)
    def lift_two(opi, i, opj, j):
        ops = [qt.qeye(2)] * nq
        ops[i] = opi
        ops[j] = opj
        return qt.tensor(ops)
    IX = [lift_one(0.5 * X, q) for q in range(nq)]
    IY = [lift_one(0.5 * Y, q) for q in range(nq)]
    Hcpl = {}
    for i in range(nq):
        for j in range(i + 1, nq):
            Hcpl[(i, j)] = 0.5 * (lift_two(X, i, X, j) + lift_two(Y, i, Y, j))
    return IX, IY, Hcpl


def pick_sequence_from_library(pulse_lib: List[dict], max_qubits: int = 2) -> List[PulseSpec]:
    rx_by_qubit: Dict[int, List[PulseSpec]] = {}
    rx_all: List[PulseSpec] = []
    iswap_list: List[PulseSpec] = []
    for entry in pulse_lib:
        if entry.get("virtual"):
            continue
        gate = str(entry.get("gate", "")).lower()
        qubits_raw = entry.get("qubits") or []
        qubits = tuple(int(q) for q in qubits_raw)
        if not qubits:
            continue
        theta = extract_theta(entry)
        width = extract_width(entry)
        if abs(theta) <= 1e-9 or width <= 0.0:
            continue
        if gate == "rx" and len(qubits) == 1:
            spec = PulseSpec("rx", theta, qubits, entry, width)
            rx_all.append(spec)
            rx_by_qubit.setdefault(qubits[0], []).append(spec)
        elif gate == "iswap" and len(qubits) == 2:
            spec = PulseSpec("iswap", theta, qubits, entry, width)
            iswap_list.append(spec)

    if not iswap_list:
        raise ValueError("no usable iSWAP pulses found in library")
    iswap_spec = iswap_list[0]

    def pick_rx(qubit: int) -> Optional[PulseSpec]:
        options = rx_by_qubit.get(qubit)
        if options:
            return options[0]
        if rx_all:
            return rx_all[0]
        return None

    rx_start = pick_rx(iswap_spec.qubits[0])
    rx_end = pick_rx(iswap_spec.qubits[1])

    seq: List[PulseSpec] = []
    if rx_start is not None:
        seq.append(rx_start)
    seq.append(iswap_spec)
    if rx_end is not None:
        seq.append(rx_end)

    if len(seq) < 2:
        raise ValueError("unable to assemble RX-iSWAP-RX sequence from pulse library")

    unique_qubits: List[int] = []
    for spec in seq:
        for q in spec.qubits:
            if q not in unique_qubits:
                unique_qubits.append(q)
    if len(unique_qubits) > max_qubits:
        raise ValueError(f"sequence spans {len(unique_qubits)} qubits, but only {max_qubits} supported")

    mapping = {q: idx for idx, q in enumerate(unique_qubits)}
    seq = [
        replace(spec, local_qubits=tuple(mapping[q] for q in spec.qubits))
        for spec in seq
    ]
    return seq


def build_controls(nq: int, logical_to_physical: Optional[Dict[int, int]] = None) -> ControlBasis:
    IX, IY, Hcpl = pauli_ops(nq)
    ctrls: List[qt.Qobj] = []
    labels: List[str] = []
    ix_index: Dict[int, int] = {}
    iy_index: Dict[int, int] = {}
    coupling_index: Dict[Tuple[int, int], int] = {}
    for q in range(nq):
        phys_q = logical_to_physical[q] if logical_to_physical and q in logical_to_physical else q
        ix_index[q] = len(ctrls)
        ctrls.append(IX[q])
        labels.append(f"I{phys_q}")
        iy_index[q] = len(ctrls)
        ctrls.append(IY[q])
        labels.append(f"Q{phys_q}")
    for key in sorted(Hcpl):
        phys_i = logical_to_physical[key[0]] if logical_to_physical and key[0] in logical_to_physical else key[0]
        phys_j = logical_to_physical[key[1]] if logical_to_physical and key[1] in logical_to_physical else key[1]
        coupling_index[key] = len(ctrls)
        ctrls.append(Hcpl[key])
        labels.append(f"J{phys_i}{phys_j}")
    return ControlBasis(ctrls=ctrls, labels=labels, ix_index=ix_index, iy_index=iy_index, coupling_index=coupling_index)


def make_rx_envelopes(spec: PulseSpec, width: float, seg_time: float, seg_slots: int) -> Tuple[np.ndarray, np.ndarray]:
    theta = spec.theta
    if seg_slots <= 0:
        return np.zeros(0, dtype=float), np.zeros(0, dtype=float)
    samples_i = extract_samples(spec.entry, "samples_i")
    samples_q = extract_samples(spec.entry, "samples_q")
    if samples_i.size == 0:
        samples_i = np.ones(2, dtype=float)
    if samples_q.size == 0:
        samples_q = np.zeros_like(samples_i)
    if samples_q.size != samples_i.size:
        samples_q = resample_to_length(samples_q, width, samples_i.size)
    if samples_i.size > 1:
        dt = width / (samples_i.size - 1)
        area_i = float(np.trapezoid(samples_i, dx=dt))
    else:
        area_i = float(samples_i[0] * width)
    if abs(area_i) < 1e-12:
        samples_i = np.ones_like(samples_i, dtype=float)
        samples_q = np.zeros_like(samples_i, dtype=float)
        if samples_i.size > 1:
            dt = width / (samples_i.size - 1)
            area_i = float(np.trapezoid(samples_i, dx=dt))
        else:
            area_i = float(samples_i[0] * width)
    scale = theta / area_i if abs(area_i) > 1e-12 else 0.0
    samples_i = samples_i * scale
    samples_q = samples_q * scale
    env_i = resample_to_length(samples_i, width, seg_slots)
    env_q = resample_to_length(samples_q, width, seg_slots)
    if seg_time > 0.0:
        if seg_slots > 1:
            dt_seg = seg_time / (seg_slots - 1)
            area_seg = float(np.trapezoid(env_i, dx=dt_seg))
        else:
            area_seg = float(env_i[0] * seg_time)
        if abs(area_seg) > 1e-12:
            correction = theta / area_seg
            env_i *= correction
            env_q *= correction
        else:
            value = theta / max(seg_time, 1e-12)
            env_i[:] = value
            env_q[:] = 0.0
    else:
        env_i[:] = 0.0
        env_q[:] = 0.0
    return env_i, env_q


def make_iswap_envelope(spec: PulseSpec, width: float, seg_time: float, seg_slots: int) -> np.ndarray:
    theta = spec.theta
    if seg_slots <= 0:
        return np.zeros(0, dtype=float)
    samples = extract_samples(spec.entry, "samples_i")
    if samples.size == 0:
        samples = np.ones(2, dtype=float)
    if samples.size > 1:
        dt = width / (samples.size - 1)
        area = float(np.trapezoid(samples, dx=dt))
    else:
        area = float(samples[0] * width)
    if abs(area) < 1e-12:
        samples = np.ones_like(samples, dtype=float)
        if samples.size > 1:
            dt = width / (samples.size - 1)
            area = float(np.trapezoid(samples, dx=dt))
        else:
            area = float(samples[0] * width)
    scale = theta / area if abs(area) > 1e-12 else 0.0
    samples = samples * scale
    env = resample_to_length(samples, width, seg_slots)
    if seg_time > 0.0:
        if seg_slots > 1:
            dt_seg = seg_time / (seg_slots - 1)
            area_seg = float(np.trapezoid(env, dx=dt_seg))
        else:
            area_seg = float(env[0] * seg_time)
        if abs(area_seg) > 1e-12:
            env *= theta / area_seg
        else:
            env[:] = theta / max(seg_time, 1e-12)
    else:
        env[:] = 0.0
    return env


def make_library_guess(seq: List[PulseSpec], basis: ControlBasis, n_ts: int, evo_time: float) -> np.ndarray:
    n_ctrls = len(basis.ctrls)
    guess = np.zeros((n_ts, n_ctrls), dtype=float)
    if not seq:
        return guess
    widths = [max(spec.width, 1e-12) for spec in seq]
    total_width = float(sum(widths))
    if total_width <= 0.0:
        total_width = float(len(seq))
        widths = [total_width / len(seq)] * len(seq)
    slots_remaining = n_ts
    cursor = 0
    for idx, spec in enumerate(seq):
        width = widths[idx]
        if idx == len(seq) - 1:
            seg_slots = slots_remaining
        else:
            tentative = int(round(width / total_width * n_ts))
            seg_slots = max(1, tentative)
            seg_slots = min(seg_slots, slots_remaining - max(len(seq) - idx - 1, 0))
        seg_slots = max(1, seg_slots)
        seg_time = evo_time * (seg_slots / max(n_ts, 1))
        sl = slice(cursor, cursor + seg_slots)
        lq = local_qubits(spec)
        if spec.gate == "rx" and len(lq) == 1:
            q = lq[0]
            ix_idx = basis.ix_index.get(q)
            iy_idx = basis.iy_index.get(q)
            env_i, env_q = make_rx_envelopes(spec, width, seg_time, seg_slots)
            if ix_idx is not None and env_i.size == seg_slots:
                guess[sl, ix_idx] += env_i
            if iy_idx is not None and env_q.size == seg_slots:
                guess[sl, iy_idx] += env_q
        elif spec.gate == "iswap" and len(lq) == 2:
            pair = tuple(sorted(lq))
            ctrl_idx = basis.coupling_index.get(pair)
            env = make_iswap_envelope(spec, width, seg_time, seg_slots)
            if ctrl_idx is not None and env.size == seg_slots:
                guess[sl, ctrl_idx] += env
        cursor += seg_slots
        slots_remaining -= seg_slots
    if cursor < n_ts:
        guess[cursor:, :] = guess[cursor - 1, :] if cursor > 0 else 0.0
    return guess


def save_pulse_plot(
    plot_path: Path,
    final_amps_real: np.ndarray,
    labels: Sequence[str],
    dt: float,
) -> Optional[Path]:
    if final_amps_real.size == 0:
        return None
    n_ts, n_ctrls = final_amps_real.shape
    if n_ts == 0 or n_ctrls == 0:
        return None
    times = np.arange(n_ts + 1, dtype=float) * dt
    times_ns = times * 1e9
    amps_step = np.vstack([final_amps_real, final_amps_real[-1:, :]])

    fig, ax = plt.subplots(figsize=(10, 6))
    for idx in range(n_ctrls):
        label = labels[idx] if idx < len(labels) else f"ctrl_{idx}"
        ax.step(times_ns, amps_step[:, idx], where="post", label=label)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Control amplitude")
    ax.set_title("Optimized merged pulse controls")
    ax.grid(True, which="both", linestyle="--", alpha=0.3)
    ax.legend(loc="upper right", fontsize="small")
    fig.tight_layout()
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot_path, dpi=200)
    plt.close(fig)
    return plot_path


def infer_dt_from_library(seq: List[PulseSpec]) -> float:
    candidates: List[float] = []
    for spec in seq:
        width = spec.width
        if width <= 0.0:
            continue
        entry = spec.entry
        for key in ("samples_i", "samples_q", "samples"):
            samples = entry.get(key) or []
            if isinstance(samples, (list, tuple)) and len(samples) > 1:
                dt = width / float(len(samples) - 1)
                if dt > 0.0:
                    candidates.append(dt)
    if candidates:
        return min(candidates)
    widths = [spec.width for spec in seq if spec.width > 0.0]
    if widths:
        return min(widths) / 20.0
    return 1e-9


def build_native_pulse_doc(seq: List[PulseSpec]) -> dict:
    pulse_library = []
    schedule = []
    current_start = 0.0
    for idx, spec in enumerate(seq):
        entry = spec.entry
        pulse_id = entry.get("id", f"pulse_{idx}")
        new_id = f"native_{idx}_{pulse_id}"
        entry_copy = copy.deepcopy(entry)
        entry_copy["id"] = new_id
        entry_copy["qubits"] = list(local_qubits(spec))
        pulse_library.append(entry_copy)
        duration = float(spec.width)
        schedule.append(
            {
                "index": idx,
                "gate": entry_copy.get("gate"),
                "qubits": list(local_qubits(spec)),
                "pulse_id": new_id,
                "start_time": current_start,
                "duration": duration,
                "parameters": entry_copy.get("parameters", {}),
                "virtual": bool(entry_copy.get("virtual")),
            }
        )
        current_start += duration
    return {"pulse_library": pulse_library, "schedule": schedule}


def simulate_native_unitary(seq: List[PulseSpec], nq: int) -> qt.Qobj:
    pulse_doc = build_native_pulse_doc(seq)
    edges: List[Tuple[int, int]] = []
    for spec in seq:
        if spec.gate == "iswap" and len(local_qubits(spec)) == 2:
            pair = tuple(sorted(local_qubits(spec)))
            if pair not in edges:
                edges.append(pair)

    columns: List[np.ndarray] = []
    for basis_index in range(2 ** nq):
        bits = [(basis_index >> (nq - 1 - k)) & 1 for k in range(nq)]
        ket = qt.tensor([qt.basis(2, bit) for bit in bits])
        (
            _,
            states,
            _t_eval,
            _dt,
            _phys_events,
            _virt_events,
            _final_phase,
            residual_phase,
            _event_records,
        ) = simulate_pulse_schedule(pulse_doc, nq, edges, rho_initial=ket)
        if states:
            state_final = states[-1]
        else:
            state_final = ket
        if not state_final.isket:
            evals, evecs = state_final.eigenstates()
            idx = int(np.argmax(np.real(evals)))
            state_final = evecs[idx]
        if residual_phase:
            U_res = phase_unitary(residual_phase, nq)
            state_final = U_res * state_final
        columns.append(state_final.full())
    U_native = np.hstack(columns)
    return qt.Qobj(U_native, dims=[[2] * nq, [2] * nq])


def optimize_candidate_sequence(
    candidate: MergeCandidate,
    seq: List[PulseSpec],
    unique_qubits: List[int],
    logical_to_physical: Dict[int, int],
    settings: OptimizationSettings,
    plot_dir: Optional[Path] = None,
) -> OptimizationResult:
    unique_qubits = list(unique_qubits)
    nq = len(unique_qubits)
    if nq == 0 or nq > 2:
        raise ValueError(f"unsupported number of qubits ({nq}) in sequence for candidate {candidate}")

    U_targ = ideal_from_sequence(seq, nq)
    U_native = simulate_native_unitary(seq, nq)
    native_fidelity = unitary_fidelity(U_targ, U_native)

    total_width = float(sum(spec.width for spec in seq))
    dt_library = infer_dt_from_library(seq)
    native_time = max(total_width, dt_library)

    if settings.evo_time_ns is not None and settings.evo_time_ns > 0.0:
        requested_time = float(settings.evo_time_ns) * 1e-9
    else:
        ratio = 0.5 if settings.compression_ratio is None else float(settings.compression_ratio)
        if ratio <= 0.0:
            raise ValueError("compression_ratio must be positive")
        requested_time = ratio * native_time
        requested_time += float(settings.compression_shift_ns) * 1e-9

    requested_time = max(dt_library, requested_time)
    n_ts = max(1, int(math.ceil(requested_time / dt_library)))
    evo_time = n_ts * dt_library
    dt_effective = dt_library
    requested_ratio = requested_time / native_time if native_time > 0.0 else 1.0
    actual_ratio = evo_time / native_time if native_time > 0.0 else 1.0

    dims = [[2] * nq, [2] * nq]
    H_d = qt.Qobj(np.zeros((2 ** nq, 2 ** nq), dtype=complex), dims=dims)
    basis = build_controls(nq, logical_to_physical=logical_to_physical)
    n_ctrls = len(basis.ctrls)

    library_guess = make_library_guess(seq, basis, n_ts, evo_time)
    guess = library_guess.copy() if settings.seed_from_library else None

    default_amp = 5.0
    amp_lbounds: List[float] = []
    amp_ubounds: List[float] = []
    user_amp = None if settings.amp_bound is None else abs(settings.amp_bound)
    guess_peak = (
        np.max(np.abs(library_guess), axis=0) if library_guess.size else np.zeros(n_ctrls)
    )
    for idx, label in enumerate(basis.labels):
        peak = float(guess_peak[idx]) if idx < guess_peak.size else 0.0
        if user_amp is not None and user_amp > 0.0:
            upper = user_amp
        else:
            upper = max(default_amp, 1.05 * peak)
        lower = 0.0 if label.startswith("J") else -upper
        amp_lbounds.append(lower)
        amp_ubounds.append(upper)

    ctrl_scale = np.array([max(1.0, abs(ub)) for ub in amp_ubounds], dtype=float)
    scaled_ctrls = [
        basis.ctrls[idx] * ctrl_scale[idx] for idx in range(n_ctrls)
    ]
    amp_lbounds_scaled = (np.asarray(amp_lbounds) / ctrl_scale).tolist()
    amp_ubounds_scaled = (np.asarray(amp_ubounds) / ctrl_scale).tolist()

    if guess is not None and guess.size:
        lower_arr = np.asarray(amp_lbounds_scaled)[np.newaxis, :]
        upper_arr = np.asarray(amp_ubounds_scaled)[np.newaxis, :]
        guess = np.clip(guess / ctrl_scale[np.newaxis, :], lower_arr, upper_arr)

    p_type = "LIN"
    optim = pulseoptim.create_pulse_optimizer(
        H_d,
        scaled_ctrls,
        qt.qeye([2] * nq),
        U_targ,
        num_tslots=n_ts,
        evo_time=evo_time,
        amp_lbound=amp_lbounds_scaled,
        amp_ubound=amp_ubounds_scaled,
        fid_err_targ=1e-4,
        min_grad=1e-20,
        max_iter=settings.max_iter,
        max_wall_time=settings.max_wall_time,
        alg="GRAPE",
        optim_method="FMIN_L_BFGS_B",
        method_params={"max_metric_corr": 20, "accuracy_factor": 1e8},
        dyn_type="UNIT",
        fid_params={"phase_option": "PSU"},
        log_level=0,
        gen_stats=False,
        init_pulse_type=p_type,
    )

    dyn = optim.dynamics
    dyn.init_timeslots()
    p_gen = optim.pulse_generator
    n_ctrls = len(basis.ctrls)
    init_amps = np.zeros((n_ts, n_ctrls), dtype=float)
    if isinstance(p_gen, pulsegen.PulseGenLinear):
        for j in range(n_ctrls):
            p_gen.lbound = amp_lbounds_scaled[j]
            p_gen.ubound = amp_ubounds_scaled[j]
            p_gen.scaling = float(j) - float(n_ctrls - 1) / 2.0
            init_amps[:, j] = p_gen.gen_pulse()
    elif getattr(p_gen, "periodic", False):
        phase_diff = math.pi / max(n_ctrls, 1)
        for j in range(n_ctrls):
            p_gen.lbound = amp_lbounds_scaled[j]
            p_gen.ubound = amp_ubounds_scaled[j]
            init_amps[:, j] = p_gen.gen_pulse(start_phase=phase_diff * j)
    else:
        for j in range(n_ctrls):
            p_gen.lbound = amp_lbounds_scaled[j]
            p_gen.ubound = amp_ubounds_scaled[j]
            init_amps[:, j] = p_gen.gen_pulse()

    if guess is not None and guess.size:
        lower_arr = np.asarray(amp_lbounds_scaled)[np.newaxis, :]
        upper_arr = np.asarray(amp_ubounds_scaled)[np.newaxis, :]
        guess = np.clip(guess, lower_arr, upper_arr)
        init_amps = 0.5 * init_amps + 0.5 * guess

    lower_bounds = np.repeat(np.asarray(amp_lbounds_scaled)[np.newaxis, :], n_ts, axis=0)
    upper_bounds = np.repeat(np.asarray(amp_ubounds_scaled)[np.newaxis, :], n_ts, axis=0)
    init_amps = np.clip(init_amps, lower_bounds, upper_bounds)

    dyn.initialize_controls(init_amps.copy())

    res = optim.run_optimization()
    optimized_fidelity = max(0.0, 1.0 - float(np.real(res.fid_err)))

    final_amps = np.array(res.final_amps)
    final_amps_physical = final_amps * ctrl_scale[np.newaxis, :]
    final_amps_real = np.real(final_amps_physical)

    control_waveforms: Dict[str, np.ndarray] = {}
    for idx, label in enumerate(basis.labels):
        control_waveforms[label] = final_amps_real[:, idx].copy()

    plot_path: Optional[Path] = None
    if plot_dir is not None:
        try:
            plot_dir = Path(plot_dir)
            plot_dir.mkdir(parents=True, exist_ok=True)
            qubit_suffix = "_".join(str(q) for q in unique_qubits)
            plot_filename = f"{candidate.gate}_{qubit_suffix}_merged_pulse.png"
            plot_candidate = plot_dir / plot_filename
            saved = save_pulse_plot(
                plot_candidate,
                final_amps_real,
                basis.labels,
                dt_effective,
            )
            if saved is not None:
                plot_path = saved
        except Exception as exc:
            print(f"Warning: failed to save pulse plot for {candidate} ({exc})", file=sys.stderr)

    return OptimizationResult(
        candidate=candidate,
        seq=seq,
        unique_qubits=unique_qubits,
        logical_to_physical=logical_to_physical,
        basis=basis,
        control_waveforms=control_waveforms,
        evo_time=float(evo_time),
        dt_effective=float(dt_effective),
        native_time=float(native_time),
        native_fidelity=float(native_fidelity),
        optimized_fidelity=float(optimized_fidelity),
        requested_ratio=float(requested_ratio),
        actual_ratio=float(actual_ratio),
        final_fid_err=float(np.real(res.fid_err)),
        n_ts=int(n_ts),
        plot_path=plot_path,
    )


def gate_sequence_descriptor(seq: Sequence[PulseSpec]) -> List[str]:
    descriptor: List[str] = []
    for spec in seq:
        gate = spec.gate
        if not descriptor or descriptor[-1] != gate:
            descriptor.append(gate)
    return descriptor


def make_gate_label(result: OptimizationResult) -> str:
    raw_name = result.candidate.gate.lower()
    gate_part = re.sub(r"[^a-z0-9]+", "_", raw_name).strip("_")
    if not gate_part:
        gate_part = "gate"
    return f"merged_{gate_part}"


def ensure_pulse_list(pulse_doc: dict) -> Tuple[List[dict], str]:
    for key in ("pulse_definitions", "pulse_library"):
        value = pulse_doc.get(key)
        if isinstance(value, list):
            return value, key
    pulse_doc["pulse_definitions"] = []
    return pulse_doc["pulse_definitions"], "pulse_definitions"


def extract_physical_qubits_from_device(device_doc: dict) -> List[int]:
    qubits: set[int] = set()
    physical = device_doc.get("physical_qubits")
    if isinstance(physical, Sequence):
        for q in physical:
            try:
                qubits.add(int(q))
            except Exception:
                continue
    metadata = device_doc.get("metadata")
    if isinstance(metadata, dict):
        meta_phys = metadata.get("physical_qubits")
        if isinstance(meta_phys, Sequence):
            for q in meta_phys:
                try:
                    qubits.add(int(q))
                except Exception:
                    continue
    return sorted(qubits)


def extract_coupling_pairs_from_device(device_doc: dict) -> List[Tuple[int, int]]:
    pairs: set[Tuple[int, int]] = set()

    def ingest(source) -> None:
        if isinstance(source, Sequence):
            for item in source:
                if isinstance(item, Sequence) and len(item) == 2:
                    try:
                        q0 = int(item[0])
                        q1 = int(item[1])
                    except Exception:
                        continue
                    ordered = tuple(sorted((q0, q1)))
                    if ordered[0] != ordered[1]:
                        pairs.add(ordered)

    ingest(device_doc.get("cx_coupling"))
    ingest(device_doc.get("coupling_map"))
    metadata = device_doc.get("metadata")
    if isinstance(metadata, dict):
        ingest(metadata.get("cx_coupling"))
        ingest(metadata.get("coupling_map"))
    return sorted(pairs)


def extract_candidate_theta(candidate: MergeCandidate, seq: Sequence[PulseSpec]) -> Optional[float]:
    for spec in seq:
        if spec.gate == candidate.gate:
            return float(spec.theta)
    return None


def synthesize_pulse_entries(result: OptimizationResult, gate_label: str) -> List[dict]:
    width = float(result.evo_time)
    dt = float(result.dt_effective)
    sequence_descriptor = gate_sequence_descriptor(result.seq)
    theta_value = extract_candidate_theta(result.candidate, result.seq)

    drive_waveforms: Dict[int, Dict[str, np.ndarray]] = {}
    coupling_waveforms: Dict[Tuple[int, int], np.ndarray] = {}
    for label, waveform in result.control_waveforms.items():
        axis, qubits = parse_control_label(label)
        if axis in ("I", "Q") and len(qubits) == 1:
            q = qubits[0]
            drive_data = drive_waveforms.setdefault(q, {})
            drive_data[axis] = waveform
        elif axis == "J" and len(qubits) == 2:
            pair = tuple(sorted(qubits))
            coupling_waveforms[pair] = waveform

    entries: List[dict] = []
    for qubit in sorted(drive_waveforms):
        axes = drive_waveforms[qubit]
        samples_i = axes.get("I")
        samples_q = axes.get("Q")
        n_samples = 0
        if samples_i is not None:
            n_samples = samples_i.size
        elif samples_q is not None:
            n_samples = samples_q.size
        if n_samples == 0:
            continue
        if samples_i is None:
            samples_i = np.zeros(n_samples, dtype=float)
        if samples_q is None:
            samples_q = np.zeros(n_samples, dtype=float)
        samples_i_list = [float(v) for v in samples_i.tolist()]
        samples_q_list = [float(v) for v in samples_q.tolist()]
        params = {
            "sequence": sequence_descriptor,
        }
        if theta_value is not None:
            params["theta_rad"] = theta_value
        entry = {
            "id": f"{gate_label}_drive_q{qubit}",
            "gate": gate_label,
            "qubits": [int(qubit)],
            "waveform_type": "arbitrary",
            "width": width,
            "parameters": params,
            "samples_i": samples_i_list,
            "samples_q": samples_q_list,
            "amplitude": float(max(np.max(np.abs(samples_i)), np.max(np.abs(samples_q)))),
        }
        entries.append(entry)

    for pair in sorted(coupling_waveforms):
        samples = coupling_waveforms[pair]
        if samples.size == 0:
            continue
        n_samples = int(samples.size)
        samples_i_list = [float(v) for v in samples.tolist()]
        samples_q_list = [0.0 for _ in range(n_samples)]
        params = {
            "sequence": sequence_descriptor,
        }
        if theta_value is not None:
            params["theta_rad"] = theta_value
        q0, q1 = pair
        base_entry = {
            "id": f"{gate_label}_coupling_q{q0}_q{q1}",
            "gate": gate_label,
            "qubits": [int(q0), int(q1)],
            "waveform_type": "arbitrary",
            "width": width,
            "parameters": params,
            "samples_i": samples_i_list,
            "samples_q": samples_q_list,
            "amplitude": float(np.max(np.abs(samples))),
        }
        entries.append(base_entry)
        if q0 != q1:
            reversed_entry = dict(base_entry)
            reversed_entry["id"] = f"{gate_label}_coupling_q{q1}_q{q0}"
            reversed_entry["qubits"] = [int(q1), int(q0)]
            entries.append(reversed_entry)

    return entries


def append_pulse_entries(pulse_doc: dict, entries: Sequence[dict]) -> None:
    library, key = ensure_pulse_list(pulse_doc)
    new_ids = {entry.get("id") for entry in entries if "id" in entry}
    if new_ids:
        library[:] = [entry for entry in library if entry.get("id") not in new_ids]
    library.extend(entries)
    if key == "pulse_library":
        pulse_doc["pulse_library"] = library
    else:
        pulse_doc["pulse_definitions"] = library


def replicate_gate_entries(
    pulse_doc: dict,
    gate_label: str,
    target_single_qubits: Sequence[int],
    target_pairs: Sequence[Tuple[int, int]],
) -> None:
    library, key = ensure_pulse_list(pulse_doc)
    existing = [entry for entry in library if entry.get("gate") == gate_label]
    base_drive = next((entry for entry in existing if len(entry.get("qubits", [])) == 1), None)
    base_coupling = next((entry for entry in existing if len(entry.get("qubits", [])) == 2), None)

    if base_drive and target_single_qubits:
        for q in target_single_qubits:
            qubits = [int(q)]
            if any(entry.get("gate") == gate_label and entry.get("qubits") == qubits for entry in library):
                continue
            new_entry = copy.deepcopy(base_drive)
            new_entry["id"] = f"{gate_label}_drive_q{q}"
            new_entry["qubits"] = qubits
            library.append(new_entry)

    if base_coupling and target_pairs:
        for pair in target_pairs:
            q0, q1 = pair
            for ordered in ((q0, q1), (q1, q0)):
                qubits = [int(ordered[0]), int(ordered[1])]
                if any(entry.get("gate") == gate_label and entry.get("qubits") == qubits for entry in library):
                    continue
                new_entry = copy.deepcopy(base_coupling)
                new_entry["id"] = f"{gate_label}_coupling_q{ordered[0]}_q{ordered[1]}"
                new_entry["qubits"] = qubits
                library.append(new_entry)

    if key == "pulse_library":
        pulse_doc["pulse_library"] = library
    else:
        pulse_doc["pulse_definitions"] = library


def update_device_config(
    device_doc: dict,
    gate_label: str,
    duration_s: float,
    error_rate: float,
    logical_gate: Optional[str] = None,
    qubits: Optional[Sequence[int]] = None,
) -> None:
    duration_s = float(duration_s)
    error_rate = float(max(0.0, error_rate))

    def ensure_mapping(doc: dict, key: str) -> Dict[str, float]:
        mapping = doc.get(key)
        if not isinstance(mapping, dict):
            mapping = {}
            doc[key] = mapping
        return mapping

    def ensure_list(doc: dict, key: str) -> List[str]:
        value = doc.get(key)
        if not isinstance(value, list):
            value = []
            doc[key] = value
        return value

    gate_lens = ensure_mapping(device_doc, "gate_lens")
    gate_lens[gate_label] = duration_s
    gate_errs = ensure_mapping(device_doc, "gate_errs")
    gate_errs[gate_label] = error_rate
    basis_gates = ensure_list(device_doc, "basis_gates")
    if gate_label not in basis_gates:
        basis_gates.append(gate_label)

    metadata = device_doc.get("metadata")
    if isinstance(metadata, dict):
        meta_basis = metadata.get("basis_gates")
        if isinstance(meta_basis, list) and gate_label not in meta_basis:
            meta_basis.append(gate_label)
        meta_gate_lens = metadata.get("gate_lens")
        if isinstance(meta_gate_lens, dict):
            meta_gate_lens[gate_label] = duration_s
        meta_gate_errs = metadata.get("gate_errs")
        if isinstance(meta_gate_errs, dict):
            meta_gate_errs[gate_label] = error_rate

    if logical_gate:
        logical_gate_l = str(logical_gate).strip().lower()
        if logical_gate_l:
            alias_entry = {
                "logical_gate": logical_gate_l,
            }
            alias_map = ensure_mapping(device_doc, "merged_gate_aliases")
            alias_map[gate_label] = alias_entry
            metadata = device_doc.get("metadata")
            if isinstance(metadata, dict):
                meta_alias = metadata.get("merged_gate_aliases")
                if not isinstance(meta_alias, dict):
                    meta_alias = {}
                    metadata["merged_gate_aliases"] = meta_alias
                meta_alias[gate_label] = dict(alias_entry)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Optimize merged pulses and augment pulse/device configuration files.",
    )
    parser.add_argument(
        "--pulse-dump",
        type=Path,
        required=True,
        help="Pulse JSON emitted by qasmtrans containing the schedule, pulse library, and critical path metadata.",
    )
    parser.add_argument(
        "--pulse-template",
        type=Path,
        required=True,
        help="Base pulse template JSON to clone and augment with merged definitions.",
    )
    parser.add_argument(
        "--device-config",
        type=Path,
        required=True,
        help="Base device configuration JSON to clone and augment with merged gate timing/error data.",
    )
    parser.add_argument(
        "--output-pulses",
        type=Path,
        required=True,
        help="Destination path for the augmented pulse template JSON.",
    )
    parser.add_argument(
        "--output-device",
        type=Path,
        required=True,
        help="Destination path for the augmented device configuration JSON.",
    )
    parser.add_argument(
        "--merge-limit",
        type=int,
        default=0,
        help="Maximum number of merge candidates to optimize (0 means no limit).",
    )
    parser.add_argument(
        "--compression-ratio",
        type=float,
        default=0.5,
        help="Fraction of the native gate duration to target (default 0.5).",
    )
    parser.add_argument(
        "--compression-shift-ns",
        type=float,
        default=0.0,
        help="Additional time shift (ns) applied after scaling by the compression ratio.",
    )
    parser.add_argument(
        "--evo-time-ns",
        type=float,
        default=None,
        help="Explicit merged gate duration in ns (overrides compression parameters).",
    )
    parser.add_argument(
        "--amp-bound",
        type=float,
        default=None,
        help="Optional absolute amplitude bound applied to all controls (auto if omitted).",
    )
    parser.add_argument(
        "--seed-from-library",
        action="store_true",
        help="Initialise the optimisation with the stitched library waveform instead of the default guess.",
    )
    parser.add_argument(
        "--min-merged-fidelity",
        type=float,
        default=0.95,
        help="Skip merged gates whose optimised fidelity falls below this value.",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        default=None,
        help="If provided, save control amplitude plots for each merged gate in this directory.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-candidate optimisation progress.",
    )
    args = parser.parse_args()

    pulse_dump_doc = load_json(args.pulse_dump)
    pulse_lib = canonicalize_pulse_library(pulse_dump_doc)
    library_by_id: Dict[str, dict] = {}
    for entry in pulse_lib:
        pulse_id = entry.get("id")
        if isinstance(pulse_id, str) and pulse_id:
            library_by_id[pulse_id] = entry
    schedule = pulse_dump_doc.get("schedule") or []
    device_template_doc = load_json(args.device_config)
    native_gates = extract_native_gates(device_template_doc)
    critical_indices = {
        int(idx) for idx in (pulse_dump_doc.get("critical_path") or {}).get("path_indices", [])
    }
    if critical_indices:
        schedule_path = [
            event for event in schedule if int(event.get("index", -1)) in critical_indices
        ]
    else:
        schedule_path = schedule

    allow_parameterized = pulse_dump_doc.get("allow_parameterized_candidates")
    if allow_parameterized is None:
        allow_parameterized = (
            (pulse_dump_doc.get("critical_path") or {}).get("allow_parameterized_candidates", True)
        )
    allow_parameterized = bool(allow_parameterized)

    def label_has_parameters(label: str) -> bool:
        if allow_parameterized:
            return False
        if not label:
            return False
        if "[" not in label:
            return False
        open_bracket = label.find("[")
        close_bracket = label.find("]", open_bracket + 1)
        if close_bracket == -1:
            return False
        return "=" in label[open_bracket:close_bracket]

    logical_instances: Dict[int, dict] = {}
    for event in schedule:
        logical_id = event.get("logical_gate_id")
        if logical_id is None:
            continue
        try:
            logical_id_int = int(logical_id)
        except (TypeError, ValueError):
            continue
        raw_label = str(event.get("logical_label", ""))
        if not raw_label or label_has_parameters(raw_label):
            continue
        label_lower = raw_label.lower()
        info = logical_instances.setdefault(
            logical_id_int,
            {
                "label": label_lower,
                "events": [],
                "qubits": set(),
                "ordered_qubits": [],
                "duration": 0.0,
                "first_index": None,
            },
        )
        qubits_list = [int(q) for q in (event.get("qubits") or [])]
        if not qubits_list:
            continue
        info["events"].append(event)
        for q in qubits_list:
            info["qubits"].add(int(q))
        if not info["ordered_qubits"] or len(qubits_list) > len(info["ordered_qubits"]):
            info["ordered_qubits"] = list(qubits_list)
        info["duration"] += float(event.get("duration", 0.0) or 0.0)
        if info["first_index"] is None:
            try:
                info["first_index"] = int(event.get("index"))
            except Exception:
                info["first_index"] = len(logical_instances)

    gate_target_qubits: Dict[str, set[Tuple[int, ...]]] = {}
    instance_scores: List[Tuple[str, Tuple[int, ...], int, float]] = []
    for logical_id, info in logical_instances.items():
        if not info["events"]:
            continue
        label = info["label"]
        if not label:
            continue
        if label_has_parameters(label):
            continue
        qubits_tuple = tuple(sorted(info["qubits"]))
        order = info.get("first_index")
        if order is None:
            order = 0.0
        instance_scores.append((label, qubits_tuple, logical_id, float(order)))
        gate_target_qubits.setdefault(label, set()).add(qubits_tuple)
    instance_scores.sort(key=lambda item: item[3])

    schedule_totals = aggregate_schedule_candidates(schedule_path)

    rx_specs = collect_single_qubit_specs(pulse_lib)
    iswap_specs = collect_two_qubit_specs(pulse_lib)
    top_gate_entries = extract_top_gate_entries(pulse_dump_doc)
    schedule_single_qubits = sorted(
        {
            int(tpl[0])
            for tuples in gate_target_qubits.values()
            for tpl in tuples
            if len(tpl) == 1
        }
    )
    schedule_pair_qubits = sorted(
        {
            tuple(sorted((int(tpl[0]), int(tpl[1]))))
            for tuples in gate_target_qubits.values()
            for tpl in tuples
            if len(tpl) == 2
        }
    )
    candidates = select_merge_candidates(
        rx_specs,
        iswap_specs,
        args.merge_limit,
        native_gates=native_gates,
        schedule_totals=schedule_totals,
        instance_scores=instance_scores,
        top_gate_entries=top_gate_entries,
        logical_instances=logical_instances,
    )
    if not candidates:
        print("No merge candidates matched the supplied pulse dump.", file=sys.stderr)
        return 1

    settings = OptimizationSettings(
        compression_ratio=args.compression_ratio,
        compression_shift_ns=args.compression_shift_ns,
        evo_time_ns=args.evo_time_ns,
        amp_bound=args.amp_bound,
        seed_from_library=args.seed_from_library,
        min_merged_fidelity=args.min_merged_fidelity,
    )

    results: List[OptimizationResult] = []
    for candidate in candidates:
        seq_info = build_sequence_for_candidate(
            candidate,
            rx_specs,
            iswap_specs,
            logical_instances=logical_instances,
            library_by_id=library_by_id,
        )
        if seq_info is None:
            if args.verbose:
                print(f"Skipping {candidate}: missing supporting pulses")
            continue
        seq, unique_qubits, logical_to_physical = seq_info
        try:
            if args.verbose:
                qubit_str = ", ".join(str(q) for q in unique_qubits)
                print(
                    f"Optimizing {candidate.gate} on qubits [{qubit_str}] "
                    f"(score={candidate.score:.6g})"
                )
            result = optimize_candidate_sequence(
                candidate,
                seq,
                unique_qubits,
                logical_to_physical,
                settings,
                args.plot_dir,
            )
            if result.optimized_fidelity < settings.min_merged_fidelity:
                if args.verbose:
                    print(
                        f"  merged fidelity {result.optimized_fidelity:.6f} "
                        f"below threshold {settings.min_merged_fidelity:.6f}; skipping"
                    )
                continue
            results.append(result)
            if args.verbose:
                print(
                    f"  fidelity {result.optimized_fidelity:.6f}, "
                    f"duration {result.evo_time * 1e9:.3f} ns "
                    f"(native {result.native_time * 1e9:.3f} ns)"
                )
        except Exception as exc:
            print(f"Failed to optimise candidate {candidate}: {exc}", file=sys.stderr)

    if not results:
        print("No merged pulses were successfully produced.", file=sys.stderr)
        return 1

    pulse_template_doc = load_json(args.pulse_template)
    pulse_output_doc = copy.deepcopy(pulse_template_doc)
    device_output_doc = copy.deepcopy(device_template_doc)

    device_single_targets = extract_physical_qubits_from_device(device_template_doc)
    device_pair_targets = extract_coupling_pairs_from_device(device_template_doc)

    num_qubits = int(
        device_template_doc.get("num_qubits")
        or device_template_doc.get("metadata", {}).get("num_qubits")
        or 0
    )

    summary_rows: List[Tuple[str, List[int], float, float, float, int]] = []
    for result in results:
        gate_label = make_gate_label(result)
        entries = synthesize_pulse_entries(result, gate_label)
        append_pulse_entries(pulse_output_doc, entries)
        gate_key = result.candidate.gate
        target_single_qubits = list(device_single_targets)
        if not target_single_qubits:
            if num_qubits > 0:
                target_single_qubits = list(range(num_qubits))
            else:
                target_single_qubits = list(schedule_single_qubits)
        if not target_single_qubits:
            target_single_qubits = sorted({spec.qubits[0] for spec in rx_specs.values()})

        target_pairs = list(device_pair_targets)
        if not target_pairs:
            if num_qubits > 0:
                target_pairs = [(i, j) for i in range(num_qubits) for j in range(i + 1, num_qubits)]
            else:
                target_pairs = list(schedule_pair_qubits)
        if not target_pairs:
            target_pairs = sorted(iswap_specs.keys())
        if not target_pairs:
            target_pairs = sorted(
                {tuple(sorted(tpl)) for tpl in gate_target_qubits.get(gate_key, set()) if len(tpl) == 2}
            )
        replicate_gate_entries(
            pulse_output_doc,
            gate_label,
            target_single_qubits,
            target_pairs if len(result.unique_qubits) > 1 else [],
        )
        error_rate = 1.0 - result.optimized_fidelity
        update_device_config(
            device_output_doc,
            gate_label,
            result.evo_time,
            error_rate,
            logical_gate=result.candidate.gate,
            qubits=result.unique_qubits,
        )
        summary_rows.append(
            (
                gate_label,
                result.unique_qubits,
                result.native_time,
                result.evo_time,
                result.optimized_fidelity,
                len(entries),
            )
        )

    args.output_pulses.parent.mkdir(parents=True, exist_ok=True)
    with args.output_pulses.open("w", encoding="utf-8") as handle:
        json.dump(pulse_output_doc, handle, indent=2)

    args.output_device.parent.mkdir(parents=True, exist_ok=True)
    with args.output_device.open("w", encoding="utf-8") as handle:
        json.dump(device_output_doc, handle, indent=2)

    print(f"Augmented pulse template written to {args.output_pulses}")
    print(f"Augmented device config written to {args.output_device}")
    print("Merged gates:")
    for gate_label, qubits, native_time, merged_time, fidelity, entry_count in summary_rows:
        qubit_text = ", ".join(str(q) for q in qubits)
        native_ns = native_time * 1e9
        merged_ns = merged_time * 1e9
        delta_ns = merged_ns - native_ns
        delta_pct = 0.0
        if native_ns > 1e-12:
            delta_pct = (delta_ns / native_ns) * 100.0
        print(
            f"  {gate_label}: qubits [{qubit_text}], "
            f"duration {merged_ns:.3f} ns (native {native_ns:.3f} ns, "
            f"delta {delta_ns:+.3f} ns, {delta_pct:+.1f}%), "
            f"fidelity {fidelity:.6f}, entries {entry_count}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
