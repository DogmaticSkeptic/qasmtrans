#!/usr/bin/env python3
"""Simulate Rigetti Ankaa pulse schedules and report circuit fidelity."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import qutip as qt


# --------------------------------------------------------------------------- #
# Utility helpers
# --------------------------------------------------------------------------- #


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def normalize_gate_label(gate: str) -> Tuple[str, float | None]:
    gate_l = gate.lower()
    if gate_l == "x":
        return "rx", math.pi
    if gate_l == "sx":
        return "rx", 0.5 * math.pi
    return gate_l, None


def canonicalize_pulse_doc(pulse_doc: dict) -> None:
    pulse_library = pulse_doc.get("pulse_library")
    if isinstance(pulse_library, list):
        for entry in pulse_library:
            gate_raw = entry.get("gate")
            gate_norm, default_theta = normalize_gate_label(str(gate_raw or ""))
            entry["gate"] = gate_norm
            if default_theta is not None:
                params = entry.setdefault("parameters", {})
                params.setdefault("theta", default_theta)

    schedule = pulse_doc.get("schedule")
    if isinstance(schedule, list):
        for event in schedule:
            gate_raw = event.get("gate")
            gate_norm, default_theta = normalize_gate_label(str(gate_raw or ""))
            event["gate"] = gate_norm
            if default_theta is not None:
                params = event.setdefault("parameters", {})
                params.setdefault("theta", default_theta)

    critical = pulse_doc.get("critical_path")
    if isinstance(critical, dict):
        top_gates = critical.get("top_gates")
        if isinstance(top_gates, list):
            for item in top_gates:
                gate_raw = item.get("gate")
                gate_norm, _ = normalize_gate_label(str(gate_raw or ""))
                item["gate"] = gate_norm
def parse_edges(device_metadata: dict, schedule: Sequence[dict], max_qubit: int) -> List[Tuple[int, int]]:
    edges: set[Tuple[int, int]] = set()
    coupling = device_metadata.get("coupling") or []
    for entry in coupling:
        try:
            a_str, b_str = entry.split("_")
            a = int(a_str)
            b = int(b_str)
            if a < max_qubit and b < max_qubit:
                edges.add(tuple(sorted((a, b))))
        except ValueError:
            continue
    if not edges:
        for item in schedule:
            qubits = item.get("qubits") or []
            if len(qubits) == 2:
                qq = tuple(sorted(int(q) for q in qubits))
                if qq[0] < max_qubit and qq[1] < max_qubit:
                    edges.add(qq)
    return sorted(edges)


def _accumulate_relaxation_map(source, target: Dict[int, float]) -> None:
    if isinstance(source, dict):
        items = source.items()
    elif isinstance(source, (list, tuple)):
        items = enumerate(source)
    else:
        return
    for key, value in items:
        try:
            idx = int(key)
        except (TypeError, ValueError):
            continue
        try:
            val = float(value)
        except (TypeError, ValueError):
            continue
        if val > 0.0:
            target[idx] = val


def extract_t1_t2_times(device_doc: dict) -> Tuple[Dict[int, float], Dict[int, float]]:
    t1_map: Dict[int, float] = {}
    t2_map: Dict[int, float] = {}
    _accumulate_relaxation_map(device_doc.get("T1"), t1_map)
    _accumulate_relaxation_map(device_doc.get("T2"), t2_map)
    metadata = device_doc.get("metadata")
    if isinstance(metadata, dict):
        _accumulate_relaxation_map(metadata.get("T1"), t1_map)
        _accumulate_relaxation_map(metadata.get("T2"), t2_map)
    return t1_map, t2_map


def evaluate_expression(expr: str) -> float:
    expr = expr.strip()
    if not expr:
        return 0.0
    safe_dict = {"pi": math.pi}
    return float(eval(expr, {"__builtins__": None}, safe_dict))


# --------------------------------------------------------------------------- #
# Pulse replay utilities
# --------------------------------------------------------------------------- #


def compute_dt_from_samples(width: float, samples_len: int, params: dict) -> float:
    if samples_len > 1 and width > 0.0:
        return width / float(samples_len - 1)
    dt_ns = params.get("dt_ns")
    if dt_ns is not None:
        return float(dt_ns) * 1e-9
    if width > 0.0:
        # default to 0.25 ns grid if no explicit sample resolution is available
        return min(width, 0.25e-9)
    return 0.25e-9


def resample_waveform(samples: Sequence[float], width: float, dt: float) -> np.ndarray:
    samples = list(samples)
    if not samples:
        n_steps = max(2, int(round(width / dt)) + 1)
        return np.zeros(n_steps, dtype=float)
    if len(samples) == 1:
        n_steps = max(2, int(round(width / dt)) + 1)
        return np.full(n_steps, float(samples[0]), dtype=float)
    raw = np.asarray(samples, dtype=float)
    source_times = np.linspace(0.0, width, raw.size, dtype=float)
    n_steps = max(2, int(round(width / dt)) + 1)
    target_times = np.linspace(0.0, width, n_steps, dtype=float)
    return np.interp(target_times, source_times, raw).astype(float)


def resample_to_length(samples: Sequence[float], width: float, length: int) -> np.ndarray:
    length = max(1, int(length))
    if not samples:
        return np.zeros(length, dtype=float)
    if len(samples) == length:
        return np.asarray(samples, dtype=float)
    if len(samples) == 1:
        return np.full(length, float(samples[0]), dtype=float)
    raw = np.asarray(samples, dtype=float)
    source_times = np.linspace(0.0, width, raw.size, dtype=float)
    if length == 1 or math.isclose(width, 0.0, abs_tol=1e-15):
        return np.full(length, raw[-1], dtype=float)
    target_times = np.linspace(0.0, width, length, dtype=float)
    return np.interp(target_times, source_times, raw).astype(float)


def rotate_iq(i_samples: np.ndarray, q_samples: np.ndarray, phase: float) -> Tuple[np.ndarray, np.ndarray]:
    if abs(phase) <= 1e-12:
        return i_samples, q_samples
    cos_p = math.cos(phase)
    sin_p = math.sin(phase)
    i_rot = cos_p * i_samples - sin_p * q_samples
    q_rot = sin_p * i_samples + cos_p * q_samples
    return i_rot, q_rot


def pulse_rotation_area(pulse: dict, duration: float) -> float:
    """Estimate the integrated rotation angle (radians) represented by a pulse entry."""
    if duration <= 0.0:
        return 0.0
    samples_i = pulse.get("samples_i") or []
    if samples_i:
        arr = np.asarray(samples_i, dtype=float)
        if arr.size >= 2:
            step = duration / max(arr.size - 1, 1)
            area = float(np.trapezoid(arr, dx=step))
        else:
            area = float(arr[0]) * duration
        return area
    amplitude = pulse.get("amplitude")
    if amplitude is not None:
        return float(amplitude) * duration
    return 0.0


def build_hamiltonian(
    time_grid: np.ndarray,
    dt: float,
    I_env: Dict[int, np.ndarray],
    Q_env: Dict[int, np.ndarray],
    J_env: Dict[Tuple[int, int], np.ndarray],
    num_qubits: int,
    edges: Iterable[Tuple[int, int]],
) -> List[List[qt.Qobj]]:
    e0 = qt.basis(2, 0)
    e1 = qt.basis(2, 1)
    s01 = e0 * e1.dag()
    s10 = s01.dag()
    X = s01 + s10
    Y = -1j * s01 + 1j * s10

    def lift_one(op: qt.Qobj, site: int, total: int) -> qt.Qobj:
        ops = [qt.qeye(2)] * total
        ops = list(ops)
        ops[site] = op
        return qt.tensor(ops)

    def lift_two(op_i: qt.Qobj, i: int, op_j: qt.Qobj, j: int, total: int) -> qt.Qobj:
        ops = [qt.qeye(2)] * total
        ops = list(ops)
        ops[i] = op_i
        ops[j] = op_j
        return qt.tensor(ops)

    def coeff_from_array(values: np.ndarray) -> callable:
        arr = np.asarray(values, dtype=float)
        inv_dt = 1.0 / dt

        def coeff(t_val: float, _args=None) -> float:
            idx = int(t_val * inv_dt)
            if idx < 0:
                idx = 0
            elif idx >= arr.size:
                idx = arr.size - 1
            return float(arr[idx])

        return coeff

    H_terms: List[List[qt.Qobj]] = []
    for qubit in range(num_qubits):
        Ii = I_env.get(qubit)
        Qi = Q_env.get(qubit)
        if Ii is not None and np.any(np.abs(Ii) > 0):
            H_terms.append([lift_one(0.5 * X, qubit, num_qubits), coeff_from_array(Ii)])
        if Qi is not None and np.any(np.abs(Qi) > 0):
            H_terms.append([lift_one(0.5 * Y, qubit, num_qubits), coeff_from_array(Qi)])
    for edge in edges:
        J_amp = J_env.get(edge)
        if J_amp is None or not np.any(np.abs(J_amp) > 0):
            continue
        i, j = edge
        exch_xy = 0.5 * (lift_two(X, i, X, j, num_qubits) + lift_two(Y, i, Y, j, num_qubits))
        H_terms.append([exch_xy, coeff_from_array(J_amp)])
    return H_terms


# --------------------------------------------------------------------------- #
# Pulse fidelity sanity checks
# --------------------------------------------------------------------------- #


def ideal_unitary_for_pulse(entry: dict, num_qubits: int) -> qt.Qobj:
    qubits = entry.get("qubits") or []
    if not qubits:
        return qt.qeye([2] * num_qubits)

    gate = str(entry.get("gate", "")).lower()
    waveform_type = str(entry.get("waveform_type", "")).lower()
    is_virtual = bool(entry.get("virtual")) or waveform_type == "virtual"

    if is_virtual:
        theta = (entry.get("parameters") or {}).get("theta")
        if theta is None:
            theta = (entry.get("calibration") or {}).get("default_phase_rad", 0.0)
        phase_map = {int(qubits[0]): float(theta or 0.0)}
        return phase_unitary(phase_map, num_qubits)

    params = entry.get("parameters") or {}
    theta_target = params.get("theta")
    if theta_target is None:
        theta_target = (entry.get("calibration") or {}).get("theta_actual_rad", 0.0)
    theta = float(theta_target)

    if len(qubits) == 1:
        phase_default = float((entry.get("calibration") or {}).get("default_phase_rad", 0.0))
        unitary = single_xy_unitary(theta, phase_default)
        return embed_single(unitary, int(qubits[0]), num_qubits)

    if len(qubits) == 2:
        unitary = two_qubit_xy_unitary(theta)
        return embed_two(unitary, int(qubits[0]), int(qubits[1]), num_qubits)

    return qt.qeye([2] * num_qubits)


def sanity_check_pulse_library(pulse_doc: dict) -> List[dict]:
    results: List[dict] = []
    for entry in pulse_doc.get("pulse_library", []):
        entry_id = entry.get("id", "")
        qubits = entry.get("qubits") or []
        waveform_type = str(entry.get("waveform_type", "")).lower()
        is_virtual = bool(entry.get("virtual")) or waveform_type == "virtual"

        if not qubits:
            results.append({"id": entry_id, "status": "skip", "reason": "no qubits"})
            continue

        num_qubits = max(int(q) for q in qubits) + 1

        if is_virtual:
            results.append({"id": entry_id, "status": "virtual"})
            continue

        width_s = float(entry.get("width", 0.0) or 0.0)
        if width_s <= 0.0:
            duration_ns = (entry.get("calibration") or {}).get("duration_ns")
            if duration_ns is not None:
                width_s = float(duration_ns) * 1e-9

        schedule_entry = {
            "index": 0,
            "gate": entry.get("gate"),
            "qubits": qubits,
            "pulse_id": entry_id,
            "start_time": 0.0,
            "duration": width_s,
            "parameters": entry.get("parameters", {}),
            "virtual": False,
        }

        pulse_doc_single = {
            "pulse_library": [entry],
            "schedule": [schedule_entry],
        }

        edges = []
        if len(qubits) == 2:
            edges = [tuple(sorted(int(q) for q in qubits))]

        rho_init = None
        if len(qubits) == 2:
            basis_states = []
            for idx in range(num_qubits):
                if idx == qubits[0]:
                    basis_states.append(qt.basis(2, 1))
                else:
                    basis_states.append(qt.basis(2, 0))
            psi10 = qt.tensor(basis_states)
            rho_init = qt.ket2dm(psi10)
        (
            rho0,
            states,
            _,
            _,
            _,
            _,
            _,
            residual_map,
            _,
        ) = simulate_pulse_schedule(pulse_doc_single, num_qubits, edges, rho_initial=rho_init)

        rho_sim_state = states[-1] if states else rho0
        rho_sim = qt.ket2dm(rho_sim_state) if rho_sim_state.isket else rho_sim_state
        if any(abs(residual_map.get(q, 0.0)) > 1e-12 for q in range(num_qubits)):
            U_res = phase_unitary(residual_map, num_qubits)
            rho_sim = U_res * rho_sim * U_res.dag()

        unitary = ideal_unitary_for_pulse(entry, num_qubits)
        rho_target = unitary * rho0 * unitary.dag()

        fidelity = float(qt.fidelity(rho_sim, rho_target))
        expected = (entry.get("calibration") or {}).get("expected_fidelity")
        results.append(
            {
                "id": entry_id,
                "status": "physical",
                "fidelity": fidelity,
                "expected": expected,
            }
        )
    return results


# --------------------------------------------------------------------------- #
# Ideal state construction from QASM
# --------------------------------------------------------------------------- #


ISWAP_MATRIX = np.array(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0j, 0.0],
        [0.0, 1.0j, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ],
    dtype=complex,
)


def embed_single(unitary: qt.Qobj, target: int, num_qubits: int) -> qt.Qobj:
    ops = [qt.qeye(2)] * num_qubits
    ops = list(ops)
    ops[target] = unitary
    return qt.tensor(ops)


def embed_two(unitary: np.ndarray, q0: int, q1: int, num_qubits: int) -> qt.Qobj:
    q0, q1 = sorted((q0, q1))
    dim = 2 ** num_qubits
    data = np.zeros((dim, dim), dtype=complex)
    for col in range(dim):
        bits_col = [(col >> (num_qubits - 1 - k)) & 1 for k in range(num_qubits)]
        for row in range(dim):
            bits_row = [(row >> (num_qubits - 1 - k)) & 1 for k in range(num_qubits)]
            match = True
            for idx in range(num_qubits):
                if idx == q0 or idx == q1:
                    continue
                if bits_row[idx] != bits_col[idx]:
                    match = False
                    break
            if not match:
                continue
            sub_col = (bits_col[q0] << 1) | bits_col[q1]
            sub_row = (bits_row[q0] << 1) | bits_row[q1]
            data[row, col] = unitary[sub_row, sub_col]
    return qt.Qobj(data, dims=[[2] * num_qubits, [2] * num_qubits])


def build_collapse_ops(
    num_qubits: int,
    t1_times: Dict[int, float] | None = None,
    t2_times: Dict[int, float] | None = None,
    relax_to_ground: bool = True,
) -> List[qt.Qobj]:
    collapse: List[qt.Qobj] = []
    t1_times = t1_times or {}
    t2_times = t2_times or {}
    for q in range(num_qubits):
        t1 = float(t1_times.get(q, 0.0) or 0.0)
        gamma1 = 0.0
        if t1 > 0.0:
            gamma1 = 1.0 / t1
            if relax_to_ground and gamma1 > 0.0:
                lower = embed_single(qt.sigmam(), q, num_qubits)
                collapse.append(math.sqrt(gamma1) * lower)
        t2 = float(t2_times.get(q, 0.0) or 0.0)
        gamma_phi = 0.0
        if t2 > 0.0:
            gamma2 = 1.0 / t2
            gamma_phi = max(0.0, gamma2 - 0.5 * gamma1)
        if gamma_phi > 0.0:
            sz = embed_single(qt.sigmaz(), q, num_qubits)
            collapse.append(math.sqrt(gamma_phi / 2.0) * sz)
    return collapse


def single_xy_unitary(theta: float, phase: float) -> qt.Qobj:
    return (-1j * theta / 2.0 * (math.cos(phase) * qt.sigmax() + math.sin(phase) * qt.sigmay())).expm()


def two_qubit_xy_unitary(theta: float) -> np.ndarray:
    X = qt.sigmax()
    Y = qt.sigmay()
    H = -1j * theta / 2.0 * (qt.tensor(X, X) + qt.tensor(Y, Y))
    return H.expm().full()


def phase_unitary(phase_map: Dict[int, float], num_qubits: int) -> qt.Qobj:
    ops = []
    for q in range(num_qubits):
        angle = phase_map.get(q, 0.0)
        if abs(angle) <= 1e-12:
            ops.append(qt.qeye(2))
        else:
            ops.append(
                qt.Qobj(
                    [
                        [np.exp(-1j * angle / 2.0), 0.0],
                        [0.0, np.exp(1j * angle / 2.0)],
                    ]
                )
            )
    return qt.tensor(ops)


def single_qubit_bloch(rho: qt.Qobj, qubit: int, num_qubits: int) -> Tuple[float, float, float]:
    """Return Bloch vector (x, y, z) for a given qubit."""
    try:
        reduced = rho.ptrace(qubit)
    except Exception:
        # Fallback if ptrace fails due to dims metadata.
        dims = [[2] * num_qubits, [2] * num_qubits]
        rho = qt.Qobj(rho.full(), dims=dims)
        reduced = rho.ptrace(qubit)
    x = float(np.real(qt.expect(qt.sigmax(), reduced)))
    y = float(np.real(qt.expect(qt.sigmay(), reduced)))
    z = float(np.real(qt.expect(qt.sigmaz(), reduced)))
    return (x, y, z)


def bloch_rotation(before: Tuple[float, float, float], after: Tuple[float, float, float]) -> Tuple[float, Tuple[float, float, float]]:
    """Return rotation angle (rad) and axis inferred from Bloch vectors."""
    b = np.array(before, dtype=float)
    a = np.array(after, dtype=float)
    norm_b = np.linalg.norm(b)
    norm_a = np.linalg.norm(a)
    if norm_b == 0.0 or norm_a == 0.0:
        return 0.0, (0.0, 0.0, 0.0)
    b_norm = b / norm_b
    a_norm = a / norm_a
    dot = float(np.clip(np.dot(b_norm, a_norm), -1.0, 1.0))
    angle = math.acos(dot)
    axis_vec = np.cross(b_norm, a_norm)
    axis_norm = np.linalg.norm(axis_vec)
    if axis_norm > 1e-9:
        axis = tuple(axis_vec / axis_norm)
    else:
        axis = (0.0, 0.0, 0.0)
    return angle, axis


def parse_qasm_operations(path: Path) -> List[Tuple[str, List[float], List[int]]]:
    operations: List[Tuple[str, List[float], List[int]]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.split("//", 1)[0].strip()
            if not line:
                continue
            if line.startswith("OPENQASM") or line.startswith("include"):
                continue
            if line.startswith("qreg") or line.startswith("creg"):
                continue
            if line.startswith("barrier") or line.startswith("measure") or line.startswith("reset"):
                continue
            if not line.endswith(";"):
                continue
            line = line[:-1].strip()
            if not line:
                continue
            name: str
            params: List[float] = []
            rest: str
            head = line.split()[0]
            if "(" in head and ")" in line:
                name = head[: head.find("(")].lower()
                param_section = line[line.find("(") + 1 : line.find(")")]
                params = [evaluate_expression(tok) for tok in param_section.split(",") if tok.strip()]
                rest = line[line.find(")") + 1 :].strip()
            else:
                parts = line.split(None, 1)
                name = parts[0].lower()
                rest = parts[1] if len(parts) > 1 else ""
            targets = [tok.strip() for tok in rest.split(",") if tok.strip()]
            qubits: List[int] = []
            for tok in targets:
                if "[" in tok and tok.endswith("]"):
                    idx = tok[tok.find("[") + 1 : -1]
                    qubits.append(int(idx))
                else:
                    raise ValueError(f"Unable to parse qubit token '{tok}' in line '{raw_line.strip()}'")
            operations.append((name, params, qubits))
    return operations


def build_ideal_state(num_qubits: int, qasm_path: Path) -> qt.Qobj:
    ops = parse_qasm_operations(qasm_path)
    state = qt.tensor([qt.basis(2, 0)] * num_qubits)
    for name, params, qubits in ops:
        if any(q >= num_qubits or q < 0 for q in qubits):
            raise ValueError(f"Gate '{name}' targets qubit(s) {qubits} outside simulated range 0..{num_qubits - 1}")
        if name == "rx":
            if len(qubits) != 1 or not params:
                raise ValueError("rx gate requires one parameter and one target qubit")
            unitary = (-1j * params[0] / 2.0 * qt.sigmax()).expm()
            state = embed_single(unitary, qubits[0], num_qubits) * state
        elif name == "x":
            if len(qubits) != 1:
                raise ValueError("x gate requires one target qubit")
            unitary = (-1j * math.pi / 2.0 * qt.sigmax()).expm()
            state = embed_single(unitary, qubits[0], num_qubits) * state
        elif name == "sx":
            if len(qubits) != 1:
                raise ValueError("sx gate requires one target qubit")
            unitary = (-1j * (0.5 * math.pi) / 2.0 * qt.sigmax()).expm()
            state = embed_single(unitary, qubits[0], num_qubits) * state
        elif name == "rz":
            if len(qubits) != 1 or not params:
                raise ValueError("rz gate requires one parameter and one target qubit")
            unitary = qt.Qobj(
                [
                    [np.exp(-1j * params[0] / 2.0), 0.0],
                    [0.0, np.exp(1j * params[0] / 2.0)],
                ],
                dims=[[2], [2]],
            )
            state = embed_single(unitary, qubits[0], num_qubits) * state
        elif name == "iswap":
            if len(qubits) != 2:
                raise ValueError("iswap gate requires two qubits")
            full = embed_two(ISWAP_MATRIX, qubits[0], qubits[1], num_qubits)
            state = full * state
        else:
            raise ValueError(f"Unsupported gate '{name}' in QASM file")
    return qt.ket2dm(state)


def build_ideal_progression(
    num_qubits: int,
    schedule: Sequence[dict],
    pulse_library: Dict[str, dict],
) -> List[Tuple[dict, qt.Qobj]]:
    phase_map: Dict[int, float] = {q: 0.0 for q in range(num_qubits)}
    state = qt.tensor([qt.basis(2, 0)] * num_qubits)
    snapshots: List[Tuple[dict, qt.Qobj]] = []

    for event in sorted(schedule, key=lambda item: float(item.get("start_time", 0.0))):
        qubits = event.get("qubits") or []
        if not qubits:
            continue
        pulse = pulse_library.get(event.get("pulse_id"))
        if pulse is None:
            continue
        waveform_type = str(pulse.get("waveform_type", "")).lower()
        is_virtual = bool(pulse.get("virtual")) or bool(event.get("virtual")) or waveform_type == "virtual"
        theta_param = None
        params = event.get("parameters")
        if isinstance(params, dict):
            theta_param = params.get("theta")
        if theta_param is None:
            theta_param = (pulse.get("parameters") or {}).get("theta")
        event_debug = dict(event)
        duration = float(event.get("duration", 0.0))
        if is_virtual:
            theta_val = float(theta_param) if theta_param is not None else 0.0
            event_debug["ideal_theta_actual"] = theta_val
            event_debug["ideal_gate_model"] = "virtual_z"
            for q in qubits:
                phase_map[int(q)] = phase_map.get(int(q), 0.0) + theta_val
            event_debug["ideal_phase_snapshot"] = {int(q): phase_map.get(int(q), 0.0) for q in qubits}
            snapshots.append((event_debug, qt.ket2dm(state)))
            continue

        if len(qubits) == 1:
            q = int(qubits[0])
            calib = pulse.get("calibration", {})
            theta_target = float(theta_param) if theta_param is not None else float(
                pulse.get("parameters", {}).get("theta", 0.0)
            )
            if abs(theta_target) <= 1e-12 and duration > 0.0:
                theta_area = pulse_rotation_area(pulse, duration)
                if abs(theta_area) > 1e-12:
                    theta_target = theta_area
            theta_actual = float(calib.get("theta_actual_rad", theta_target))
            frame_phase = phase_map.get(q, 0.0)
            unitary = single_xy_unitary(theta_target, frame_phase)
            state = embed_single(unitary, q, num_qubits) * state
            event_debug["ideal_theta_actual"] = theta_target
            event_debug["calibration_theta_actual"] = theta_actual
            event_debug["ideal_phase_snapshot"] = {int(q): phase_map.get(int(q), 0.0) for q in range(num_qubits)}
            event_debug["ideal_gate_model"] = "single_xy"
            event_debug["ideal_drive_axis"] = {
                "cos": math.cos(frame_phase),
                "sin": math.sin(frame_phase),
            }
        elif len(qubits) == 2:
            q0, q1 = (int(qubits[0]), int(qubits[1]))
            calib = pulse.get("calibration", {})
            theta_target = float(theta_param) if theta_param is not None else float(
                pulse.get("parameters", {}).get("theta", 0.25 * math.pi)
            )
            if (theta_param is None and not calib) or abs(theta_target) <= 1e-12:
                theta_area = pulse_rotation_area(pulse, duration)
                if abs(theta_area) > 1e-12:
                    theta_target = theta_area
            theta_actual = float(calib.get("theta_actual_rad", theta_target))
            state = embed_two(two_qubit_xy_unitary(theta_target), q0, q1, num_qubits) * state
            event_debug["ideal_theta_actual"] = theta_target
            event_debug["calibration_theta_actual"] = theta_actual
            event_debug["ideal_phase_snapshot"] = {int(q): phase_map.get(int(q), 0.0) for q in range(num_qubits)}
            event_debug["ideal_gate_model"] = "two_xy"
        else:
            raise ValueError(f"Unsupported pulse targeting {len(qubits)} qubits")

        snapshots.append((event_debug, qt.ket2dm(state)))

    return snapshots


# --------------------------------------------------------------------------- #
# Pulse schedule simulation
# --------------------------------------------------------------------------- #


def simulate_pulse_schedule(
    pulse_doc: dict,
    num_qubits: int,
    edges: Iterable[Tuple[int, int]],
    rho_initial: qt.Qobj | None = None,
    t1_times: Dict[int, float] | None = None,
    t2_times: Dict[int, float] | None = None,
    relax_to_ground: bool = True,
) -> Tuple[
    qt.Qobj,
    List[qt.Qobj],
    np.ndarray,
    float,
    List[dict],
    List[dict],
    Dict[int, float],
    Dict[int, float],
    List[dict],
]:
    pulse_library = {entry["id"]: entry for entry in pulse_doc.get("pulse_library", [])}
    schedule = sorted(pulse_doc.get("schedule", []), key=lambda item: float(item.get("start_time", 0.0)))
    phase_map: Dict[int, float] = {q: 0.0 for q in range(num_qubits)}
    residual_phase: Dict[int, float] = {q: 0.0 for q in range(num_qubits)}

    edges_set = {tuple(sorted(edge)) for edge in edges}

    if not schedule:
        rho0 = rho_initial if rho_initial is not None else qt.ket2dm(qt.tensor([qt.basis(2, 0)] * num_qubits))
        return rho0, [], np.array([], dtype=float), 0.25e-9, [], {q: 0.0 for q in range(num_qubits)}

    event_info: List[dict] = []
    physical_events: List[dict] = []
    event_records: List[dict] = []
    virtual_events: List[dict] = []
    dt_base: float | None = None

    for event in schedule:
        qubits = event.get("qubits") or []
        if not qubits:
            continue
        start = float(event.get("start_time", 0.0))
        duration = float(event.get("duration", 0.0))
        pulse_id = event.get("pulse_id")
        if not pulse_id:
            continue
        if any(int(q) >= num_qubits for q in qubits):
            raise ValueError(
                f"Pulse '{pulse_id}' targets qubit(s) {qubits} outside simulated range 0..{num_qubits - 1}"
            )
        pulse = pulse_library.get(pulse_id)
        if pulse is None:
            raise ValueError(f"Pulse '{pulse_id}' referenced in schedule but missing from pulse library")
        waveform_type = str(pulse.get("waveform_type", "")).lower()
        is_virtual = bool(pulse.get("virtual")) or bool(event.get("virtual")) or waveform_type == "virtual"
        params = pulse.get("parameters") or {}

        samples_len = max(len(pulse.get("samples_i", []) or []), len(pulse.get("samples_q", []) or []))
        dt_event = None
        if not is_virtual and duration > 0.0:
            dt_event = compute_dt_from_samples(duration, samples_len, params)
            if dt_event <= 0.0:
                dt_event = 0.25e-9
            dt_base = dt_event if dt_base is None else min(dt_base, dt_event)
        event_info.append(
            {
                "event": event,
                "pulse": pulse,
                "qubits": [int(q) for q in qubits],
                "start": start,
                "duration": duration,
                "dt_event": dt_event,
                "is_virtual": is_virtual,
            }
        )

    if dt_base is None:
        dt_base = 0.25e-9

    # Determine total number of samples required on global grid.
    max_end_step = 0
    for info in event_info:
        if info["is_virtual"] or info["duration"] <= 0.0:
            continue
        start_idx = int(round(info["start"] / dt_base))
        length = max(1, int(round(info["duration"] / dt_base)))
        max_end_step = max(max_end_step, start_idx + length)

    if max_end_step == 0:
        rho0 = qt.ket2dm(qt.tensor([qt.basis(2, 0)] * num_qubits))
        phase_snapshot = {q: phase_map.get(q, 0.0) for q in range(num_qubits)}
        residual_snapshot = {q: residual_phase.get(q, 0.0) for q in range(num_qubits)}
        return (
            rho0,
            [],
            np.array([], dtype=float),
            dt_base,
            physical_events,
            virtual_events,
            phase_snapshot,
            residual_snapshot,
            event_records,
        )

    total_steps = max_end_step
    I_env = {q: np.zeros(total_steps, dtype=float) for q in range(num_qubits)}
    Q_env = {q: np.zeros(total_steps, dtype=float) for q in range(num_qubits)}
    J_env = {edge: np.zeros(total_steps, dtype=float) for edge in edges_set}

    for info in event_info:
        qubits = info["qubits"]
        start_time = info["start"]
        duration = info["duration"]
        pulse = info["pulse"]
        phase_before_full = {int(q): phase_map.get(int(q), 0.0) for q in range(num_qubits)}
        if info["is_virtual"]:
            theta = None
            params = info["event"].get("parameters")
            if isinstance(params, dict):
                theta = params.get("theta")
            if theta is None:
                theta = (pulse.get("parameters") or {}).get("theta")
            phase_before = {int(q): phase_map.get(int(q), 0.0) for q in qubits}
            theta_val = float(theta) if theta is not None else None
            if theta_val is not None:
                for q in qubits:
                    phase_map[q] = phase_map.get(int(q), 0.0) + theta_val
                    residual_phase[q] = residual_phase.get(int(q), 0.0) + theta_val
            phase_after = {int(q): phase_map.get(int(q), 0.0) for q in qubits}
            residual_after = {int(q): residual_phase.get(int(q), 0.0) for q in qubits}
            virtual_events.append(
                {
                    "event": info["event"],
                    "qubits": qubits,
                    "start": start_time,
                    "duration": duration,
                    "theta": theta_val,
                    "phase_before": phase_before,
                    "phase_after": phase_after,
                    "residual_after": residual_after,
                    "event_index": int(info["event"].get("index", -1)),
                    "phase_snapshot": {int(q): phase_map.get(int(q), 0.0) for q in range(num_qubits)},
                }
            )
            event_records.append(
                {
                    "event": info["event"],
                    "is_virtual": True,
                    "theta": theta_val,
                    "phase_before": phase_before_full,
                    "phase_after": {int(q): phase_map.get(int(q), 0.0) for q in range(num_qubits)},
                    "phase_snapshot": {int(q): phase_map.get(int(q), 0.0) for q in range(num_qubits)},
                    "residual_snapshot": {int(q): residual_phase.get(int(q), 0.0) for q in range(num_qubits)},
                    "state_index": None,
                }
            )
            continue

        if duration <= 0.0:
            continue

        start_idx = int(round(start_time / dt_base))
        length = max(1, int(round(duration / dt_base)))
        end_idx = start_idx + length
        if end_idx > total_steps:
            extend = end_idx - total_steps
            for q in range(num_qubits):
                I_env[q] = np.concatenate([I_env[q], np.zeros(extend, dtype=float)])
                Q_env[q] = np.concatenate([Q_env[q], np.zeros(extend, dtype=float)])
            for edge in list(J_env.keys()):
                J_env[edge] = np.concatenate([J_env[edge], np.zeros(extend, dtype=float)])
            total_steps = end_idx

        samples_i = pulse.get("samples_i", []) or []
        samples_q = pulse.get("samples_q", []) or []
        amplitude = float(pulse.get("amplitude", 0.0))

        if len(qubits) == 1:
            target = qubits[0]
            if samples_i or samples_q:
                I_seg = resample_to_length(samples_i, duration, length)
                Q_seg = resample_to_length(samples_q, duration, length)
            else:
                I_seg = np.full(length, amplitude, dtype=float)
                Q_seg = np.zeros(length, dtype=float)
            I_seg, Q_seg = rotate_iq(I_seg, Q_seg, phase_map.get(target, 0.0))
            I_env[target][start_idx:end_idx] = I_seg
            Q_env[target][start_idx:end_idx] = Q_seg
        elif len(qubits) == 2:
            pair = tuple(sorted(qubits))
            if samples_i:
                J_seg = resample_to_length(samples_i, duration, length)
            else:
                J_seg = np.full(length, amplitude, dtype=float)
            if pair not in edges_set:
                edges_set.add(pair)
                J_env[pair] = np.zeros(total_steps, dtype=float)
            elif pair not in J_env:
                J_env[pair] = np.zeros(total_steps, dtype=float)
            J_env[pair][start_idx:end_idx] = J_seg
        else:
            raise ValueError(f"Unsupported pulse targeting {len(qubits)} qubits")

        info["start_idx"] = start_idx
        info["end_idx"] = end_idx
        info["phase_snapshot"] = {q: phase_map.get(q, 0.0) for q in range(num_qubits)}
        event_records.append(
            {
                "event": info["event"],
                "is_virtual": False,
                "theta": None,
                "phase_before": phase_before_full,
                "phase_after": {q: phase_map.get(q, 0.0) for q in range(num_qubits)},
                "phase_snapshot": info["phase_snapshot"],
                "residual_snapshot": {int(q): residual_phase.get(int(q), 0.0) for q in range(num_qubits)},
                "state_index": end_idx,
            }
        )
        physical_events.append(info)

    time_grid = np.arange(total_steps, dtype=float) * dt_base
    if time_grid.size == 0:
        rho0 = rho_initial if rho_initial is not None else qt.ket2dm(qt.tensor([qt.basis(2, 0)] * num_qubits))
        phase_snapshot = {q: phase_map.get(q, 0.0) for q in range(num_qubits)}
        residual_snapshot = {q: residual_phase.get(q, 0.0) for q in range(num_qubits)}
        return (
            rho0,
            [],
            np.array([], dtype=float),
            dt_base,
            physical_events,
            virtual_events,
            phase_snapshot,
            residual_snapshot,
            event_records,
        )
    t_eval = np.concatenate([time_grid, [time_grid[-1] + dt_base]])
    H_terms = build_hamiltonian(time_grid, dt_base, I_env, Q_env, J_env, num_qubits, edges_set)
    rho0 = rho_initial if rho_initial is not None else qt.ket2dm(qt.tensor([qt.basis(2, 0)] * num_qubits))
    collapse_ops: List[qt.Qobj] = []
    if t1_times or t2_times:
        collapse_ops = build_collapse_ops(
            num_qubits,
            t1_times=t1_times,
            t2_times=t2_times,
            relax_to_ground=relax_to_ground,
        )

    if not H_terms and not collapse_ops:
        phase_snapshot = {q: phase_map.get(q, 0.0) for q in range(num_qubits)}
        residual_snapshot = {q: residual_phase.get(q, 0.0) for q in range(num_qubits)}
        return (
            rho0,
            [],
            t_eval,
            dt_base,
            physical_events,
            virtual_events,
            phase_snapshot,
            residual_snapshot,
            event_records,
        )

    opts = qt.Options(method="bdf", rtol=1e-7, atol=1e-9, nsteps=300000, max_step=dt_base, progress_bar=None, store_states=True)
    result = qt.mesolve(H_terms, rho0, t_eval, c_ops=collapse_ops, e_ops=[], options=opts)
    phase_snapshot = {q: phase_map.get(q, 0.0) for q in range(num_qubits)}
    residual_snapshot = {q: residual_phase.get(q, 0.0) for q in range(num_qubits)}
    return rho0, result.states, t_eval, dt_base, physical_events, virtual_events, phase_snapshot, residual_snapshot, event_records


# --------------------------------------------------------------------------- #
# Command-line interface
# --------------------------------------------------------------------------- #


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pulse", required=True, type=Path, help="Pulse JSON produced by QASMTrans (-p output)")
    parser.add_argument("--device", required=True, type=Path, help="Calibrated device JSON used for the run")
    parser.add_argument("--qasm", required=True, type=Path, help="Original QASM circuit file")
    parser.add_argument("--output", type=Path, help="Optional JSON report path")
    parser.add_argument("-v", "--verbose", type=int, default=0, help="Verbosity level (>=1 prints per-pulse fidelities)")
    parser.add_argument("--no-relax-to-ground", action="store_true", help="Do not replenish qubits to ground state during T1 relaxation.")
    parser.add_argument("--no-t1t2", action="store_true", help="Ignore T1/T2 data from the device file (unitary-only evolution).")
    args = parser.parse_args()

    pulse_doc = load_json(args.pulse)
    canonicalize_pulse_doc(pulse_doc)
    device_doc = load_json(args.device)
    metadata = device_doc.get("metadata", {})
    num_qubits = int(device_doc.get("num_qubits") or metadata.get("num_qubits") or 0)
    if num_qubits <= 0:
        raise ValueError("Unable to determine number of qubits from device file")

    schedule = pulse_doc.get("schedule", [])
    pulse_library = {entry["id"]: entry for entry in pulse_doc.get("pulse_library", [])}
    ideal_progression_overall = build_ideal_progression(num_qubits, schedule, pulse_library)
    edges = parse_edges(metadata, schedule, num_qubits)

    pulse_checks = sanity_check_pulse_library(pulse_doc)
    if pulse_checks:
        print("Pulse sanity checks:")
        tolerance = 5e-3
        for check in pulse_checks:
            if check["status"] == "virtual":
                print(f"  {check['id']:>20s}: virtual pulse (skipped)")
                continue
            if check["status"] == "skip":
                print(f"  {check['id']:>20s}: skipped ({check.get('reason', 'unhandled')})")
                continue
            fidelity = check.get("fidelity")
            expected = check.get("expected")
            delta = None
            if fidelity is not None and expected is not None:
                delta = fidelity - expected
            print(
                f"  {check['id']:>20s}: fidelity={fidelity:.6f}"
                + (f" expected={expected:.6f} delta={delta:+.3e}" if expected is not None else " expected=unknown")
            )
            if expected is not None and delta is not None and abs(delta) > tolerance:
                print(f"    WARNING: pulse '{check['id']}' deviates from expected fidelity by {delta:+.3e}")

    if args.no_t1t2:
        t1_map, t2_map = {}, {}
    else:
        t1_map, t2_map = extract_t1_t2_times(device_doc)

    (
        rho0,
        states,
        t_eval,
        dt_base,
        physical_events,
        virtual_events,
        final_phase_map,
        residual_phase_map,
        event_records,
    ) = simulate_pulse_schedule(
        pulse_doc,
        num_qubits,
        edges,
        t1_times=t1_map,
        t2_times=t2_map,
        relax_to_ground=not args.no_relax_to_ground,
    )
    rho_sim_state = states[-1] if states else rho0
    rho_sim = qt.ket2dm(rho_sim_state) if rho_sim_state.isket else rho_sim_state
    if ideal_progression_overall:
        rho_ideal = ideal_progression_overall[-1][1]
    else:
        rho_ideal = qt.ket2dm(qt.tensor([qt.basis(2, 0)] * num_qubits))
    fidelity = float(qt.fidelity(rho_sim, rho_ideal))

    print("Pulse simulation complete.")
    print(f"  Qubits           : {num_qubits}")
    print(f"  Pulses simulated : {len(schedule)}")
    print(f"  Fidelity         : {fidelity:.6f}")

    if args.verbose and event_records:
        ideal_progression = ideal_progression_overall
        print("\nPer-event detail (debug):")

        rho_latest_physical = rho0
        rho_logical = rho0
        rho_prev_phys = rho0
        ideal_prev_dm = qt.ket2dm(qt.tensor([qt.basis(2, 0)] * num_qubits))
        for record, (ideal_meta, ideal_dm) in zip(event_records, ideal_progression):
            event = record["event"]
            gate_name = event.get("gate", "?")
            qubits = event.get("qubits", [])
            start_time = float(event.get("start_time", 0.0))
            ideal_theta = ideal_meta.get("ideal_theta_actual")
            ideal_phase_snapshot = ideal_meta.get("ideal_phase_snapshot")
            ideal_gate_model = ideal_meta.get("ideal_gate_model")

            if record["is_virtual"]:
                theta_val = float(record.get("theta") or 0.0)
                theta_str = f"{theta_val:+.6f}"
                fidelity_lab = qt.fidelity(rho_latest_physical, ideal_dm)
                fidelity_event = qt.fidelity(rho_logical, ideal_dm)
                print(
                    f"    t={start_time:>10.3e}s  gate={gate_name+'(v)'}  qubits={qubits}  "
                    f"theta={theta_str}  fidelity_logical={fidelity_event:.6f}  fidelity_lab={fidelity_lab:.6f}"
                )
                ideal_prev_dm = ideal_dm
                continue

            state_idx = record.get("state_index")
            if states and state_idx is not None and state_idx < len(states):
                state_obj = states[state_idx]
            else:
                state_obj = rho_latest_physical
            rho_phys = qt.ket2dm(state_obj) if state_obj.isket else state_obj
            rho_latest_physical = rho_phys

            rho_logical = rho_phys

            time_val = float(t_eval[state_idx]) if t_eval.size and state_idx is not None else start_time
            fidelity_logical = qt.fidelity(rho_logical, ideal_dm)
            fidelity_lab = qt.fidelity(rho_phys, ideal_dm)
            extra = []
            if ideal_theta is not None:
                extra.append(f"ideal_theta={ideal_theta:+.6f}")
            if ideal_phase_snapshot:
                # Print compact snapshot keyed by qubit index.
                phase_str = "{" + ", ".join(f"{q}:{val:+.4f}" for q, val in sorted(ideal_phase_snapshot.items())) + "}"
                extra.append(f"ideal_phase={phase_str}")
            if ideal_gate_model:
                extra.append(f"model={ideal_gate_model}")
            if len(qubits) == 1:
                q = int(qubits[0])
                bloch_phys = single_qubit_bloch(rho_phys, q, num_qubits)
                bloch_phys_before = single_qubit_bloch(rho_prev_phys, q, num_qubits)
                bloch_ideal = single_qubit_bloch(ideal_dm, q, num_qubits)
                bloch_ideal_before = single_qubit_bloch(ideal_prev_dm, q, num_qubits)
                rot_angle_phys, rot_axis_phys = bloch_rotation(bloch_phys_before, bloch_phys)
                rot_angle_ideal, rot_axis_ideal = bloch_rotation(bloch_ideal_before, bloch_ideal)
                extra.append(
                    "bloch_phys=({:+.3f},{:+.3f},{:+.3f})".format(*bloch_phys)
                )
                extra.append(
                    "bloch_ideal=({:+.3f},{:+.3f},{:+.3f})".format(*bloch_ideal)
                )
                extra.append(
                    "bloch_phys_before=({:+.3f},{:+.3f},{:+.3f})".format(*bloch_phys_before)
                )
                extra.append(
                    "bloch_ideal_before=({:+.3f},{:+.3f},{:+.3f})".format(*bloch_ideal_before)
                )
                extra.append(
                    "rot_phys_angle={:+.3f}".format(rot_angle_phys)
                )
                extra.append(
                    "rot_ideal_angle={:+.3f}".format(rot_angle_ideal)
                )
                if any(abs(val) > 1e-6 for val in rot_axis_phys):
                    extra.append(
                        "rot_phys_axis=({:+.3f},{:+.3f},{:+.3f})".format(*rot_axis_phys)
                    )
                if any(abs(val) > 1e-6 for val in rot_axis_ideal):
                    extra.append(
                        "rot_ideal_axis=({:+.3f},{:+.3f},{:+.3f})".format(*rot_axis_ideal)
                    )
            extrastr = ("  " + "  ".join(extra)) if extra else ""
            print(
                f"    t={time_val:>10.3e}s  gate={gate_name:<7} qubits={qubits}  "
                f"fidelity_logical={fidelity_logical:.6f}  fidelity_lab={fidelity_lab:.6f}{extrastr}"
            )
            rho_prev_phys = rho_phys
            ideal_prev_dm = ideal_dm

        rho_readout = rho_sim
        fidelity_readout = qt.fidelity(rho_readout, rho_ideal)
        print(
            f"    t={float(t_eval[-1]) if t_eval.size else 0.0:>10.3e}s  gate=readout  qubits=all  fidelity={fidelity_readout:.6f}"
        )

        final_fidelity_debug = fidelity_readout
        print(f"\nFinal fidelity (after residual phases): {final_fidelity_debug:.6f}")


    if args.output:
        report = {
            "pulse_file": str(args.pulse),
            "device_file": str(args.device),
            "qasm_file": str(args.qasm),
            "num_qubits": num_qubits,
            "num_pulses": len(schedule),
            "fidelity": fidelity,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2)
        print(f"Wrote report to {args.output}")


if __name__ == "__main__":
    main()
