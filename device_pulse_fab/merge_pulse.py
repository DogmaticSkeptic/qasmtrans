#!/usr/bin/env python3
import argparse
import copy
import json
import math
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pulse_library", type=Path, required=True)
    parser.add_argument(
        "--compression_ratio",
        type=float,
        default=None,
        help="Fraction of the native gate duration to target (defaults to 0.5 if unspecified)",
    )
    parser.add_argument(
        "--compression_shift_ns",
        type=float,
        default=0.0,
        help="Additional time shift (ns) applied after scaling by the compression ratio",
    )
    parser.add_argument(
        "--evo_time_ns",
        type=float,
        default=None,
        help="Explicit merged gate duration in ns (overrides compression parameters)",
    )
    parser.add_argument("--coeffs", type=int, default=6)
    parser.add_argument("--amp_bound", type=float, default=None, help="Optional absolute amplitude bound applied to all controls (auto if omitted)")
    parser.add_argument("--seed_from_library", action="store_true")
    parser.add_argument("--out_json", type=Path)
    args = parser.parse_args()

    pulse_doc = load_json(args.pulse_library)
    pulse_lib = canonicalize_pulse_library(pulse_doc)
    seq = pick_sequence_from_library(pulse_lib, max_qubits=2)
    unique_qubits = []
    for spec in seq:
        for q in spec.qubits:
            if q not in unique_qubits:
                unique_qubits.append(q)
    nq = len(unique_qubits)
    if nq != 2:
        raise ValueError(f"expected two distinct qubits in sequence, found {nq}")
    logical_to_physical = {idx: q for idx, q in enumerate(unique_qubits)}
    U_targ = ideal_from_sequence(seq, nq)
    U_native = simulate_native_unitary(seq, nq)
    native_fidelity = unitary_fidelity(U_targ, U_native)
    print(f"Native sequence fidelity: {native_fidelity:.6f}")
    if native_fidelity < 0.995:
        raise ValueError(
            f"Native pulse sequence fidelity below threshold: {native_fidelity:.6f}"
        )

    total_width = float(sum(spec.width for spec in seq))
    dt_library = infer_dt_from_library(seq)
    native_time = max(total_width, dt_library)

    if args.evo_time_ns is not None and args.evo_time_ns > 0.0:
        requested_time = float(args.evo_time_ns) * 1e-9
    else:
        ratio = 0.5 if args.compression_ratio is None else float(args.compression_ratio)
        if ratio <= 0.0:
            raise ValueError("compression_ratio must be positive")
        requested_time = ratio * native_time
        requested_time += float(args.compression_shift_ns) * 1e-9

    requested_time = max(dt_library, requested_time)
    n_ts = max(1, int(math.ceil(requested_time / dt_library)))
    evo_time = n_ts * dt_library
    dt_effective = dt_library
    requested_ratio = requested_time / native_time
    actual_ratio = evo_time / native_time

    dims = [[2] * nq, [2] * nq]
    H_d = qt.Qobj(np.zeros((2 ** nq, 2 ** nq), dtype=complex), dims=dims)
    basis = build_controls(nq, logical_to_physical=logical_to_physical)
    n_ctrls = len(basis.ctrls)

    library_guess = make_library_guess(seq, basis, n_ts, evo_time)
    guess = library_guess.copy() if args.seed_from_library else None

    default_amp = 5.0
    amp_lbounds: List[float] = []
    amp_ubounds: List[float] = []
    user_amp = None if args.amp_bound is None else abs(args.amp_bound)
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

    print(
        "Optimizing merged pulse: "
        f"native_time={native_time * 1e9:.3f} ns, "
        f"requested_time={requested_time * 1e9:.3f} ns (ratio={requested_ratio:.3f}), "
        f"actual_time={evo_time * 1e9:.3f} ns (ratio={actual_ratio:.3f}), "
        f"dt={dt_library * 1e9:.3f} ns, n_ts={n_ts}"
    )

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
        max_iter=400,
        max_wall_time=120,
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

    guess_snapshot = None
    if guess is not None:
        guess = np.asarray(guess, dtype=float)
        guess_snapshot = guess.copy()
        init_amps = 0.5 * init_amps + 0.5 * guess

    lower_bounds = np.repeat(np.asarray(amp_lbounds_scaled)[np.newaxis, :], n_ts, axis=0)
    upper_bounds = np.repeat(np.asarray(amp_ubounds_scaled)[np.newaxis, :], n_ts, axis=0)
    init_amps = np.clip(init_amps, lower_bounds, upper_bounds)

    dyn.initialize_controls(init_amps.copy())

    res = optim.run_optimization()
    optimized_fidelity = max(0.0, 1.0 - float(np.real(res.fid_err)))
    print(f"Optimized pulse fidelity: {optimized_fidelity:.6f}")

    final_amps = np.array(res.final_amps)
    final_amps_physical = final_amps * ctrl_scale[np.newaxis, :]
    final_amps_real = np.real(final_amps_physical)
    final_amps_imag = np.imag(final_amps_physical)
    U_mat = U_targ.full()
    U_real = np.real(U_mat)
    U_imag = np.imag(U_mat)

    plot_path = None
    try:
        library_stem = args.pulse_library.stem if args.pulse_library else "pulse"
        plot_filename = f"{library_stem}_merged_pulse.png"
        plot_path_candidate = Path(__file__).resolve().parent / plot_filename
        plot_saved = save_pulse_plot(
            plot_path_candidate,
            final_amps_real,
            basis.labels,
            dt_effective,
        )
        if plot_saved is not None:
            plot_path = plot_saved
    except Exception as exc:
        print(f"Warning: failed to save pulse plot ({exc})", file=sys.stderr)

    out = {
        "sequence": [
            {
                "gate": spec.gate,
                "theta": float(spec.theta),
                "qubits": list(spec.qubits),
                "local_qubits": list(local_qubits(spec)),
            }
            for spec in seq
        ],
        "physical_qubits": unique_qubits,
        "logical_to_physical": {int(k): int(v) for k, v in logical_to_physical.items()},
        "target_unitary": {
            "real": U_real.tolist(),
            "imag": U_imag.tolist(),
        },
        "evo_time_s": float(evo_time),
        "dt_library_s": float(dt_library),
        "dt_effective_s": float(dt_effective),
        "total_library_time_s": float(total_width),
        "num_tslots": int(n_ts),
        "final_fid_err": float(np.real(res.fid_err)),
        "final_amps": {
            "real": final_amps_real.tolist(),
            "imag": final_amps_imag.tolist(),
        },
        "ctrls_order": basis.labels,
        "control_scaling": ctrl_scale.tolist(),
        "native_fidelity": float(native_fidelity),
        "optimized_fidelity": float(optimized_fidelity),
        "amp_lbounds": amp_lbounds,
        "amp_ubounds": amp_ubounds,
        "native_total_time_s": float(native_time),
        "requested_evo_time_s": float(requested_time),
        "compression_ratio_requested": float(requested_ratio),
        "compression_ratio_actual": float(actual_ratio),
    }
    if plot_path is not None:
        out["pulse_plot"] = str(plot_path)
    if guess_snapshot is not None:
        guess_phys = guess_snapshot * ctrl_scale[np.newaxis, :]
        guess_real = np.real(guess_phys)
        guess_imag = np.imag(guess_phys)
        out["initial_guess"] = {
            "real": guess_real.tolist(),
            "imag": guess_imag.tolist(),
        }

    if args.out_json:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        with args.out_json.open("w", encoding="utf-8") as h:
            json.dump(out, h, indent=2)
    else:
        print(f"Native fidelity: {native_fidelity:.6f}")
        print(f"Optimized fidelity: {optimized_fidelity:.6f}")
        print(f"Final fidelity error: {float(np.real(res.fid_err)):.6e}")
        print(f"Compression ratio (requested/actual): {requested_ratio:.3f}/{actual_ratio:.3f}")
        if plot_path is not None:
            print(f"Pulse plot saved to: {plot_path}")


if __name__ == "__main__":
    main()
