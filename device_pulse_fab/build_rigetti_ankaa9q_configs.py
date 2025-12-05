#!/usr/bin/env python3
"""Build Rigetti Ankaa-9Q configs and pulses from devicelib metadata."""
from __future__ import annotations

import argparse
import copy
import json
import math
import re
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple

import numpy as np
import qutip as qt

ROOT = Path(__file__).resolve().parents[1]
DEVICELIB_PATH = ROOT / "devicelib" / "rigetti_Ankaa_9Q_3_n9.json"
OUT_DEVICE_PATH = ROOT / "data" / "devices" / "rigetti_ankaa9q3_device.json"
OUT_PULSE_PATH = ROOT / "data" / "devices" / "rigetti_ankaa9q3_pulses.json"
DT_NS = 0.25
SIGMA_FRAC = 0.2
FIDELITY_TOL = 5e-3
RX_ANGLES = [
    ("pi", math.pi),
    ("neg_pi", -math.pi),
    ("pi_over_2", 0.5 * math.pi),
    ("neg_pi_over_2", -0.5 * math.pi),
]
EPS = 1e-12


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def dump_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def sanitize_metadata(value):
    if isinstance(value, dict):
        return {k: sanitize_metadata(v) for k, v in value.items()}
    if isinstance(value, list):
        return [sanitize_metadata(v) for v in value]
    if value is None:
        return 0
    return value


INDEX_PATTERN = re.compile(r"\d+")


def remap_key_with_indices(key: str, mapping: Dict[int, int]) -> str | None:
    result = []
    last = 0
    for match in INDEX_PATTERN.finditer(key):
        idx = int(match.group())
        if idx not in mapping:
            return None
        result.append(key[last:match.start()])
        result.append(str(mapping[idx]))
        last = match.end()
    result.append(key[last:])
    return "".join(result)


def filter_devicelib(devicelib_doc: dict, selected_qubits: List[int], mapping: Dict[int, int]) -> dict:
    selected_set = set(selected_qubits)
    filtered = copy.deepcopy(devicelib_doc)

    filtered["num_qubits"] = len(selected_qubits)
    filtered["coupling"] = []
    for edge in devicelib_doc.get("coupling", []):
        parts = str(edge).split("_")
        try:
            mapped = [mapping[int(p)] for p in parts if int(p) in selected_set]
        except KeyError:
            continue
        if len(mapped) == len(parts):
            filtered["coupling"].append("_".join(str(v) for v in mapped))

    def remap_gate_map(map_obj: dict) -> dict:
        remapped = {}
        for key, value in map_obj.items():
            new_key = remap_key_with_indices(key, mapping)
            if new_key is not None:
                remapped[new_key.lower()] = value
        return remapped

    def remap_index_map(map_obj: dict) -> dict:
        return {
            str(mapping[int(key)]): value
            for key, value in map_obj.items()
            if int(key) in mapping
        }

    for field in ("gate_lens", "gate_errs"):
        if field in filtered and isinstance(filtered[field], dict):
            filtered[field] = remap_gate_map(filtered[field])

    for field in ("T1", "T2", "readout_length", "prob_meas0_prep1", "prob_meas1_prep0", "freq"):
        if field in filtered and isinstance(filtered[field], dict):
            filtered[field] = remap_index_map(filtered[field])

    return filtered


def drop_two_qubit_gate_entries(map_obj: dict) -> dict:
    """Remove two-qubit gate metadata (e.g., iswap0_1) when coupling is disabled."""
    return {k: v for k, v in map_obj.items() if "_" not in str(k)}


def build_device_document(devicelib_doc: dict) -> Tuple[dict, Dict[str, float], Dict[str, float]]:
    meta_keys = [
        "data_source",
        "last_update_date",
        "online_date",
        "anharmonicity",
        "processor",
        "name",
        "version",
        "num_qubits",
        "basis_gates",
        "T1",
        "T2",
        "freq",
        "readout_length",
        "prob_meas0_prep1",
        "prob_meas1_prep0",
        "gate_lens",
        "gate_errs",
        "coupling",
    ]
    metadata = {key: sanitize_metadata(devicelib_doc.get(key)) for key in meta_keys}
    basis_lower = [str(gate).lower() for gate in metadata.get("basis_gates", [])]
    gate_lens_raw = metadata.get("gate_lens", {})
    gate_errs_raw = metadata.get("gate_errs", {})
    gate_lens = {k.lower(): float(v) for k, v in gate_lens_raw.items()}
    gate_errs = {k.lower(): float(v) for k, v in gate_errs_raw.items()}
    device_doc = {
        "metadata": metadata,
        "name": "rigetti_ankaa9q3_calibrated",
        "version": metadata.get("last_update_date", ""),
        "num_qubits": int(metadata.get("num_qubits", 0)),
        "basis_gates": basis_lower,
        "gate_lens": gate_lens,
        "gate_errs": gate_errs,
        "cx_coupling": metadata.get("coupling", []),
        "physical_qubits": list(range(int(metadata.get("num_qubits", 0)))),
        "readout_length": metadata.get("readout_length", {}),
        "prob_meas0_prep1": metadata.get("prob_meas0_prep1", {}),
        "prob_meas1_prep0": metadata.get("prob_meas1_prep0", {}),
    }
    return device_doc, gate_lens, gate_errs


def gaussian_envelope(duration_ns: float, dt_ns: float, sigma_frac: float) -> Tuple[np.ndarray, np.ndarray]:
    dt = dt_ns * 1e-9
    duration = duration_ns * 1e-9
    samples = int(round(duration / dt)) + 1
    times = np.arange(samples, dtype=float) * dt
    sigma = sigma_frac * duration
    center = 0.5 * duration
    envelope = np.exp(-0.5 * ((times - center) / sigma) ** 2)
    return times, envelope


def rx_waveform(theta: float, duration_ns: float, dt_ns: float, sigma_frac: float) -> Tuple[np.ndarray, np.ndarray, float, float]:
    times, envelope = gaussian_envelope(duration_ns, dt_ns, sigma_frac)
    dt = dt_ns * 1e-9
    # numpy 1.24 lacks np.trapezoid; np.trapz is equivalent here
    area = float(np.trapz(envelope, dx=dt))
    if area <= EPS:
        raise ValueError("Gaussian envelope area is zero; check dt_ns and duration_ns")
    amplitude = theta / area
    samples = amplitude * envelope
    return times, samples, amplitude, sigma_frac * duration_ns


def simulate_single_qubit_rx(
    theta: float,
    duration_ns: float,
    dt_ns: float,
    sigma_frac: float,
    theta_target: float | None = None,
) -> float:
    times, samples, _, _ = rx_waveform(theta, duration_ns, dt_ns, sigma_frac)
    coeff = lambda t, _args=None: float(np.interp(t, times, samples, left=samples[0], right=samples[-1]))
    H = [0.5 * qt.sigmax(), coeff]
    rho0 = qt.ket2dm(qt.basis(2, 0))
    dt = dt_ns * 1e-9
    options = qt.Options(nsteps=300000, atol=1e-9, rtol=1e-7, max_step=dt)
    result = qt.mesolve(H, rho0, times, c_ops=[], options=options)
    rho_f = result.states[-1]
    if rho_f.isket:
        rho_f = qt.ket2dm(rho_f)
    theta_ref = theta if theta_target is None else float(theta_target)
    target_state = (-1j * theta_ref / 2.0 * qt.sigmax()).expm() * qt.basis(2, 0)
    rho_target = qt.ket2dm(target_state)
    fidelity = float(qt.fidelity(rho_f, rho_target))
    return fidelity


def simulate_iswap(
    theta: float,
    duration_ns: float,
    dt_ns: float,
    theta_target: float | None = None,
) -> float:
    duration = duration_ns * 1e-9
    dt = dt_ns * 1e-9
    samples = int(round(duration / dt)) + 1
    times = np.arange(samples, dtype=float) * dt
    amplitude = theta / duration if duration > EPS else 0.0
    coeff = lambda _t, _args=None: amplitude
    sx = qt.sigmax()
    sy = qt.sigmay()
    H_xy = 0.5 * (qt.tensor(sx, sx) + qt.tensor(sy, sy))
    H = [H_xy, coeff]
    rho0 = qt.ket2dm(qt.tensor(qt.basis(2, 1), qt.basis(2, 0)))
    options = qt.Options(nsteps=300000, atol=1e-9, rtol=1e-7, max_step=dt)
    result = qt.mesolve(H, rho0, times, c_ops=[], options=options)
    rho_f = result.states[-1]
    if rho_f.isket:
        rho_f = qt.ket2dm(rho_f)
    theta_ref = theta if theta_target is None else float(theta_target)
    U_target = (-1j * theta_ref * H_xy).expm()
    rho_target = U_target * rho0 * U_target.dag()
    fidelity = float(qt.fidelity(rho_f, rho_target))
    return fidelity


def tune_scale_to_fidelity(
    evaluator: Callable[[float], float],
    target_fidelity: float,
    tolerance: float,
    max_scale: float = 5.0,
) -> Tuple[float, float]:
    cache: Dict[float, float] = {}

    def eval_cached(scale: float) -> float:
        if scale not in cache:
            cache[scale] = evaluator(scale)
        return cache[scale]

    best_scale: float | None = None
    best_fidelity: float | None = None

    def update_best(scale: float) -> None:
        nonlocal best_scale, best_fidelity
        fidelity = eval_cached(scale)
        if best_scale is None or abs(fidelity - target_fidelity) < abs(best_fidelity - target_fidelity):
            best_scale = scale
            best_fidelity = fidelity

    update_best(0.0)
    f_low = eval_cached(0.0) - target_fidelity
    update_best(1.0)
    f_high = eval_cached(1.0) - target_fidelity
    if abs(f_high) <= tolerance:
        return 1.0, eval_cached(1.0)

    low = 0.0
    high = 1.0
    while f_low * f_high > 0.0 and high <= max_scale:
        high *= 1.5
        update_best(high)
        f_high = eval_cached(high) - target_fidelity

    if f_low * f_high > 0.0:
        if best_scale is None or best_fidelity is None:
            raise ValueError("Unable to bracket fidelity target")
        return best_scale, best_fidelity

    for _ in range(80):
        mid = 0.5 * (low + high)
        update_best(mid)
        f_mid = eval_cached(mid) - target_fidelity
        if abs(f_mid) <= tolerance:
            return mid, eval_cached(mid)
        if f_low * f_mid <= 0.0:
            high = mid
            f_high = f_mid
        else:
            low = mid
            f_low = f_mid

    if best_scale is None or best_fidelity is None:
        raise ValueError("Failed to calibrate fidelity")
    return best_scale, best_fidelity


def calibrate_single_qubit_theta(
    theta_nominal: float,
    target_fidelity: float,
    duration_ns: float,
    dt_ns: float,
    sigma_frac: float,
    tolerance: float,
) -> Tuple[float, float]:
    target = max(0.0, min(1.0, float(target_fidelity)))
    if target >= 1.0 - 1e-9:
        theta_eff = theta_nominal
        fidelity = simulate_single_qubit_rx(theta_eff, duration_ns, dt_ns, sigma_frac, theta_target=theta_nominal)
        return theta_eff, fidelity

    eff_tolerance = min(tolerance, max(1e-6, 0.05 * (1.0 - target)))

    def evaluator(scale: float) -> float:
        theta_actual = theta_nominal * scale
        return simulate_single_qubit_rx(theta_actual, duration_ns, dt_ns, sigma_frac, theta_target=theta_nominal)

    scale, fidelity = tune_scale_to_fidelity(evaluator, target, eff_tolerance)
    return theta_nominal * scale, fidelity


def calibrate_iswap_theta(
    theta_nominal: float,
    target_fidelity: float,
    duration_ns: float,
    dt_ns: float,
    tolerance: float,
) -> Tuple[float, float]:
    target = max(0.0, min(1.0, float(target_fidelity)))
    if target >= 1.0 - 1e-9:
        theta_eff = theta_nominal
        fidelity = simulate_iswap(theta_eff, duration_ns, dt_ns, theta_target=theta_nominal)
        return theta_eff, fidelity

    eff_tolerance = min(tolerance, max(1e-6, 0.05 * (1.0 - target)))

    def evaluator(scale: float) -> float:
        theta_actual = theta_nominal * scale
        return simulate_iswap(theta_actual, duration_ns, dt_ns, theta_target=theta_nominal)

    scale, fidelity = tune_scale_to_fidelity(evaluator, target, eff_tolerance)
    return theta_nominal * scale, fidelity


def build_pulse_library(
    metadata: dict,
    gate_lens: Dict[str, float],
    gate_errs: Dict[str, float],
    dt_ns: float,
    sigma_frac: float,
    tolerance: float,
    ideal: bool = False,
) -> Tuple[List[dict], Dict[str, float]]:
    num_qubits = int(metadata.get("num_qubits", 0))
    pulses: List[dict] = []
    fidelity_map: Dict[str, float] = {}
    for q in range(num_qubits):
        duration_s = float(gate_lens.get(f"rx{q}", 0.0))
        duration_ns = duration_s * 1e9
        gate_err = float(gate_errs.get(f"rx{q}", 0.0))
        target_fidelity = 1.0 if ideal else max(0.0, 1.0 - gate_err)
        for label, theta in RX_ANGLES:
            if ideal:
                theta_eff = theta
                fidelity = 1.0
            else:
                theta_eff, fidelity = calibrate_single_qubit_theta(
                    theta,
                    target_fidelity,
                    duration_ns,
                    dt_ns,
                    sigma_frac,
                    tolerance,
                )
            _, samples, amplitude, sigma_ns = rx_waveform(theta_eff, duration_ns, dt_ns, sigma_frac)
            fidelity_map[f"rx_q{q}_{label}"] = fidelity
            calibration = {
                "theta_actual_rad": float(theta_eff),
                "duration_ns": float(duration_ns),
                "sigma_ns": float(sigma_ns),
                "omega_peak_rad_per_s": float(amplitude),
                "default_phase_rad": 0.0,
                "expected_fidelity": fidelity,
            }
            pulses.append(
                {
                    "id": f"rx_q{q}_{label}",
                    "gate": "rx",
                    "qubits": [q],
                    "shape": "arbitrary",
                    "waveform_type": "arbitrary",
                    "width": duration_s,
                    "amplitude": float(np.max(np.abs(samples))),
                    "samples_i": samples.tolist(),
                    "samples_q": [0.0] * len(samples),
                    "parameters": {"theta": float(theta)},
                    "calibration": calibration,
                    "note": f"rx pulse on qubit {q}; Gaussian envelope sigma_frac={sigma_frac}",
                }
            )
        rz_duration_s = float(gate_lens.get(f"rz{q}", 0.0))
        pulses.append(
            {
                "id": f"rz_q{q}",
                "gate": "rz",
                "qubits": [q],
                    "shape": "virtual",
                    "waveform_type": "virtual",
                    "width": rz_duration_s,
                    "amplitude": 0.0,
                    "virtual": True,
                    "calibration": {
                        "duration_ns": float(rz_duration_s * 1e9),
                    },
                    "note": f"virtual frame change for qubit {q}",
                }
        )
    theta_nominal = 0.25 * math.pi
    for coupling in metadata.get("coupling", []):
        left, right = coupling.split("_")
        qi = int(left)
        qj = int(right)
        key = f"iswap{qi}_{qj}"
        duration_s = float(gate_lens.get(key, gate_lens.get(key.lower(), 0.0)))
        duration_ns = duration_s * 1e9
        gate_err = float(gate_errs.get(key, gate_errs.get(key.lower(), 0.0)))
        target_fidelity = 1.0 if ideal else max(0.0, 1.0 - gate_err)
        if ideal:
            theta_eff = theta_nominal
            fidelity = 1.0
        else:
            theta_eff, fidelity = calibrate_iswap_theta(
                theta_nominal,
                target_fidelity,
                duration_ns,
                dt_ns,
                tolerance,
            )
        amplitude = theta_eff / (duration_s if duration_s > EPS else 1.0)
        samples = np.full(int(round(duration_ns / dt_ns)) + 1, amplitude, dtype=float)
        fidelity_map[f"iswap_q{qi}_q{qj}"] = fidelity
        calibration = {
            "theta_actual_rad": float(theta_eff),
            "duration_ns": float(duration_ns),
            "J_amp_rad_per_s": float(amplitude),
            "J_area_rad": float(theta_eff),
            "expected_fidelity": fidelity,
        }
        pulses.append(
            {
                "id": f"iswap_q{qi}_q{qj}",
                "gate": "iswap",
                "qubits": [qi, qj],
                "shape": "flat_top",
                "waveform_type": "flat_top",
                "width": duration_s,
                "amplitude": float(amplitude),
                "samples_i": samples.tolist(),
                "samples_q": [0.0] * len(samples),
                "parameters": {"theta": float(theta_nominal)},
                "calibration": calibration,
                "note": f"sqrt(iSWAP) on qubits {qi}-{qj}; constant exchange {amplitude:0.6e} rad/s",
            }
        )
    return pulses, fidelity_map


def compose_pulse_document(metadata: dict, pulses: Iterable[dict], name: str) -> dict:
    return {
        "name": name,
        "version": metadata.get("last_update_date", ""),
        "num_qubits": int(metadata.get("num_qubits", 0)),
        "basis_gates": [str(gate).lower() for gate in metadata.get("basis_gates", [])],
        "pulse_definitions": list(pulses),
    }


def print_fidelity_table(pulse_entries: Iterable[dict]) -> None:
    rows: List[Tuple[str, str, float | None, float | None, float]] = []
    for entry in pulse_entries:
        params = entry.get("parameters", {})
        calib = entry.get("calibration", {})
        fidelity = calib.get("expected_fidelity")
        if fidelity is None:
            continue
        theta = params.get("theta")
        duration_ns = calib.get("duration_ns")
        qubit_label = ",".join(str(q) for q in entry.get("qubits", []))
        rows.append((entry.get("id", ""), qubit_label, theta, duration_ns, float(fidelity)))
    if not rows:
        return
    headers = ["pulse_id", "qubits", "theta_rad", "duration_ns", "expected_fidelity"]
    widths = [len(h) for h in headers]
    for pulse_id, qubits, theta, duration_ns, fidelity in rows:
        widths[0] = max(widths[0], len(pulse_id))
        widths[1] = max(widths[1], len(qubits))
        widths[2] = max(widths[2], len(f"{theta:.6f}" if theta is not None else "-"))
        widths[3] = max(widths[3], len(f"{duration_ns:.3f}" if duration_ns is not None else "-"))
        widths[4] = max(widths[4], len(f"{fidelity:.6f}"))
    fmt = "  ".join(f"{{:{w}}}" for w in widths)
    print(fmt.format(*headers))
    print("-" * (sum(widths) + 8))
    for pulse_id, qubits, theta, duration_ns, fidelity in sorted(rows):
        theta_str = f"{theta:.6f}" if theta is not None else "-"
        dur_str = f"{duration_ns:.3f}" if duration_ns is not None else "-"
        print(fmt.format(pulse_id, qubits, theta_str, dur_str, f"{fidelity:.6f}"))


def build_pulse_index(pulses: Iterable[dict]) -> Dict[str, dict]:
    return {entry["id"]: entry for entry in pulses}


def extract_waveform(entry: dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    width = float(entry.get("width", 0.0))
    samples_i = np.asarray(entry.get("samples_i", []), dtype=float)
    samples_q = np.asarray(entry.get("samples_q", []), dtype=float)
    n = max(samples_i.size, samples_q.size)
    if n < 2:
        n = 2
    if samples_i.size == 0:
        samples_i = np.zeros(n, dtype=float)
    elif samples_i.size != n:
        samples_i = np.resize(samples_i, n)
    if samples_q.size == 0:
        samples_q = np.zeros(n, dtype=float)
    elif samples_q.size != n:
        samples_q = np.resize(samples_q, n)
    times = np.linspace(0.0, width, n, dtype=float)
    dt = times[1] - times[0] if times.size > 1 else (width if width > 0.0 else DT_NS * 1e-9)
    return times, samples_i, samples_q, dt


def pulse_id_for(kind: str, q: int) -> str:
    if kind == "rx_pi":
        return f"rx_q{q}_pi"
    if kind == "rx_pi2":
        return f"rx_q{q}_pi_over_2"
    if kind == "rx_mpi2":
        return f"rx_q{q}_neg_pi_over_2"
    raise ValueError(f"Unsupported pulse kind '{kind}' for qubit {q}")


LEGACY_NQ = 2
LEGACY_EDGE = (0, 1)


def single_qubit_ops():
    e0 = qt.basis(2, 0)
    e1 = qt.basis(2, 1)
    s01 = e0 * e1.dag()
    s10 = s01.dag()
    X01 = s01 + s10
    Y01 = -1j * s01 + 1j * s10
    Z01 = e0 * e0.dag() - e1 * e1.dag()
    return e0, e1, s01, s10, X01, Y01, Z01


def lift_one(op, site, N):
    ops = [qt.qeye(2)] * N
    ops = list(ops)
    ops[site] = op
    return qt.tensor(ops)


def lift_two(op_i, i, op_j, j, N):
    ops = [qt.qeye(2)] * N
    ops = list(ops)
    ops[i] = op_i
    ops[j] = op_j
    return qt.tensor(ops)


def coeff_from_array(arr, dt):
    y = np.asarray(arr, dtype=float)
    N = y.size
    inv_dt = 1.0 / dt

    def c(t, args=None):
        idx = int(t * inv_dt)
        if idx < 0:
            idx = 0
        elif idx >= N:
            idx = N - 1
        return float(y[idx])

    return c


def build_H(t, dt, I_env, Q_env, J_env, N, edges):
    _, _, s01, s10, X, Y, _ = single_qubit_ops()
    H = []
    for q in range(N):
        Ii = I_env.get(q, np.zeros_like(t))
        Qi = Q_env.get(q, np.zeros_like(t))
        if np.any(Ii):
            H.append([lift_one(X, q, N), coeff_from_array(0.5 * Ii, dt)])
        if np.any(Qi):
            H.append([lift_one(Y, q, N), coeff_from_array(0.5 * Qi, dt)])
    for e in edges:
        i, j = e
        Jij = J_env.get(e, np.zeros_like(t))
        if np.any(Jij):
            exch = lift_two(s10, i, s01, j, N) + lift_two(s01, i, s10, j, N)
            H.append([exch, coeff_from_array(Jij, dt)])
    return H


def simulate_legacy(t, H, rho0, dt, c_ops=None):
    opts = {
        "method": "bdf",
        "rtol": 1e-7,
        "atol": 1e-9,
        "nsteps": 300000,
        "max_step": dt,
        "progress_bar": False,
        "store_states": True,
    }
    return qt.mesolve(H, rho0, t, c_ops=c_ops or [], e_ops=[], options=opts)


def make_time_grid(dt, duration):
    n = int(np.round(duration / dt)) + 1
    return np.arange(n, dtype=float) * dt


def gaussian_centered(t, duration, sigma):
    t0 = 0.5 * duration
    return np.exp(-0.5 * ((t - t0) / sigma) ** 2)


def omega_gaussian_for_rotation(theta, duration, dt, sigma_frac):
    t = make_time_grid(dt, duration)
    sigma = sigma_frac * duration
    g = gaussian_centered(t, duration, sigma)
    area = float(np.trapz(g, dx=dt))
    A = theta / area
    return t, A * g, A, sigma


def iq_from_omega_phase(omega_t, phase):
    I = omega_t * np.cos(phase)
    Q = omega_t * np.sin(phase)
    return I, Q


def run_rx(
    theta,
    duration_ns,
    target,
    dt_ns=DT_NS,
    sigma_frac=SIGMA_FRAC,
    phase=0.0,
    pulse_entry: dict | None = None,
):
    dt = dt_ns * 1e-9
    duration = duration_ns * 1e-9
    if pulse_entry is None:
        t, omega, A_peak, sigma = omega_gaussian_for_rotation(theta, duration, dt, sigma_frac)
        Iq = {q: np.zeros_like(t) for q in range(LEGACY_NQ)}
        Qq = {q: np.zeros_like(t) for q in range(LEGACY_NQ)}
        I, Q = iq_from_omega_phase(omega, phase)
        Iq[target] = I
        Qq[target] = Q
        J = {LEGACY_EDGE: np.zeros_like(t)}
        params = {
            "type": "gaussian",
            "target_qubit": int(target),
            "theta_rad": float(theta),
            "theta_actual_rad": float(theta),
            "duration_ns": float(duration_ns),
            "dt_ns": float(dt_ns),
            "sigma_ns": float(sigma * 1e9),
            "omega_peak_rad_per_s": float(A_peak),
            "I_peak_rad_per_s": float(A_peak),
            "default_phase_rad": float(phase),
            "expected_fidelity": None,
        }
    else:
        params_entry = dict(pulse_entry.get("parameters", {}))
        calib_entry = dict(pulse_entry.get("calibration", {}))
        theta = float(params_entry.get("theta", theta))
        theta_actual = float(calib_entry.get("theta_actual_rad", theta))
        duration = float(calib_entry.get("duration_ns", params_entry.get("duration_ns", duration_ns))) * 1e-9
        duration_ns = duration * 1e9
        t, samples_i, samples_q, dt = extract_waveform(pulse_entry)
        if abs(phase) > EPS:
            cos_p = math.cos(phase)
            sin_p = math.sin(phase)
            samples_i, samples_q = cos_p * samples_i - sin_p * samples_q, sin_p * samples_i + cos_p * samples_q
        Iq = {q: np.zeros_like(samples_i) for q in range(LEGACY_NQ)}
        Qq = {q: np.zeros_like(samples_q) for q in range(LEGACY_NQ)}
        Iq[target] = samples_i
        Qq[target] = samples_q
        J = {LEGACY_EDGE: np.zeros_like(samples_i)}
        params = {
            "type": pulse_entry.get("waveform_type", "arbitrary"),
            "target_qubit": int(target),
            "theta_rad": float(theta),
            "theta_actual_rad": float(theta_actual),
            "duration_ns": float(duration_ns),
            "dt_ns": float(dt * 1e9),
            "sigma_ns": float(calib_entry.get("sigma_ns", SIGMA_FRAC * duration_ns)),
            "omega_peak_rad_per_s": float(calib_entry.get("omega_peak_rad_per_s", np.max(np.abs(samples_i)))),
            "I_peak_rad_per_s": float(calib_entry.get("omega_peak_rad_per_s", np.max(np.abs(samples_i)))),
            "default_phase_rad": float(calib_entry.get("default_phase_rad", phase)),
            "expected_fidelity": calib_entry.get("expected_fidelity"),
        }
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    H = build_H(t, dt, Iq, Qq, J, LEGACY_NQ, [LEGACY_EDGE])
    res = simulate_legacy(t, H, rho0, dt)
    return t, Iq, Qq, J, res, params


def run_iswap(duration_ns, dt_ns=DT_NS, pulse_entry: dict | None = None):
    if pulse_entry is None:
        dt = dt_ns * 1e-9
        duration = duration_ns * 1e-9
        t = make_time_grid(dt, duration)
        J_area = 0.25 * np.pi
        J_amp = J_area / duration
        J = {LEGACY_EDGE: np.full_like(t, J_amp)}
        params = {
            "type": "square_iswap",
            "edge": [int(LEGACY_EDGE[0]), int(LEGACY_EDGE[1])],
            "duration_ns": float(duration_ns),
            "dt_ns": float(dt_ns),
            "J_amp_rad_per_s": float(J_amp),
            "J_area_rad": float(J_area),
            "theta_actual_rad": float(J_area),
        }
    else:
        params_entry = dict(pulse_entry.get("parameters", {}))
        calib_entry = dict(pulse_entry.get("calibration", {}))
        duration = float(calib_entry.get("duration_ns", params_entry.get("duration_ns", duration_ns))) * 1e-9
        duration_ns = duration * 1e9
        t, samples_i, samples_q, dt = extract_waveform(pulse_entry)
        J = {LEGACY_EDGE: samples_i}
        theta_actual = float(calib_entry.get("theta_actual_rad", np.trapz(samples_i, dx=dt)))
        params = {
            "type": pulse_entry.get("waveform_type", "flat_top"),
            "edge": [int(LEGACY_EDGE[0]), int(LEGACY_EDGE[1])],
            "duration_ns": float(duration_ns),
            "dt_ns": float(dt * 1e9),
            "J_amp_rad_per_s": float(calib_entry.get("J_amp_rad_per_s", np.max(np.abs(samples_i)))),
            "J_area_rad": float(calib_entry.get("J_area_rad", np.trapz(samples_i, dx=dt))),
            "theta": float(calib_entry.get("theta_actual_rad", theta_actual)),
            "theta_actual_rad": float(theta_actual),
            "expected_fidelity": calib_entry.get("expected_fidelity"),
        }
    dt = dt if pulse_entry is not None else dt_ns * 1e-9
    t = make_time_grid(dt, duration) if pulse_entry is None else t
    Iq = {0: np.zeros_like(t), 1: np.zeros_like(t)}
    Qq = {0: np.zeros_like(t), 1: np.zeros_like(t)}
    H = build_H(t, dt, Iq, Qq, J, LEGACY_NQ, [LEGACY_EDGE])
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e1, e0]))
    res = simulate_legacy(t, H, rho0, dt)
    return t, Iq, Qq, J, res, params


def run_xy_entangler(theta, duration_ns, dt_ns=DT_NS):
    dt = dt_ns * 1e-9
    duration = duration_ns * 1e-9
    t = make_time_grid(dt, duration)
    J_area = float(theta)
    J_amp = J_area / duration
    J = {LEGACY_EDGE: np.full_like(t, J_amp)}
    Iq = {0: np.zeros_like(t), 1: np.zeros_like(t)}
    Qq = {0: np.zeros_like(t), 1: np.zeros_like(t)}
    H = build_H(t, dt, Iq, Qq, J, LEGACY_NQ, [LEGACY_EDGE])
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    res = simulate_legacy(t, H, rho0, dt)
    params = {
        "type": "square_xy",
        "edge": [int(LEGACY_EDGE[0]), int(LEGACY_EDGE[1])],
        "duration_ns": float(duration_ns),
        "dt_ns": float(dt_ns),
        "J_amp_rad_per_s": float(J_amp),
        "J_area_rad": float(J_area),
    }
    return t, Iq, Qq, J, res, params


def rz_virtual_update(phase_map, q, zeta):
    phase_map[q] = phase_map.get(q, 0.0) + float(zeta)


def rz_unitary(zeta):
    z = float(zeta)
    return qt.Qobj(np.array([[np.exp(-1j * z / 2.0), 0.0], [0.0, np.exp(1j * z / 2.0)]]))


def ry_unitary(theta):
    _, _, _, _, _, Y, _ = single_qubit_ops()
    return (-1j * theta / 2.0 * Y).expm()


def su2_from_zyz(alpha, beta, gamma):
    return rz_unitary(alpha) * ry_unitary(beta) * rz_unitary(gamma)


def xy_unitary(theta, phi):
    _, _, _, _, X, Y, _ = single_qubit_ops()
    return (-1j * theta / 2.0 * (np.cos(phi) * X + np.sin(phi) * Y)).expm()


def fidelity_dm(rho, sigma):
    if rho.isket:
        rho = qt.ket2dm(rho)
    if sigma.isket:
        sigma = qt.ket2dm(sigma)
    return float(qt.fidelity(rho, sigma))


def simulate_sequence_on_q0(blocks, pulses_map, dt_ns=DT_NS, sigma_frac=SIGMA_FRAC):
    dt = dt_ns * 1e-9
    phase = {0: 0.0, 1: 0.0}
    t_all = np.array([], dtype=float)
    I_all = {0: np.zeros(0), 1: np.zeros(0)}
    Q_all = {0: np.zeros(0), 1: np.zeros(0)}
    J_all = {LEGACY_EDGE: np.zeros(0)}
    for kind, args in blocks:
        if kind == "rz":
            rz_virtual_update(phase, 0, args["zeta"])
        elif kind in ("rx_pi2", "rx_mpi2", "rx_pi"):
            pulse_id = pulse_id_for(kind, 0)
            entry = pulses_map[pulse_id]
            theta = float(entry.get("parameters", {}).get("theta", np.pi / 2.0))
            duration_ns = float(entry.get("calibration", {}).get("duration_ns", args.get("duration_ns", 25.0)))
            t, Iq, Qq, Jq, _, _ = run_rx(theta, duration_ns, 0, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[0], pulse_entry=entry)
            dt_local = t[1] - t[0] if t.size > 1 else dt
            t_blk = t if t_all.size == 0 else t_all[-1] + dt_local + t
            t_all = np.concatenate([t_all, t_blk]) if t_all.size else t_blk
            for q in [0, 1]:
                I_all[q] = np.concatenate([I_all[q], Iq[q]]) if I_all[q].size else Iq[q]
                Q_all[q] = np.concatenate([Q_all[q], Qq[q]]) if Q_all[q].size else Qq[q]
            J_all[LEGACY_EDGE] = np.concatenate([J_all[LEGACY_EDGE], Jq[LEGACY_EDGE]]) if J_all[LEGACY_EDGE].size else Jq[LEGACY_EDGE]
    dt_seq = dt if t_all.size < 2 else t_all[1] - t_all[0]
    H = build_H(t_all, dt_seq, I_all, Q_all, J_all, LEGACY_NQ, [LEGACY_EDGE]) if t_all.size else None
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    if H is None:
        rho_f = rho0
    else:
        res = simulate_legacy(t_all, H, rho0, dt_seq)
        rho_f = res.states[-1]
    if rho_f.isket:
        rho_f = qt.ket2dm(rho_f)
    return t_all, I_all, Q_all, J_all, rho_f


def print_params_table(rx_pi_params, rx_pi2_params, rz_rule, iswap_params):
    headers = ["gate", "duration_ns", "sigma_ns", "peak_rad_per_s", "theta_actual_rad", "expected_fidelity", "default_phase_rad"]
    rows = [
        [
            "rx_pi",
            rx_pi_params["duration_ns"],
            rx_pi_params["sigma_ns"],
            rx_pi_params["omega_peak_rad_per_s"],
            rx_pi_params["theta_actual_rad"],
            rx_pi_params.get("expected_fidelity") or "-",
            rx_pi_params["default_phase_rad"],
        ],
        [
            "rx_pi_over_2",
            rx_pi2_params["duration_ns"],
            rx_pi2_params["sigma_ns"],
            rx_pi2_params["omega_peak_rad_per_s"],
            rx_pi2_params["theta_actual_rad"],
            rx_pi2_params.get("expected_fidelity") or "-",
            rx_pi2_params["default_phase_rad"],
        ],
        ["rz_virtual", "-", "-", "-", "-", "-", rz_rule["calibration_reference_phase_rad"]],
        [
            "iswap",
            iswap_params["duration_ns"],
            "-",
            iswap_params["J_amp_rad_per_s"],
            iswap_params.get("theta_actual_rad", "-"),
            iswap_params.get("expected_fidelity") or "-",
            "-",
        ],
    ]
    colw = [max(len(str(h)), max(len(str(r[i])) for r in rows)) for i, h in enumerate(headers)]
    line = "  ".join(str(headers[i]).ljust(colw[i]) for i in range(len(headers)))
    print(line)
    print("-" * len(line))
    for r in rows:
        print("  ".join(str(r[i]).ljust(colw[i]) for i in range(len(headers))))


def step_fidelities_matched(blocks, pulses_map, dt_ns=DT_NS, sigma_frac=SIGMA_FRAC):
    dt = dt_ns * 1e-9
    phase = {0: 0.0, 1: 0.0}
    t_all = np.array([], dtype=float)
    I_all = {0: np.zeros(0), 1: np.zeros(0)}
    Q_all = {0: np.zeros(0), 1: np.zeros(0)}
    J_all = {LEGACY_EDGE: np.zeros(0)}
    U_cum = qt.qeye(2)
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    out = []
    for kind, args in blocks:
        if kind == "rz":
            zeta = float(args["zeta"])
            rz_virtual_update(phase, 0, zeta)
        elif kind in ("rx_pi2", "rx_mpi2", "rx_pi"):
            pulse_id = pulse_id_for(kind, 0)
            entry = pulses_map[pulse_id]
            theta = float(entry.get("parameters", {}).get("theta", np.pi / 2.0))
            duration_ns = float(entry.get("calibration", {}).get("duration_ns", args.get("duration_ns", 25.0)))
            t, Iq, Qq, Jq, _, _ = run_rx(theta, duration_ns, 0, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[0], pulse_entry=entry)
            dt_local = t[1] - t[0] if t.size > 1 else dt
            t_blk = t if t_all.size == 0 else t_all[-1] + dt_local + t
            t_all = t_blk if t_all.size == 0 else np.concatenate([t_all, t_blk])
            for q in [0, 1]:
                I_all[q] = Iq[q] if I_all[q].size == 0 else np.concatenate([I_all[q], Iq[q]])
                Q_all[q] = Qq[q] if Q_all[q].size == 0 else np.concatenate([Q_all[q], Qq[q]])
            J_all[LEGACY_EDGE] = Jq[LEGACY_EDGE] if J_all[LEGACY_EDGE].size == 0 else np.concatenate([J_all[LEGACY_EDGE], Jq[LEGACY_EDGE]])
            U_cum = xy_unitary(theta, phase[0]) * U_cum
        dt_seq = dt if t_all.size < 2 else t_all[1] - t_all[0]
        H = build_H(t_all, dt_seq, I_all, Q_all, J_all, LEGACY_NQ, [LEGACY_EDGE]) if t_all.size else None
        rho_sim = rho0 if H is None else simulate_legacy(t_all, H, rho0, dt_seq).states[-1]
        U_full = qt.tensor([U_cum, qt.qeye(2)])
        rho_id = U_full * rho0 * U_full.dag()
        out.append((kind, fidelity_dm(rho_sim, rho_id)))
    return out


def ideal_unitary_from_blocks(blocks, pulses_map):
    phase = 0.0
    U = qt.qeye(2)
    for kind, args in blocks:
        if kind == "rz":
            phase += float(args["zeta"])
        elif kind == "rx_pi2":
            entry = pulses_map[pulse_id_for("rx_pi2", 0)]
            theta = float(entry.get("parameters", {}).get("theta", np.pi / 2.0))
            U = xy_unitary(theta, phase) * U
        elif kind == "rx_mpi2":
            entry = pulses_map[pulse_id_for("rx_mpi2", 0)]
            theta = float(entry.get("parameters", {}).get("theta", -np.pi / 2.0))
            U = xy_unitary(theta, phase) * U
        elif kind == "rx_pi":
            entry = pulses_map[pulse_id_for("rx_pi", 0)]
            theta = float(entry.get("parameters", {}).get("theta", np.pi))
            U = xy_unitary(theta, phase) * U
    return U


def U_local(q, U2):
    return qt.tensor([U2, qt.qeye(2)]) if q == 0 else qt.tensor([qt.qeye(2), U2])


def U_xy(theta):
    _, _, _, _, X, Y, Z = single_qubit_ops()
    XX = lift_two(X, 0, X, 1, 2)
    YY = lift_two(Y, 0, Y, 1, 2)
    Hxy = 0.5 * (XX + YY)
    return (-1j * float(theta) * Hxy).expm()


def synth_su2_blocks(q, alpha, beta, gamma, dur_ns):
    return [
        ("rz", {"q": int(q), "zeta": float(alpha)}),
        ("rx_pi2", {"q": int(q), "duration_ns": float(dur_ns)}),
        ("rz", {"q": int(q), "zeta": float(beta)}),
        ("rx_mpi2", {"q": int(q), "duration_ns": float(dur_ns)}),
        ("rz", {"q": int(q), "zeta": float(gamma)}),
    ]


def simulate_sequence_2q(blocks, pulses_map, dt_ns=DT_NS, sigma_frac=SIGMA_FRAC):
    dt = dt_ns * 1e-9
    phase = {0: 0.0, 1: 0.0}
    t_all = np.array([], dtype=float)
    I_all = {0: np.zeros(0), 1: np.zeros(0)}
    Q_all = {0: np.zeros(0), 1: np.zeros(0)}
    J_all = {LEGACY_EDGE: np.zeros(0)}
    for kind, args in blocks:
        if kind == "rz":
            rz_virtual_update(phase, int(args["q"]), float(args["zeta"]))
        elif kind in ["rx_pi2", "rx_mpi2", "rx_pi"]:
            q = int(args["q"])
            entry = pulses_map[pulse_id_for(kind, q)]
            theta = float(entry.get("parameters", {}).get("theta", np.pi / 2.0))
            duration_ns = float(entry.get("calibration", {}).get("duration_ns", args.get("duration_ns", 25.0)))
            t, Iq, Qq, Jq, _, _ = run_rx(theta, duration_ns, q, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[q], pulse_entry=entry)
            dt_local = t[1] - t[0] if t.size > 1 else dt
            t_blk = t if t_all.size == 0 else t_all[-1] + dt_local + t
            t_all = t_blk if t_all.size == 0 else np.concatenate([t_all, t_blk])
            for qq in [0, 1]:
                I_all[qq] = Iq[qq] if I_all[qq].size == 0 else np.concatenate([I_all[qq], Iq[qq]])
                Q_all[qq] = Qq[qq] if Q_all[qq].size == 0 else np.concatenate([Q_all[qq], Qq[qq]])
            J_all[LEGACY_EDGE] = Jq[LEGACY_EDGE] if J_all[LEGACY_EDGE].size == 0 else np.concatenate([J_all[LEGACY_EDGE], Jq[LEGACY_EDGE]])
        elif kind == "iswap":
            qi, qj = args["edge"]
            entry = pulses_map[f"iswap_q{qi}_q{qj}"]
            theta = float(entry.get("parameters", {}).get("theta", 0.25 * np.pi))
            duration_ns = float(entry.get("calibration", {}).get("duration_ns", args.get("duration_ns", 50.0)))
            t, Iq, Qq, Jq, _, _ = run_iswap(duration_ns, dt_ns=dt_ns, pulse_entry=entry)
            dt_local = t[1] - t[0] if t.size > 1 else dt
            t_blk = t if t_all.size == 0 else t_all[-1] + dt_local + t
            t_all = t_blk if t_all.size == 0 else np.concatenate([t_all, t_blk])
            for qq in [0, 1]:
                I_all[qq] = Iq[qq] if I_all[qq].size == 0 else np.concatenate([I_all[qq], Iq[qq]])
                Q_all[qq] = Qq[qq] if Q_all[qq].size == 0 else np.concatenate([Q_all[qq], Qq[qq]])
            J_all[LEGACY_EDGE] = Jq[LEGACY_EDGE] if J_all[LEGACY_EDGE].size == 0 else np.concatenate([J_all[LEGACY_EDGE], Jq[LEGACY_EDGE]])
    dt_seq = dt if t_all.size < 2 else t_all[1] - t_all[0]
    H = build_H(t_all, dt_seq, I_all, Q_all, J_all, LEGACY_NQ, [LEGACY_EDGE]) if t_all.size else None
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    rho_f = rho0 if H is None else simulate_legacy(t_all, H, rho0, dt_seq).states[-1]
    if rho_f.isket:
        rho_f = qt.ket2dm(rho_f)
    return t_all, I_all, Q_all, J_all, rho_f


def ideal_unitary_from_blocks_2q(blocks, pulses_map):
    phase = {0: 0.0, 1: 0.0}
    U = qt.qeye([2, 2])
    for kind, args in blocks:
        if kind == "rz":
            q = int(args["q"])
            phase[q] += float(args["zeta"])
        elif kind in ["rx_pi2", "rx_mpi2", "rx_pi"]:
            q = int(args["q"])
            entry = pulses_map[pulse_id_for(kind, q)]
            ang = float(entry.get("parameters", {}).get("theta", np.pi / 2.0))
            Uloc = xy_unitary(ang, phase[q])
            U = U_local(q, Uloc) * U
        elif kind == "iswap":
            qi, qj = args["edge"]
            entry = pulses_map[f"iswap_q{qi}_q{qj}"]
            theta = float(entry.get("parameters", {}).get("theta", 0.25 * np.pi))
            U = U_xy(theta) * U
    return U


def random_zyz():
    return np.random.uniform(-np.pi, np.pi), np.random.uniform(0.0, np.pi), np.random.uniform(-np.pi, np.pi)


def random_su4_blocks(pulses_map, depth, single_dur_ns=25.0, ent_p=0.4):
    blocks = []
    for _ in range(depth):
        if np.random.rand() < ent_p:
            entry = pulses_map["iswap_q0_q1"]
            theta = float(entry.get("parameters", {}).get("theta", 0.25 * np.pi))
            duration_ns = float(entry.get("calibration", {}).get("duration_ns", 50.0))
            blocks.append(("iswap", {"edge": [0, 1], "theta": theta, "duration_ns": duration_ns}))
        else:
            q = int(np.random.choice([0, 1]))
            a, b, c = random_zyz()
            blocks.extend(synth_su2_blocks(q, a, b, c, single_dur_ns))
    return blocks


def test_random_su4(pulses_map, num_tests=5, depth=4, dt_ns=DT_NS, sigma_frac=SIGMA_FRAC):
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    out = []
    for _ in range(num_tests):
        blocks = random_su4_blocks(pulses_map, depth)
        _, _, _, _, rho_sim = simulate_sequence_2q(blocks, pulses_map, dt_ns=dt_ns, sigma_frac=sigma_frac)
        U_id = ideal_unitary_from_blocks_2q(blocks, pulses_map)
        rho_tar = U_id * rho0 * U_id.dag()
        F = fidelity_dm(rho_sim, rho_tar)
        out.append((blocks, F))
    return out


def run_legacy_validations(pulses):
    pulses_map = build_pulse_index(pulses)
    if "iswap_q0_q1" not in pulses_map:
        print("\nLegacy calibration diagnostics skipped: no two-qubit pulses present.")
        return
    entry_rx_pi_q0 = pulses_map["rx_q0_pi"]
    theta_rx_pi = float(entry_rx_pi_q0.get("parameters", {}).get("theta", np.pi))
    dur_rx_pi = float(entry_rx_pi_q0.get("calibration", {}).get("duration_ns", 25.0))
    t_rx_pi, I_rx_pi, Q_rx_pi, J_rx_pi, res_rx_pi, p_rx_pi = run_rx(theta_rx_pi, dur_rx_pi, 0, pulse_entry=entry_rx_pi_q0)

    entry_rx_pi2_q1 = pulses_map["rx_q1_pi_over_2"]
    theta_rx_pi2 = float(entry_rx_pi2_q1.get("parameters", {}).get("theta", np.pi / 2.0))
    dur_rx_pi2 = float(entry_rx_pi2_q1.get("calibration", {}).get("duration_ns", 25.0))
    t_rx_pi2, I_rx_pi2, Q_rx_pi2, J_rx_pi2, res_rx_pi2, p_rx_pi2 = run_rx(theta_rx_pi2, dur_rx_pi2, 1, pulse_entry=entry_rx_pi2_q1)

    entry_iswap = pulses_map["iswap_q0_q1"]
    dur_iswap = float(entry_iswap.get("calibration", {}).get("duration_ns", 50.0))
    t_sw, I_sw, Q_sw, J_sw, res_sw, p_sw = run_iswap(dur_iswap, pulse_entry=entry_iswap)

    entry_ry_q0 = pulses_map["rx_q0_pi_over_2"]
    entry_rx_pi2_q0 = entry_ry_q0
    theta_ry = float(entry_ry_q0.get("parameters", {}).get("theta", np.pi / 2.0))
    dur_ry = float(entry_ry_q0.get("calibration", {}).get("duration_ns", 25.0))
    t_ry, I_ry, Q_ry, J_ry, res_ry, p_ry = run_rx(theta_ry, dur_ry, 0, phase=np.pi / 2.0, pulse_entry=entry_ry_q0)

    print("\nLegacy calibration diagnostics:")
    e0, e1, *_ = single_qubit_ops()
    U_rx = xy_unitary(theta_rx_pi, 0.0)
    rho_target_rx = qt.ket2dm(qt.tensor([U_rx * e0, e0]))
    F_rx = fidelity_dm(res_rx_pi.states[-1], rho_target_rx)
    ry_state_q1 = xy_unitary(theta_rx_pi2, 0.0) * e0
    rx_pi2_state_q1 = (ry_state_q1).unit()
    rho_target_rx2 = qt.ket2dm(qt.tensor([e0, rx_pi2_state_q1]))
    F_rx2 = fidelity_dm(res_rx_pi2.states[-1], rho_target_rx2)
    psi10 = qt.tensor([e1, e0])
    theta_iswap = float(entry_iswap.get("parameters", {}).get("theta", 0.25 * np.pi))
    U_is = U_xy(theta_iswap)
    rho_target_sw = U_is * qt.ket2dm(psi10) * U_is.dag()
    F_sw = fidelity_dm(res_sw.states[-1], rho_target_sw)
    U_ry = ry_unitary(theta_ry)
    rho_target_ry = qt.ket2dm(qt.tensor([U_ry * e0, e0]))
    F_ry = fidelity_dm(res_ry.states[-1], rho_target_ry)
    alpha = np.random.uniform(-np.pi, np.pi)
    beta = np.random.uniform(0.0, np.pi)
    gamma = np.random.uniform(-np.pi, np.pi)
    entry_rx_mpi2_q0 = pulses_map["rx_q0_neg_pi_over_2"]
    dur_pi2_q0 = float(entry_rx_pi2_q0.get("calibration", {}).get("duration_ns", 25.0))
    dur_mpi2_q0 = float(entry_rx_mpi2_q0.get("calibration", {}).get("duration_ns", 25.0))
    blocks = [
        ("rz", {"zeta": float(alpha)}),
        ("rx_pi2", {"duration_ns": dur_pi2_q0}),
        ("rz", {"zeta": float(beta)}),
        ("rx_mpi2", {"duration_ns": dur_mpi2_q0}),
        ("rz", {"zeta": float(gamma)}),
    ]
    _, _, _, _, rho_su2 = simulate_sequence_on_q0(blocks, pulses_map, dt_ns=DT_NS, sigma_frac=SIGMA_FRAC)
    rho0q = qt.ket2dm(qt.tensor([e0, e0]))
    U_blocks = ideal_unitary_from_blocks(blocks, pulses_map)
    U_full_seq = qt.tensor([U_blocks, qt.qeye(2)])
    rho_seq_target = U_full_seq * rho0q * U_full_seq.dag()
    F_su2_frame_virtual = fidelity_dm(rho_su2, rho_seq_target)
    U_ZYZ = su2_from_zyz(alpha, beta, gamma)
    U_full_ZYZ = qt.tensor([U_ZYZ, qt.qeye(2)])
    rho_ZYZ = U_full_ZYZ * rho0q * U_full_ZYZ.dag()
    F_su2_physical_ZYZ = fidelity_dm(rho_su2, rho_ZYZ)
    rz_rule = {
        "type": "virtual_phase_advance",
        "rule": "I=Omega*cos(phi), Q=Omega*sin(phi); Rz(zeta): phi:=phi+zeta",
        "calibration_reference_phase_rad": 0.0,
    }
    print_params_table(p_rx_pi, p_rx_pi2, rz_rule, p_sw)
    print(f"Fidelity rx_pi on q0 to |1>: {F_rx}")
    print(f"Fidelity rx_pi_over_2 on q1 to Rx(pi/2)|0>: {F_rx2}")
    print("Fidelity sqrt(iSWAP) on |10> to (|10|-i|01>)/sqrt2:", F_sw)
    print("Fidelity ry_pi_over_2 on q0 to Ry(pi/2)|0>:", F_ry)
    step_F = step_fidelities_matched(blocks, pulses_map, dt_ns=DT_NS, sigma_frac=SIGMA_FRAC)
    cum = []
    for kind, Fk in step_F:
        cum.append(kind)
        print("Fidelity after " + " + ".join(cum) + ":", Fk)
    print("Fidelity SU2 sequence against frame-virtual target:", F_su2_frame_virtual)
    print("Fidelity SU2 sequence against physical ZYZ target:", F_su2_physical_ZYZ)
    rnd_tests = test_random_su4(pulses_map, num_tests=5, depth=4, dt_ns=DT_NS, sigma_frac=SIGMA_FRAC)
    for idx, (blk, F) in enumerate(rnd_tests, start=1):
        print("Random SU4 test", idx, "fidelity:", F)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__ or "")
    parser.add_argument("--devicelib", type=Path, default=DEVICELIB_PATH, help="Source Rigetti devicelib JSON")
    parser.add_argument("--num-qubits", type=int, help="Number of qubits to include (<= devicelib size)")
    parser.add_argument("--qubits", type=str, help="Comma-separated list of physical qubit indices to include")
    parser.add_argument("--out-device", type=Path, help="Output path for generated device JSON")
    parser.add_argument("--out-pulses", type=Path, help="Output path for generated pulse JSON")
    parser.add_argument(
        "--no-coupling",
        action="store_true",
        help="Drop all two-qubit couplings and related iswap pulses (single-qubit only).",
    )
    parser.add_argument(
        "--ideal-gates",
        action="store_true",
        help="Generate idealized pulse/device configs with unity fidelities and zero gate errors.",
    )
    args = parser.parse_args()

    devicelib_doc = load_json(args.devicelib)
    source_qubits = int(devicelib_doc.get("num_qubits", 0))
    if source_qubits <= 0:
        raise ValueError("Devicelib does not specify a positive number of qubits")

    if args.qubits:
        try:
            selected_qubits = [int(tok.strip()) for tok in args.qubits.split(",") if tok.strip()]
        except ValueError as exc:
            raise ValueError("--qubits must be a comma-separated list of integers") from exc
        if not selected_qubits:
            raise ValueError("--qubits list cannot be empty")
        if len(set(selected_qubits)) != len(selected_qubits):
            raise ValueError("--qubits list contains duplicates")
        if any(q < 0 or q >= source_qubits for q in selected_qubits):
            raise ValueError(f"--qubits entries must be within [0, {source_qubits - 1}]")
    else:
        if args.num_qubits is not None:
            if args.num_qubits <= 0 or args.num_qubits > source_qubits:
                raise ValueError(f"--num-qubits must be between 1 and {source_qubits}")
            selected_qubits = list(range(args.num_qubits))
        else:
            selected_qubits = list(range(source_qubits))

    mapping = {old: new for new, old in enumerate(selected_qubits)}
    filtered_devicelib = filter_devicelib(devicelib_doc, selected_qubits, mapping)
    if args.no_coupling:
        filtered_devicelib["coupling"] = []
        for field in ("gate_lens", "gate_errs"):
            if field in filtered_devicelib and isinstance(filtered_devicelib[field], dict):
                filtered_devicelib[field] = drop_two_qubit_gate_entries(filtered_devicelib[field])

    device_doc, gate_lens, gate_errs = build_device_document(filtered_devicelib)

    if args.ideal_gates:
        zero_errs = {key: 0.0 for key in gate_errs}
        gate_errs = zero_errs
        device_doc["gate_errs"] = zero_errs
        metadata = device_doc.setdefault("metadata", {})
        meta_errs = metadata.get("gate_errs")
        if isinstance(meta_errs, dict):
            metadata["gate_errs"] = {key: 0.0 for key in meta_errs}
        else:
            metadata["gate_errs"] = zero_errs.copy()

    num_qubits = len(selected_qubits)
    calibrated_name = f"rigetti_ankaa{num_qubits}q_calibrated"
    device_doc["name"] = calibrated_name
    device_doc["num_qubits"] = num_qubits
    device_doc["physical_qubits"] = list(range(num_qubits))
    device_doc["cx_coupling"] = filtered_devicelib.get("coupling", [])

    metadata = device_doc.get("metadata", {})
    metadata["num_qubits"] = num_qubits
    metadata["coupling"] = filtered_devicelib.get("coupling", [])
    metadata["basis_gates"] = filtered_devicelib.get("basis_gates", metadata.get("basis_gates", []))
    if args.no_coupling:
        device_doc["basis_gates"] = [g for g in device_doc.get("basis_gates", []) if "iswap" not in str(g).lower()]
        metadata["basis_gates"] = [g for g in metadata.get("basis_gates", []) if "iswap" not in str(g).lower()]

    pulses, fidelity_map = build_pulse_library(
        device_doc["metadata"],
        gate_lens,
        gate_errs,
        DT_NS,
        SIGMA_FRAC,
        FIDELITY_TOL,
        ideal=args.ideal_gates,
    )
    pulse_doc = compose_pulse_document(device_doc["metadata"], pulses, calibrated_name)

    default_prefix = f"rigetti_ankaa{num_qubits}q"
    if args.out_device is not None:
        out_device_path = args.out_device
    elif num_qubits == source_qubits and DEVICELIB_PATH == args.devicelib:
        out_device_path = OUT_DEVICE_PATH
    else:
        out_device_path = ROOT / "data" / "devices" / f"{default_prefix}_device.json"

    if args.out_pulses is not None:
        out_pulse_path = args.out_pulses
    elif num_qubits == source_qubits and DEVICELIB_PATH == args.devicelib:
        out_pulse_path = OUT_PULSE_PATH
    else:
        out_pulse_path = ROOT / "data" / "devices" / f"{default_prefix}_pulses.json"

    dump_json(out_device_path, device_doc)
    dump_json(out_pulse_path, pulse_doc)
    print(f"Wrote device config to {out_device_path}")
    print(f"Wrote pulse template to {out_pulse_path}")
    print("Fidelity summary:")
    for key, fidelity in sorted(fidelity_map.items()):
        print(f"  {key}: {fidelity:.6f}")
    print("\nPulse Fidelity Table:")
    print_fidelity_table(pulses)
    if num_qubits >= 2:
        run_legacy_validations(pulses)


if __name__ == "__main__":
    main()
