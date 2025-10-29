#!/usr/bin/env python3
"""Produce RX/iSWAP pulse plots and a GRAPE-merged gate with smooth ramps and masked edges."""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from qutip_qtrl import pulseoptim


DEFAULT_DT_NS = 0.25
DEFAULT_SIGMA_FRAC = 0.2
DEFAULT_ISWAP_RAMP_FRAC = 0.1
DEFAULT_RX_DURATION_NS = 25.0
DEFAULT_ISWAP_DURATION_NS = 120.0
DEFAULT_MERGED_DURATION_NS = 120.0
DEFAULT_RX_THETA = 0.5 * math.pi
DEFAULT_ISWAP_THETA = 0.25 * math.pi
DEFAULT_MAX_ITER = 300
DEFAULT_MAX_WALL = 120
DEFAULT_SEED = None
DEFAULT_SINE_MODES = 5


def _time_axis_ns(slots: int, dt_ns: float) -> np.ndarray:
    slots = max(1, int(slots))
    return np.arange(slots, dtype=float) * dt_ns


def _to_ghz(amps_rad_s: np.ndarray) -> np.ndarray:
    return amps_rad_s / (2.0 * math.pi * 1e9)


def generate_rx_gaussian(theta: float, duration_ns: float, dt_ns: float, sigma_frac: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    num_slots = max(2, int(round(duration_ns / dt_ns)))
    times_ns = _time_axis_ns(num_slots, dt_ns)
    times_s = times_ns * 1e-9
    sigma = sigma_frac * duration_ns * 1e-9
    center = 0.5 * duration_ns * 1e-9
    envelope = np.exp(-0.5 * ((times_s - center) / sigma) ** 2)
    dt = dt_ns * 1e-9
    area = float(np.trapezoid(envelope, dx=dt))
    if abs(area) <= 1e-18:
        raise ValueError("RX Gaussian envelope area is negligible; adjust duration or dt.")
    samples = (theta / area) * envelope
    return times_ns, samples.astype(float), np.zeros_like(samples)


def generate_iswap_flattop(theta: float, duration_ns: float, dt_ns: float, ramp_frac: float) -> tuple[np.ndarray, np.ndarray]:
    num_slots = max(2, int(round(duration_ns / dt_ns)))
    times_ns = _time_axis_ns(num_slots, dt_ns)
    times_s = times_ns * 1e-9
    duration_s = duration_ns * 1e-9
    sigma = max(1e-12, ramp_frac * duration_s)
    offset = ramp_frac * duration_s
    erf = math.erf
    # Smooth turn-on/off using error-function ramps
    envelope = 0.5 * (1.0 + np.vectorize(erf)((times_s - offset) / (math.sqrt(2.0) * sigma)))
    envelope *= 0.5 * (1.0 + np.vectorize(erf)(((duration_s - offset) - times_s) / (math.sqrt(2.0) * sigma)))
    dt = dt_ns * 1e-9
    area = float(np.trapezoid(envelope, dx=dt))
    if abs(area) <= 1e-18:
        raise ValueError("iSWAP envelope area is negligible; adjust ramp parameters.")
    amplitude = theta / area
    samples = amplitude * envelope
    return times_ns, samples.astype(float)


def resample_waveform(samples: np.ndarray, source_duration_ns: float, target_slots: int, target_dt_ns: float) -> np.ndarray:
    target_slots = max(1, int(target_slots))
    if samples.size == 0:
        return np.zeros(target_slots, dtype=float)
    if target_slots == samples.size:
        return samples.astype(float)
    times_src = np.linspace(0.0, source_duration_ns, samples.size, dtype=float)
    times_tgt = np.linspace(0.0, source_duration_ns, target_slots, dtype=float)
    return np.interp(times_tgt, times_src, samples, left=0.0, right=0.0).astype(float)


def random_sine_series(n_slots: int, n_ctrls: int, rng: np.random.Generator, modes: int, scale: float) -> np.ndarray:
    if n_slots <= 0:
        return np.zeros((0, n_ctrls), dtype=float)
    modes = max(1, modes)
    x = np.linspace(0.0, math.pi, n_slots + 2, dtype=float)[1:-1]
    waves = np.zeros((n_slots, n_ctrls), dtype=float)
    for ctrl in range(n_ctrls):
        wave = np.zeros(n_slots, dtype=float)
        for mode in range(1, modes + 1):
            coeff = rng.normal(scale=scale / mode)
            wave += coeff * np.sin(mode * x)
        waves[:, ctrl] = wave
    return waves


def sine_basis(n_slots: int, modes: int) -> np.ndarray:
    modes = max(1, modes)
    if n_slots <= 0:
        return np.zeros((0, modes), dtype=float)
    x = np.linspace(0.0, math.pi, n_slots + 2, dtype=float)[1:-1]
    basis = np.zeros((n_slots, modes), dtype=float)
    for mode in range(1, modes + 1):
        basis[:, mode - 1] = np.sin(mode * x)
    return basis


def fit_to_sine_series(data: np.ndarray, modes: int) -> np.ndarray:
    if data.size == 0:
        return data
    n_slots, n_ctrls = data.shape
    basis = sine_basis(n_slots, modes)
    if basis.size == 0:
        return data
    result = np.zeros_like(data)
    for ctrl in range(n_ctrls):
        coeffs, _, _, _ = np.linalg.lstsq(basis, data[:, ctrl], rcond=None)
        result[:, ctrl] = basis @ coeffs
    return result


def pauli_controls() -> tuple[list[qt.Qobj], list[str]]:
    sx = qt.sigmax()
    sy = qt.sigmay()
    ident = qt.qeye(2)
    ix = 0.5 * qt.tensor(sx, ident)
    iy = 0.5 * qt.tensor(sy, ident)
    coupling = 0.5 * (qt.tensor(sx, sx) + qt.tensor(sy, sy))
    return [ix, iy, coupling], ["I0", "Q0", "J01"]


def sequential_unitary(rx_theta: float, iswap_theta: float) -> qt.Qobj:
    sx = qt.sigmax()
    sy = qt.sigmay()
    ident = qt.qeye(2)
    U_rx = qt.tensor((-1j * rx_theta * 0.5 * sx).expm(), ident)
    U_iswap = (-1j * iswap_theta * 0.5 * (qt.tensor(sx, sx) + qt.tensor(sy, sy))).expm()
    return U_iswap * U_rx


def unitary_from_controls(controls: np.ndarray, control_ops: list[qt.Qobj], dt: float) -> qt.Qobj:
    U = qt.qeye([2, 2])
    for row in controls:
        H = qt.Qobj(np.zeros((4, 4), dtype=complex), dims=[[2, 2], [2, 2]])
        for amp, op in zip(row, control_ops):
            H = H + float(amp) * op
        U = (-1j * H * dt).expm() * U
    return U


def run_grape_attempt(
    target: qt.Qobj,
    dt_ns: float,
    merged_duration_ns: float,
    amp_bound: float | None,
    max_iter: int,
    max_wall_time: int,
    rx_amp_ref: float,
    iswap_amp_ref: float,
    rx_wave: np.ndarray,
    rx_duration_ns: float,
    iswap_wave: np.ndarray,
    iswap_duration_ns: float,
    sine_modes: int,
    seed: int | None,
) -> tuple[np.ndarray, float, float, list[str]]:
    ops, labels = pauli_controls()
    n_ctrls = len(ops)
    n_ts = max(2, int(round(merged_duration_ns / dt_ns)))
    evo_time = merged_duration_ns * 1e-9
    dt = dt_ns * 1e-9

    rng = np.random.default_rng(seed)
    base = np.zeros((n_ts, n_ctrls), dtype=float)
    modes = min(max(1, sine_modes), n_ts)
    total_native_ns = rx_duration_ns + iswap_duration_ns
    scale = merged_duration_ns / total_native_ns if total_native_ns > 0 else 1.0
    if n_ts > 1:
        rx_alloc_ns = max(dt_ns, rx_duration_ns * scale)
        rx_slots = max(1, min(n_ts - 1, int(round(rx_alloc_ns / dt_ns))))
        iswap_slots = max(1, n_ts - rx_slots)
        rx_resampled = resample_waveform(rx_wave, rx_duration_ns, rx_slots, dt_ns)
        iswap_resampled = resample_waveform(iswap_wave, iswap_duration_ns, iswap_slots, dt_ns)
        base[:rx_slots, 0] = rx_resampled
        base[rx_slots:rx_slots + iswap_slots, 2] = iswap_resampled[:iswap_slots]
    else:
        base[0, 0] = rx_wave[0] if rx_wave.size else 0.0

    noise_scale = max(rx_amp_ref, iswap_amp_ref, 1.0)
    sine_component = random_sine_series(n_ts, n_ctrls, rng, modes, noise_scale)
    guess_core = base + 0.2 * sine_component

    amp_lbounds: list[float] = []
    amp_ubounds: list[float] = []
    user_bound = None if amp_bound is None else abs(float(amp_bound))
    for label in labels:
        if label.startswith("J"):
            ref = max(iswap_amp_ref, 1e-3)
            upper = user_bound if user_bound else 10.0 * ref
            amp_lbounds.append(-upper)
            amp_ubounds.append(upper)
        else:
            ref = max(rx_amp_ref, 1e-3)
            upper = user_bound if user_bound else 10.0 * ref
            amp_lbounds.append(-upper)
            amp_ubounds.append(upper)

    optimizer = pulseoptim.create_pulse_optimizer(
        qt.Qobj(np.zeros((4, 4), dtype=complex), dims=[[2, 2], [2, 2]]),
        ops,
        qt.qeye([2, 2]),
        target,
        num_tslots=n_ts,
        evo_time=evo_time,
        amp_lbound=amp_lbounds,
        amp_ubound=amp_ubounds,
        fid_err_targ=1e-6,
        max_iter=max_iter,
        max_wall_time=max_wall_time,
        dyn_type="UNIT",
        optim_method="FMIN_L_BFGS_B",
    )

    optimizer.dynamics.initialize_controls(guess_core.copy())
    result = optimizer.run_optimization()
    print(f"GRAPE reported fid_err: {float(np.real(result.fid_err)):.6e}")
    controls = fit_to_sine_series(np.array(result.final_amps, dtype=float), modes)

    U = unitary_from_controls(controls, ops, dt)
    overlap = (target.dag() * U).tr()
    fidelity = float(abs(overlap) / target.shape[0])
    return controls, dt, fidelity, labels


def run_grape(
    target: qt.Qobj,
    dt_ns: float,
    merged_duration_ns: float,
    amp_bound: float | None,
    max_iter: int,
    max_wall_time: int,
    seed: int | None,
    restarts: int,
    rx_amp_ref: float,
    iswap_amp_ref: float,
    rx_wave: np.ndarray,
    rx_duration_ns: float,
    iswap_wave: np.ndarray,
    iswap_duration_ns: float,
    sine_modes: int,
) -> tuple[np.ndarray, float, float, list[str]]:
    rng = np.random.default_rng(seed)
    best_controls: np.ndarray | None = None
    best_dt = dt_ns * 1e-9
    best_fidelity = -1.0
    labels: list[str] | None = None
    for attempt in range(max(1, restarts)):
        attempt_seed = rng.integers(0, 2**32 - 1)
        controls, dt, fidelity, labels_local = run_grape_attempt(
            target,
            dt_ns,
            merged_duration_ns,
            amp_bound,
            max_iter,
            max_wall_time,
            rx_amp_ref,
            iswap_amp_ref,
            rx_wave,
            rx_duration_ns,
            iswap_wave,
            iswap_duration_ns,
            sine_modes,
            attempt_seed,
        )
        if fidelity > best_fidelity:
            best_fidelity = fidelity
            best_controls = controls
            best_dt = dt
            labels = labels_local
    assert best_controls is not None and labels is not None
    return best_controls, best_dt, best_fidelity, labels


def plot_single(
    times_ns: np.ndarray,
    samples: np.ndarray,
    output_path: Path,
    fmt: str,
    dpi: int,
    x_ticks: list[float],
) -> None:
    plt.figure(figsize=(7.5, 5.0))
    amps = _to_ghz(samples)
    max_abs = float(np.max(np.abs(amps))) if amps.size else 0.0
    denom = max(max_abs, 1e-12)
    scaled = 0.8 * (amps / denom)
    plt.plot(times_ns, scaled)
    base_size = float(plt.rcParams.get("font.size", 10.0))
    label_size = base_size * 2.5
    plt.xlabel("Time (ns)", fontweight="bold", fontsize=label_size)
    plt.ylabel("Amplitude (GHz)", fontweight="bold", fontsize=label_size)
    plt.tick_params(axis="both", labelsize=label_size)
    plt.grid(alpha=0.3)
    plt.ylim(-0.5, 1.0)
    plt.yticks([-0.5, 0.0, 0.5, 1.0], labels=["-0.5", "0", "0.5", "1"], fontsize=label_size)
    if x_ticks:
        plt.xticks(x_ticks, [str(int(t)) for t in x_ticks], fontsize=label_size)
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.tight_layout()
    plt.savefig(output_path.with_suffix(f".{fmt}"), dpi=dpi, format=fmt)
    plt.close()


def plot_controls(
    times_ns: np.ndarray,
    controls: np.ndarray,
    output_path: Path,
    fmt: str,
    dpi: int,
    x_ticks: list[float],
) -> None:
    plt.figure(figsize=(7.5, 5.0))
    amps = _to_ghz(controls)
    max_abs = float(np.max(np.abs(amps))) if amps.size else 0.0
    denom = max(max_abs, 1e-12)
    scaled = 0.8 * (amps / denom)
    for idx in range(scaled.shape[1]):
        plt.plot(times_ns, scaled[:, idx])
    base_size = float(plt.rcParams.get("font.size", 10.0))
    label_size = base_size * 2.5
    plt.xlabel("Time (ns)", fontweight="bold", fontsize=label_size)
    plt.ylabel("Amplitude (GHz)", fontweight="bold", fontsize=label_size)
    plt.tick_params(axis="both", labelsize=label_size)
    plt.grid(alpha=0.3)
    plt.ylim(-0.5, 1.0)
    plt.yticks([-0.5, 0.0, 0.5, 1.0], labels=["-0.5", "0", "0.5", "1"], fontsize=label_size)
    if x_ticks:
        plt.xticks(x_ticks, [str(int(t)) for t in x_ticks], fontsize=label_size)
    plt.axhline(0.0, color="black", linewidth=0.8)
    plt.tight_layout()
    plt.savefig(output_path.with_suffix(f".{fmt}"), dpi=dpi, format=fmt)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot RX/iSWAP pulses and a GRAPE-merged pulse with smooth ramps and masked edges.",
    )
    parser.add_argument("--dt-ns", type=float, default=DEFAULT_DT_NS, help="Time step per slot (ns).")
    parser.add_argument("--sigma-frac", type=float, default=DEFAULT_SIGMA_FRAC, help="Gaussian sigma fraction for RX.")
    parser.add_argument("--iswap-ramp-frac", type=float, default=DEFAULT_ISWAP_RAMP_FRAC, help="Ramp fraction for iSWAP Gaussian edges.")
    parser.add_argument("--rx-duration-ns", type=float, default=DEFAULT_RX_DURATION_NS, help="RX pulse duration (ns).")
    parser.add_argument("--iswap-duration-ns", type=float, default=DEFAULT_ISWAP_DURATION_NS, help="iSWAP pulse duration (ns).")
    parser.add_argument("--merged-duration-ns", type=float, default=DEFAULT_MERGED_DURATION_NS, help="Desired merged gate duration (ns).")
    parser.add_argument("--rx-theta", type=float, default=DEFAULT_RX_THETA, help="RX rotation angle (rad).")
    parser.add_argument("--iswap-theta", type=float, default=DEFAULT_ISWAP_THETA, help="iSWAP exchange angle (rad).")
    parser.add_argument("--amp-bound", type=float, help="Optional absolute amplitude bound for all controls.")
    parser.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER, help="GRAPE maximum iterations.")
    parser.add_argument("--max-wall-time", type=int, default=DEFAULT_MAX_WALL, help="GRAPE wall-time limit (s).")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Random seed for reproducibility.")
    parser.add_argument("--restarts", type=int, default=3, help="Number of random GRAPE restarts (best result kept).")
    parser.add_argument("--sine-modes", type=int, default=DEFAULT_SINE_MODES, help="Number of sine modes for initialisation and projection.")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent, help="Directory for output figures.")
    parser.add_argument("--output-format", type=str, default="pdf", choices=["pdf", "png"], help="Figure output format.")
    parser.add_argument("--dpi", type=int, default=200, help="Figure DPI (for raster formats).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.seed is not None:
        np.random.seed(args.seed)

    times_rx_ns, rx_i_samples, _ = generate_rx_gaussian(
        args.rx_theta,
        args.rx_duration_ns,
        args.dt_ns,
        args.sigma_frac,
    )
    times_iswap_ns, iswap_samples = generate_iswap_flattop(
        args.iswap_theta,
        args.iswap_duration_ns,
        args.dt_ns,
        args.iswap_ramp_frac,
    )

    target = sequential_unitary(args.rx_theta, args.iswap_theta)

    rx_amp_ref = float(np.max(np.abs(rx_i_samples))) if rx_i_samples.size else 1.0
    iswap_amp_ref = float(np.max(np.abs(iswap_samples))) if iswap_samples.size else 1.0

    controls, dt, fidelity, labels = run_grape(
        target,
        args.dt_ns,
        args.merged_duration_ns,
        args.amp_bound,
        args.max_iter,
        args.max_wall_time,
        args.seed,
        args.restarts,
        rx_amp_ref,
        iswap_amp_ref,
        rx_i_samples,
        args.rx_duration_ns,
        iswap_samples,
        args.iswap_duration_ns,
        args.sine_modes,
    )

    times_ctrl_ns = _time_axis_ns(controls.shape[0], args.dt_ns)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    fmt = args.output_format
    dpi = args.dpi
    plot_single(times_rx_ns, rx_i_samples, output_dir / "rx_pulse", fmt, dpi, [0.0, args.rx_duration_ns])
    plot_single(times_iswap_ns, iswap_samples, output_dir / "iswap_pulse", fmt, dpi, [0.0, args.iswap_duration_ns / 2.0, args.iswap_duration_ns])
    plot_controls(times_ctrl_ns, controls, output_dir / "grape_combined_pulse", fmt, dpi, [0.0, args.merged_duration_ns / 2.0, args.merged_duration_ns])

    print(f"Plots written to {output_dir.resolve()}")
    print(f"Masked GRAPE fidelity: {fidelity:.6f}")
    print(f"Effective dt: {dt * 1e9:.6f} ns")


if __name__ == "__main__":
    main()
