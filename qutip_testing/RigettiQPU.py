import json
import argparse
import numpy as np
import qutip as qt
import math
import sys

def grid_edges_3x3():
    edges = []
    for r in range(3):
        for c in range(3):
            q = 3 * r + c
            if c < 2:
                edges.append((q, q + 1))
            if r < 2:
                edges.append((q, q + 3))
    return [tuple(sorted(e)) for e in edges]

def single_qubit_ops():
    e0 = qt.basis(2, 0)
    e1 = qt.basis(2, 1)
    s01 = e0 * e1.dag()
    s10 = s01.dag()
    X = s01 + s10
    Y = -1j * s01 + 1j * s10
    return e0, e1, s01, s10, X, Y

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

def zero_array(N):
    return np.zeros(N, dtype=float)

def resample_to_length(x, L):
    if x.size == L:
        return x.copy()
    if x.size == 0:
        return np.zeros(L, dtype=float)
    idx = np.linspace(0, x.size - 1, L)
    return np.interp(idx, np.arange(x.size), x)

def gaussian_centered(width, dt):
    n = int(round(width / dt))
    n = max(n, 1)
    t = np.arange(n, dtype=float) * dt
    t0 = 0.5 * width
    sigma = max(1e-16, width * 0.2)
    g = np.exp(-0.5 * ((t - t0) / sigma) ** 2)
    return g

def parse_calibration(path):
    with open(path, "r") as f:
        cal = json.load(f)
    dt_ns = cal.get("global", {}).get("dt_ns", 0.25)
    iswap = cal.get("two_qubit", {}).get("iswap_01", {})
    J_area = float(iswap.get("J_area_rad", math.pi * 0.25))
    rz_rule = cal.get("virtual_z", {"type": "virtual_phase_advance"})
    return dt_ns, J_area, rz_rule

def parse_program(path):
    with open(path, "r") as f:
        prog = json.load(f)
    backend = prog.get("backend", {})
    pulse_lib = prog.get("pulse_library", [])
    schedule = prog.get("schedule", [])
    return backend, pulse_lib, schedule

def build_pulse_library_map(pulse_lib):
    out = {}
    for p in pulse_lib:
        pid = p.get("id")
        if pid is None:
            continue
        out[pid] = p
    return out

def infer_dt_from_samples(pulse_lib):
    for p in pulse_lib:
        si = np.asarray(p.get("samples_i", []), dtype=float)
        sq = np.asarray(p.get("samples_q", []), dtype=float)
        sj = np.asarray(p.get("samples_j", []), dtype=float)
        width = float(p.get("width", 0.0))
        L = 0
        if si.size:
            L = si.size
        elif sq.size:
            L = sq.size
        elif sj.size:
            L = sj.size
        if L > 1 and width > 0.0:
            return width / float(L)
    return None

def total_duration_from_schedule(schedule, fallback=0.0):
    tmax = 0.0
    for s in schedule:
        t0 = float(s.get("start_time", 0.0))
        dur = float(s.get("duration", 0.0))
        tmax = max(tmax, t0 + dur)
    return tmax if tmax > 0.0 else fallback

def compile_envelopes(backend, pulse_map, schedule, dt, N, num_qubits=9):
    I = {q: zero_array(N) for q in range(num_qubits)}
    Q = {q: zero_array(N) for q in range(num_qubits)}
    edges = grid_edges_3x3()
    J = {e: zero_array(N) for e in edges}
    phase = {q: 0.0 for q in range(num_qubits)}
    for s in sorted(schedule, key=lambda x: int(x.get("index", 0))):
        pid = s.get("pulse_id")
        if pid is None or pid not in pulse_map:
            continue
        spec = pulse_map[pid]
        gate = str(s.get("gate", spec.get("gate", ""))).lower()
        qubits = s.get("qubits", spec.get("qubits", []))
        t0 = float(s.get("start_time", 0.0))
        dur = float(s.get("duration", spec.get("width", 0.0)))
        k0 = max(0, int(round(t0 / dt)))
        k1 = min(N, int(round((t0 + dur) / dt)))
        if k1 <= k0:
            continue
        wf_type = str(spec.get("waveform_type", spec.get("shape", ""))).lower()
        pars_all = {}
        pars_all.update(spec.get("parameters", {}) if isinstance(spec.get("parameters", {}), dict) else {})
        pars_all.update(s.get("parameters", {}) if isinstance(s.get("parameters", {}), dict) else {})
        if len(qubits) == 1:
            q = int(qubits[0])
            if q < 0 or q >= num_qubits:
                raise ValueError("qubit index out of range")
            if wf_type == "virtual" or spec.get("virtual", False) or gate == "rz":
                zeta = float(pars_all.get("zeta", pars_all.get("zeta_rad", 0.0)))
                phase[q] = phase.get(q, 0.0) + zeta
                continue
            si = np.asarray(spec.get("samples_i", []), dtype=float)
            sq = np.asarray(spec.get("samples_q", []), dtype=float)
            if si.size or sq.size:
                Li = si.size if si.size else sq.size
                Li = max(Li, 1)
                segL = k1 - k0
                if si.size == 0:
                    si = np.zeros(Li, dtype=float)
                if sq.size == 0:
                    sq = np.zeros(Li, dtype=float)
                si = resample_to_length(si, segL)
                sq = resample_to_length(sq, segL)
                I[q][k0:k1] += si
                Q[q][k0:k1] += sq
            else:
                if wf_type == "gaussian":
                    width = float(spec.get("width", dur))
                    g = gaussian_centered(width, dt)
                    area = float(pars_all.get("theta", pars_all.get("theta_rad", 0.0)))
                    if area == 0.0:
                        amp = float(spec.get("amplitude", 0.0))
                        g = amp * np.ones_like(g)
                        area = float(np.trapezoid(g, dx=dt))
                    else:
                        g = g * (area / max(1e-32, float(np.trapezoid(g, dx=dt))))
                    segL = k1 - k0
                    g = resample_to_length(g, segL)
                    phi = float(pars_all.get("phase", pars_all.get("phase_rad", 0.0))) + phase[q]
                    I[q][k0:k1] += g * np.cos(phi)
                    Q[q][k0:k1] += g * np.sin(phi)
                else:
                    amp = float(spec.get("amplitude", pars_all.get("amplitude", 0.0)))
                    segL = k1 - k0
                    phi = float(pars_all.get("phase", pars_all.get("phase_rad", 0.0))) + phase[q]
                    g = amp * np.ones(segL, dtype=float)
                    I[q][k0:k1] += g * np.cos(phi)
                    Q[q][k0:k1] += g * np.sin(phi)
        elif len(qubits) == 2:
            i, j = sorted((int(qubits[0]), int(qubits[1])))
            if (i, j) not in J:
                continue
            sj = np.asarray(spec.get("samples_j", []), dtype=float)
            if sj.size:
                segL = k1 - k0
                sj = resample_to_length(sj, segL)
                J[(i, j)][k0:k1] += sj
            else:
                theta = float(pars_all.get("theta", pars_all.get("theta_rad", 0.0)))
                amp = float(spec.get("amplitude", pars_all.get("amplitude", 0.0)))
                segL = k1 - k0
                if theta != 0.0 and dur > 0.0:
                    J[(i, j)][k0:k1] += (theta / dur)
                else:
                    J[(i, j)][k0:k1] += amp
    return I, Q, J

def build_hamiltonian(t, dt, I, Q, J, Nq):
    _, _, s01, s10, X, Y = single_qubit_ops()
    H = []
    for q in range(Nq):
        Ii = I[q]
        Qi = Q[q]
        if np.any(Ii):
            H.append([lift_one(X, q, Nq), coeff_from_array(0.5 * Ii, dt)])
        if np.any(Qi):
            H.append([lift_one(Y, q, Nq), coeff_from_array(0.5 * Qi, dt)])
    for e, arr in J.items():
        if np.any(arr):
            i, j = e
            exch = lift_two(s10, i, s01, j, Nq) + lift_two(s01, i, s10, j, Nq)
            H.append([exch, coeff_from_array(arr, dt)])
    return H

def initial_state_ket(spec, Nq):
    if spec is None:
        return qt.tensor([qt.basis(2, 0)] * Nq)
    if isinstance(spec, str):
        if len(spec) != Nq or any(c not in "01" for c in spec):
            raise ValueError("invalid initial_state string")
        kets = [qt.basis(2, int(c)) for c in spec]
        return qt.tensor(kets)
    if isinstance(spec, list):
        if len(spec) != Nq:
            raise ValueError("invalid initial_state list")
        kets = [qt.basis(2, int(v)) for v in spec]
        return qt.tensor(kets)
    raise ValueError("invalid initial_state")

def run_pulse_and_statevector(H, rho0, t, dt, store_states=False):
    opts = {"method": "bdf", "rtol": 1e-7, "atol": 1e-9, "nsteps": 300000, "max_step": dt, "progress_bar": False, "store_states": False}
    res_me = qt.mesolve(H, rho0, t, c_ops=[], e_ops=[], options=opts)
    psi0 = rho0 if rho0.isket else rho0.purify()
    res_se = qt.sesolve(H, psi0, t, e_ops=[], args=None, options={"method": "bdf", "rtol": 1e-7, "atol": 1e-9, "nsteps": 300000, "max_step": dt, "progress_bar": False, "store_states": False})
    rho_me = res_me.states[-1]
    if rho_me.isket:
        rho_me = qt.ket2dm(rho_me)
    psi_f = res_se.states[-1]
    if not psi_f.isket:
        psi_f = psi_f.purify()
    rho_sv = qt.ket2dm(psi_f)
    F = float(qt.fidelity(rho_me, rho_sv))
    return F, rho_me, rho_sv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--program_json", required=True)
    ap.add_argument("--calib_json", required=True)
    ap.add_argument("--initial_state", default=None)
    ap.add_argument("--out_json", default=None)
    args = ap.parse_args()
    dt_ns_cal, J_area_cal, rz_rule = parse_calibration(args.calib_json)
    backend, pulse_lib, schedule = parse_program(args.program_json)
    pulse_map = build_pulse_library_map(pulse_lib)
    dt_inf = infer_dt_from_samples(pulse_lib)
    dt_ns = float(backend.get("dt_ns", dt_ns_cal))
    if not dt_ns and dt_inf is not None:
        dt_ns = dt_inf * 1e9
    if not dt_ns:
        dt_ns = 0.25
    dt = dt_ns * 1e-9
    T_backend = float(backend.get("total_duration", 0.0))
    T = total_duration_from_schedule(schedule, fallback=T_backend)
    if T <= 0.0:
        print("0.0", flush=True)
        return
    N = int(round(T / dt)) + 1
    t = np.arange(N, dtype=float) * dt
    I, Q, J = compile_envelopes(backend, pulse_map, schedule, dt, N, num_qubits=9)
    H = build_hamiltonian(t, dt, I, Q, J, 9)
    psi0 = initial_state_ket(args.initial_state, 9) if args.initial_state is not None else initial_state_ket(None, 9)
    rho0 = qt.ket2dm(psi0)
    F, rho_me, rho_sv = run_pulse_and_statevector(H, rho0, t, dt, store_states=False)
    if args.out_json:
        out = {
            "run": {"dt_ns": dt_ns, "duration_s": T, "num_qubits": 9},
            "backend": backend,
            "fidelity_mesolve_vs_sesolve": F
        }
        with open(args.out_json, "w") as f:
            json.dump(out, f, indent=2)
    print(F, flush=True)

if __name__ == "__main__":
    main()

