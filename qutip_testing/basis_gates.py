import json
import numpy as np
import qutip as qt

CALIB_PATH = "calibration_params.json"
OUT_JSON_PATH = "pulse_params_2q.json"
NQ = 2
EDGE = (0, 1)

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

def simulate(t, H, rho0, dt):
    opts = {"method": "bdf", "rtol": 1e-7, "atol": 1e-9, "nsteps": 300000, "max_step": dt, "progress_bar": False, "store_states": True}
    return qt.mesolve(H, rho0, t, c_ops=[], e_ops=[], options=opts)

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
    area = float(np.trapezoid(g, dx=dt))
    A = theta / area
    return t, A * g, A, sigma

def iq_from_omega_phase(omega_t, phase):
    I = omega_t * np.cos(phase)
    Q = omega_t * np.sin(phase)
    return I, Q

def run_rx(theta, duration_ns, target, dt_ns=0.25, sigma_frac=0.2, phase=0.0):
    dt = dt_ns * 1e-9
    duration = duration_ns * 1e-9
    t, omega, A_peak, sigma = omega_gaussian_for_rotation(theta, duration, dt, sigma_frac)
    Iq = {q: np.zeros_like(t) for q in range(NQ)}
    Qq = {q: np.zeros_like(t) for q in range(NQ)}
    I, Q = iq_from_omega_phase(omega, phase)
    Iq[target] = I
    Qq[target] = Q
    J = {EDGE: np.zeros_like(t)}
    H = build_H(t, dt, Iq, Qq, J, NQ, [EDGE])
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    res = simulate(t, H, rho0, dt)
    params = {
        "type": "gaussian",
        "target_qubit": int(target),
        "theta_rad": float(theta),
        "duration_ns": float(duration_ns),
        "dt_ns": float(dt_ns),
        "sigma_ns": float(sigma * 1e9),
        "omega_peak_rad_per_s": float(A_peak),
        "I_peak_rad_per_s": float(A_peak),
        "area_omega_rad": float(np.trapezoid(omega, dx=dt)),
        "default_phase_rad": float(phase)
    }
    return t, Iq, Qq, J, res, params

def run_iswap(duration_ns, dt_ns=0.25):
    dt = dt_ns * 1e-9
    duration = duration_ns * 1e-9
    t = make_time_grid(dt, duration)
    J_area = 0.25 * np.pi
    J_amp = J_area / duration
    J = {EDGE: np.full_like(t, J_amp)}
    Iq = {0: np.zeros_like(t), 1: np.zeros_like(t)}
    Qq = {0: np.zeros_like(t), 1: np.zeros_like(t)}
    H = build_H(t, dt, Iq, Qq, J, NQ, [EDGE])
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e1, e0]))
    res = simulate(t, H, rho0, dt)
    params = {
        "type": "square_iswap",
        "edge": [int(EDGE[0]), int(EDGE[1])],
        "duration_ns": float(duration_ns),
        "dt_ns": float(dt_ns),
        "J_amp_rad_per_s": float(J_amp),
        "J_area_rad": float(J_area)
    }
    return t, Iq, Qq, J, res, params

def run_xy_entangler(theta, duration_ns, dt_ns=0.25):
    dt = dt_ns * 1e-9
    duration = duration_ns * 1e-9
    t = make_time_grid(dt, duration)
    J_area = float(theta)
    J_amp = J_area / duration
    J = {EDGE: np.full_like(t, J_amp)}
    Iq = {0: np.zeros_like(t), 1: np.zeros_like(t)}
    Qq = {0: np.zeros_like(t), 1: np.zeros_like(t)}
    H = build_H(t, dt, Iq, Qq, J, NQ, [EDGE])
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    res = simulate(t, H, rho0, dt)
    params = {
        "type": "square_xy",
        "edge": [int(EDGE[0]), int(EDGE[1])],
        "duration_ns": float(duration_ns),
        "dt_ns": float(dt_ns),
        "J_amp_rad_per_s": float(J_amp),
        "J_area_rad": float(J_area)
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

def simulate_sequence_on_q0(blocks, dt_ns=0.25, sigma_frac=0.2):
    dt = dt_ns * 1e-9
    phase = {0: 0.0, 1: 0.0}
    t_all = np.array([], dtype=float)
    I_all = {0: np.zeros(0), 1: np.zeros(0)}
    Q_all = {0: np.zeros(0), 1: np.zeros(0)}
    J_all = {EDGE: np.zeros(0)}
    for kind, args in blocks:
        if kind == "rz":
            rz_virtual_update(phase, 0, args["zeta"])
        elif kind == "rx_pi2":
            t, Iq, Qq, Jq, _, _ = run_rx(np.pi / 2.0, args["duration_ns"], 0, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[0])
            t_blk = t if t_all.size == 0 else t_all[-1] + dt + t
            t_all = np.concatenate([t_all, t_blk]) if t_all.size else t_blk
            for q in [0, 1]:
                I_all[q] = np.concatenate([I_all[q], Iq[q]]) if I_all[q].size else Iq[q]
                Q_all[q] = np.concatenate([Q_all[q], Qq[q]]) if Q_all[q].size else Qq[q]
            J_all[EDGE] = np.concatenate([J_all[EDGE], Jq[EDGE]]) if J_all[EDGE].size else Jq[EDGE]
        elif kind == "rx_mpi2":
            t, Iq, Qq, Jq, _, _ = run_rx(-np.pi / 2.0, args["duration_ns"], 0, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[0])
            t_blk = t if t_all.size == 0 else t_all[-1] + dt + t
            t_all = np.concatenate([t_all, t_blk]) if t_all.size else t_blk
            for q in [0, 1]:
                I_all[q] = np.concatenate([I_all[q], Iq[q]]) if I_all[q].size else Iq[q]
                Q_all[q] = np.concatenate([Q_all[q], Qq[q]]) if Q_all[q].size else Qq[q]
            J_all[EDGE] = np.concatenate([J_all[EDGE], Jq[EDGE]]) if J_all[EDGE].size else Jq[EDGE]
        elif kind == "rx_pi":
            t, Iq, Qq, Jq, _, _ = run_rx(np.pi, args["duration_ns"], 0, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[0])
            t_blk = t if t_all.size == 0 else t_all[-1] + dt + t
            t_all = np.concatenate([t_all, t_blk]) if t_all.size else t_blk
            for q in [0, 1]:
                I_all[q] = np.concatenate([I_all[q], Iq[q]]) if I_all[q].size else Iq[q]
                Q_all[q] = np.concatenate([Q_all[q], Qq[q]]) if Q_all[q].size else Qq[q]
            J_all[EDGE] = np.concatenate([J_all[EDGE], Jq[EDGE]]) if J_all[EDGE].size else Jq[EDGE]
    H = build_H(t_all, dt, I_all, Q_all, J_all, NQ, [EDGE])
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    res = simulate(t_all, H, rho0, dt)
    rho_f = res.states[-1]
    if rho_f.isket:
        rho_f = qt.ket2dm(rho_f)
    return t_all, I_all, Q_all, J_all, rho_f

def print_params_table(rx_pi_params, rx_pi2_params, rz_rule, iswap_params):
    headers = ["gate", "duration_ns", "sigma_ns", "peak_rad_per_s", "area_rad", "default_phase_rad"]
    rows = []
    rows.append(["rx_pi", rx_pi_params["duration_ns"], rx_pi_params["sigma_ns"], rx_pi_params["omega_peak_rad_per_s"], rx_pi_params["area_omega_rad"], rx_pi_params["default_phase_rad"]])
    rows.append(["rx_pi_over_2", rx_pi2_params["duration_ns"], rx_pi2_params["sigma_ns"], rx_pi2_params["omega_peak_rad_per_s"], rx_pi2_params["area_omega_rad"], rx_pi2_params["default_phase_rad"]])
    rows.append(["rz_virtual", "-", "-", "-", "-", rz_rule["calibration_reference_phase_rad"]])
    rows.append(["iswap", iswap_params["duration_ns"], "-", iswap_params["J_amp_rad_per_s"], iswap_params["J_area_rad"], "-"])
    colw = [max(len(str(h)), max(len(str(r[i])) for r in rows)) for i, h in enumerate(headers)]
    line = "  ".join(str(headers[i]).ljust(colw[i]) for i in range(len(headers)))
    print(line)
    print("-" * len(line))
    for r in rows:
        print("  ".join(str(r[i]).ljust(colw[i]) for i in range(len(headers))))

def step_fidelities_matched(blocks, dt_ns=0.25, sigma_frac=0.2):
    dt = dt_ns * 1e-9
    phase = {0: 0.0, 1: 0.0}
    t_all = np.array([], dtype=float)
    I_all = {0: np.zeros(0), 1: np.zeros(0)}
    Q_all = {0: np.zeros(0), 1: np.zeros(0)}
    J_all = {EDGE: np.zeros(0)}
    U_cum = qt.qeye(2)
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    out = []
    for kind, args in blocks:
        if kind == "rz":
            zeta = float(args["zeta"])
            rz_virtual_update(phase, 0, zeta)
        elif kind == "rx_pi2":
            t, Iq, Qq, Jq, _, _ = run_rx(np.pi / 2.0, args["duration_ns"], 0, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[0])
            t_blk = t if t_all.size == 0 else t_all[-1] + dt + t
            t_all = t_blk if t_all.size == 0 else np.concatenate([t_all, t_blk])
            for q in [0, 1]:
                I_all[q] = Iq[q] if I_all[q].size == 0 else np.concatenate([I_all[q], Iq[q]])
                Q_all[q] = Qq[q] if Q_all[q].size == 0 else np.concatenate([Q_all[q], Qq[q]])
            J_all[EDGE] = Jq[EDGE] if J_all[EDGE].size == 0 else np.concatenate([J_all[EDGE], Jq[EDGE]])
            U_cum = xy_unitary(np.pi / 2.0, phase[0]) * U_cum
        elif kind == "rx_mpi2":
            t, Iq, Qq, Jq, _, _ = run_rx(-np.pi / 2.0, args["duration_ns"], 0, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[0])
            t_blk = t if t_all.size == 0 else t_all[-1] + dt + t
            t_all = t_blk if t_all.size == 0 else np.concatenate([t_all, t_blk])
            for q in [0, 1]:
                I_all[q] = Iq[q] if I_all[q].size == 0 else np.concatenate([I_all[q], Iq[q]])
                Q_all[q] = Qq[q] if Q_all[q].size == 0 else np.concatenate([Q_all[q], Qq[q]])
            J_all[EDGE] = Jq[EDGE] if J_all[EDGE].size == 0 else np.concatenate([J_all[EDGE], Jq[EDGE]])
            U_cum = xy_unitary(-np.pi / 2.0, phase[0]) * U_cum
        elif kind == "rx_pi":
            t, Iq, Qq, Jq, _, _ = run_rx(np.pi, args["duration_ns"], 0, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[0])
            t_blk = t if t_all.size == 0 else t_all[-1] + dt + t
            t_all = t_blk if t_all.size == 0 else np.concatenate([t_all, t_blk])
            for q in [0, 1]:
                I_all[q] = Iq[q] if I_all[q].size == 0 else np.concatenate([I_all[q], Iq[q]])
                Q_all[q] = Qq[q] if Q_all[q].size == 0 else np.concatenate([Q_all[q], Qq[q]])
            J_all[EDGE] = Jq[EDGE] if J_all[EDGE].size == 0 else np.concatenate([J_all[EDGE], Jq[EDGE]])
            U_cum = xy_unitary(np.pi, phase[0]) * U_cum
        H = build_H(t_all, dt, I_all, Q_all, J_all, NQ, [EDGE]) if t_all.size else None
        rho_sim = rho0 if H is None else simulate(t_all, H, rho0, dt).states[-1]
        U_full = qt.tensor([U_cum, qt.qeye(2)])
        rho_id = U_full * rho0 * U_full.dag()
        out.append((kind, fidelity_dm(rho_sim, rho_id)))
    return out

def ideal_unitary_from_blocks(blocks):
    phase = 0.0
    U = qt.qeye(2)
    for kind, args in blocks:
        if kind == "rz":
            phase += float(args["zeta"])
        elif kind == "rx_pi2":
            U = xy_unitary(np.pi / 2.0, phase) * U
        elif kind == "rx_mpi2":
            U = xy_unitary(-np.pi / 2.0, phase) * U
        elif kind == "rx_pi":
            U = xy_unitary(np.pi, phase) * U
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
    return [("rz", {"q": int(q), "zeta": float(alpha)}),
            ("rx_pi2", {"q": int(q), "duration_ns": float(dur_ns)}),
            ("rz", {"q": int(q), "zeta": float(beta)}),
            ("rx_mpi2", {"q": int(q), "duration_ns": float(dur_ns)}),
            ("rz", {"q": int(q), "zeta": float(gamma)})]

def simulate_sequence_2q(blocks, dt_ns=0.25, sigma_frac=0.2):
    dt = dt_ns * 1e-9
    phase = {0: 0.0, 1: 0.0}
    t_all = np.array([], dtype=float)
    I_all = {0: np.zeros(0), 1: np.zeros(0)}
    Q_all = {0: np.zeros(0), 1: np.zeros(0)}
    J_all = {EDGE: np.zeros(0)}
    for kind, args in blocks:
        if kind == "rz":
            rz_virtual_update(phase, int(args["q"]), float(args["zeta"]))
        elif kind in ["rx_pi2", "rx_mpi2", "rx_pi"]:
            q = int(args["q"])
            ang = np.pi / 2.0 if kind == "rx_pi2" else (-np.pi / 2.0 if kind == "rx_mpi2" else np.pi)
            t, Iq, Qq, Jq, _, _ = run_rx(ang, float(args["duration_ns"]), q, dt_ns=dt_ns, sigma_frac=sigma_frac, phase=phase[q])
            t_blk = t if t_all.size == 0 else t_all[-1] + dt + t
            t_all = t_blk if t_all.size == 0 else np.concatenate([t_all, t_blk])
            for qq in [0, 1]:
                I_all[qq] = Iq[qq] if I_all[qq].size == 0 else np.concatenate([I_all[qq], Iq[qq]])
                Q_all[qq] = Qq[qq] if Q_all[qq].size == 0 else np.concatenate([Q_all[qq], Qq[qq]])
            J_all[EDGE] = Jq[EDGE] if J_all[EDGE].size == 0 else np.concatenate([J_all[EDGE], Jq[EDGE]])
        elif kind == "xy":
            theta = float(args["theta"])
            dur = float(args["duration_ns"])
            t, Iq, Qq, Jq, _, _ = run_xy_entangler(theta, dur, dt_ns=dt_ns)
            t_blk = t if t_all.size == 0 else t_all[-1] + dt + t
            t_all = t_blk if t_all.size == 0 else np.concatenate([t_all, t_blk])
            for qq in [0, 1]:
                I_all[qq] = Iq[qq] if I_all[qq].size == 0 else np.concatenate([I_all[qq], Iq[qq]])
                Q_all[qq] = Qq[qq] if Q_all[qq].size == 0 else np.concatenate([Q_all[qq], Qq[qq]])
            J_all[EDGE] = Jq[EDGE] if J_all[EDGE].size == 0 else np.concatenate([J_all[EDGE], Jq[EDGE]])
    H = build_H(t_all, dt, I_all, Q_all, J_all, NQ, [EDGE]) if t_all.size else None
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    rho_f = rho0 if H is None else simulate(t_all, H, rho0, dt).states[-1]
    if rho_f.isket:
        rho_f = qt.ket2dm(rho_f)
    return t_all, I_all, Q_all, J_all, rho_f

def ideal_unitary_from_blocks_2q(blocks):
    phase = {0: 0.0, 1: 0.0}
    U = qt.qeye([2, 2])
    for kind, args in blocks:
        if kind == "rz":
            q = int(args["q"])
            phase[q] += float(args["zeta"])
        elif kind in ["rx_pi2", "rx_mpi2", "rx_pi"]:
            q = int(args["q"])
            ang = np.pi / 2.0 if kind == "rx_pi2" else (-np.pi / 2.0 if kind == "rx_mpi2" else np.pi)
            Uloc = xy_unitary(ang, phase[q])
            U = U_local(q, Uloc) * U
        elif kind == "xy":
            theta = float(args["theta"])
            U = U_xy(theta) * U
    return U

def random_zyz():
    return np.random.uniform(-np.pi, np.pi), np.random.uniform(0.0, np.pi), np.random.uniform(-np.pi, np.pi)

def random_su4_blocks(depth, single_dur_ns=25.0, ent_dur_ns=40.0, ent_p=0.4):
    blocks = []
    for k in range(depth):
        if np.random.rand() < ent_p:
            theta = np.random.uniform(0.0, 0.5 * np.pi)
            blocks.append(("xy", {"theta": float(theta), "duration_ns": float(ent_dur_ns)}))
        else:
            q = int(np.random.choice([0, 1]))
            a, b, c = random_zyz()
            blocks.extend(synth_su2_blocks(q, a, b, c, single_dur_ns))
    return blocks

def test_random_su4(num_tests=5, depth=4, dt_ns=0.25, sigma_frac=0.2):
    e0, e1, *_ = single_qubit_ops()
    rho0 = qt.ket2dm(qt.tensor([e0, e0]))
    out = []
    for k in range(num_tests):
        blocks = random_su4_blocks(depth)
        _, _, _, _, rho_sim = simulate_sequence_2q(blocks, dt_ns=dt_ns, sigma_frac=sigma_frac)
        U_id = ideal_unitary_from_blocks_2q(blocks)
        rho_tar = U_id * rho0 * U_id.dag()
        F = fidelity_dm(rho_sim, rho_tar)
        out.append((blocks, F))
    return out

def main():
    t_rx_pi, I_rx_pi, Q_rx_pi, J_rx_pi, res_rx_pi, p_rx_pi = run_rx(np.pi, 25.0, 0, dt_ns=0.25, sigma_frac=0.2, phase=0.0)
    t_rx_pi2, I_rx_pi2, Q_rx_pi2, J_rx_pi2, res_rx_pi2, p_rx_pi2 = run_rx(np.pi / 2.0, 25.0, 1, dt_ns=0.25, sigma_frac=0.2, phase=0.0)
    t_sw, I_sw, Q_sw, J_sw, res_sw, p_sw = run_iswap(50.0, dt_ns=0.25)
    theta_ry = np.pi / 2.0
    t_ry, I_ry, Q_ry, J_ry, res_ry, p_ry = run_rx(theta_ry, 25.0, 0, dt_ns=0.25, sigma_frac=0.2, phase=np.pi / 2.0)
    e0, e1, *_ = single_qubit_ops()
    rho_target_rx = qt.ket2dm(qt.tensor([e1, e0]))
    F_rx = fidelity_dm(res_rx_pi.states[-1], rho_target_rx)
    rx_pi2_state_q1 = (e0 - 1j * e1).unit()
    rho_target_rx2 = qt.ket2dm(qt.tensor([e0, rx_pi2_state_q1]))
    F_rx2 = fidelity_dm(res_rx_pi2.states[-1], rho_target_rx2)
    psi10 = qt.tensor([e1, e0])
    psi01 = qt.tensor([e0, e1])
    bell_is = (psi10 - 1j * psi01).unit()
    rho_target_sw = qt.ket2dm(bell_is)
    F_sw = fidelity_dm(res_sw.states[-1], rho_target_sw)
    U_ry = ry_unitary(theta_ry)
    rho_target_ry = qt.ket2dm(qt.tensor([U_ry * e0, e0]))
    F_ry = fidelity_dm(res_ry.states[-1], rho_target_ry)
    alpha = np.random.uniform(-np.pi, np.pi)
    beta = np.random.uniform(0.0, np.pi)
    gamma = np.random.uniform(-np.pi, np.pi)
    a2, b2, c2 = float(alpha), float(beta), float(gamma)
    blocks = [("rz", {"zeta": a2}), ("rx_pi2", {"duration_ns": 25.0}), ("rz", {"zeta": b2}), ("rx_mpi2", {"duration_ns": 25.0}), ("rz", {"zeta": c2})]
    t_seq, I_seq, Q_seq, J_seq, rho_su2 = simulate_sequence_on_q0(blocks, dt_ns=0.25, sigma_frac=0.2)
    rho0q = qt.ket2dm(qt.tensor([e0, e0]))
    U_blocks = ideal_unitary_from_blocks(blocks)
    U_full_seq = qt.tensor([U_blocks, qt.qeye(2)])
    rho_seq_target = U_full_seq * rho0q * U_full_seq.dag()
    F_su2_frame_virtual = fidelity_dm(rho_su2, rho_seq_target)
    U_ZYZ = su2_from_zyz(alpha, beta, gamma)
    U_full_ZYZ = qt.tensor([U_ZYZ, qt.qeye(2)])
    rho_ZYZ = U_full_ZYZ * rho0q * U_full_ZYZ.dag()
    F_su2_physical_ZYZ = fidelity_dm(rho_su2, rho_ZYZ)
    rz_rule = {"type": "virtual_phase_advance", "rule": "I=Omega*cos(phi), Q=Omega*sin(phi); Rz(zeta): phi:=phi+zeta", "calibration_reference_phase_rad": 0.0}
    print_params_table(p_rx_pi, p_rx_pi2, rz_rule, p_sw)
    print("Fidelity rx_pi on q0 to |1>:", F_rx)
    print("Fidelity rx_pi_over_2 on q1 to Rx(pi/2)|0>:", F_rx2)
    print("Fidelity sqrt(iSWAP) on |10> to (|10|-i|01>)/sqrt2:", F_sw)
    print("Fidelity ry_pi_over_2 on q0 to Ry(pi/2)|0>:", F_ry)
    step_F = step_fidelities_matched(blocks, dt_ns=0.25, sigma_frac=0.2)
    cum = []
    for k, Fk in step_F:
        cum.append(k)
        print("Fidelity after " + " + ".join(cum) + ":", Fk)
    print("Fidelity SU2 sequence against frame-virtual target:", F_su2_frame_virtual)
    print("Fidelity SU2 sequence against physical ZYZ target:", F_su2_physical_ZYZ)
    rnd_tests = test_random_su4(num_tests=5, depth=4, dt_ns=0.25, sigma_frac=0.2)
    for idx, (blk, F) in enumerate(rnd_tests, start=1):
        print("Random SU4 test", idx, "fidelity:", F)
    out = {
        "global": {"dt_ns": 0.25, "num_qubits": NQ, "edge_01": [int(EDGE[0]), int(EDGE[1])]},
        "basis_single_qubit": {
            "rx_pi_25ns": {
                "duration_ns": float(p_rx_pi["duration_ns"]),
                "sigma_ns": float(p_rx_pi["sigma_ns"]),
                "omega_peak_rad_per_s": float(p_rx_pi["omega_peak_rad_per_s"]),
                "area_omega_rad": float(p_rx_pi["area_omega_rad"]),
                "default_phase_rad": float(p_rx_pi["default_phase_rad"])
            },
            "rx_pi_over_2_25ns": {
                "duration_ns": float(p_rx_pi2["duration_ns"]),
                "sigma_ns": float(p_rx_pi2["sigma_ns"]),
                "omega_peak_rad_per_s": float(p_rx_pi2["omega_peak_rad_per_s"]),
                "area_omega_rad": float(p_rx_pi2["area_omega_rad"]),
                "default_phase_rad": float(p_rx_pi2["default_phase_rad"])
            },
            "ry_pi_over_2_25ns_q0": {
                "duration_ns": float(p_ry["duration_ns"]),
                "sigma_ns": float(p_ry["sigma_ns"]),
                "omega_peak_rad_per_s": float(p_ry["omega_peak_rad_per_s"]),
                "area_omega_rad": float(p_ry["area_omega_rad"]),
                "default_phase_rad": float(p_ry["default_phase_rad"])
            }
        },
        "virtual_z": rz_rule,
        "single_qubit": {
            "rx_pi_q0": p_rx_pi,
            "rx_pi_over_2_q1": p_rx_pi2,
            "ry_pi_over_2_q0": p_ry,
            "rz_virtual": {"type": "virtual_phase_advance"}
        },
        "two_qubit": {
            "iswap_01": p_sw
        }
    }
    with open(OUT_JSON_PATH, "w") as f:
        json.dump(out, f, indent=2)

if __name__ == "__main__":
    main()
