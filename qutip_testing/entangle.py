import json
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

PARAMS_JSON = "pulse_params_2q.json"
CALIB_JSON = "calibration_params.json"
NQ = 2
EDGE = (0, 1)

def load_calibration(path):
    with open(path, "r") as f:
        C = json.load(f)
    T1 = float(C["dissipation"]["T1_s"])
    T2 = float(C["dissipation"]["T2_s"])
    return {"T1": T1, "T2": T2}

def load_params(path):
    with open(path, "r") as f:
        P = json.load(f)
    p_pi = P["single_qubit"]["pi_q0"]
    p_sw = P["two_qubit"]["iswap_01"]
    dt = float(p_pi["dt_ns"]) * 1e-9
    return p_pi, p_sw, dt

def single_qubit_ops():
    e0 = qt.basis(2, 0)
    e1 = qt.basis(2, 1)
    s01 = e0 * e1.dag()
    s10 = s01.dag()
    X01 = s01 + s10
    Y01 = -1j * s01 + 1j * s10
    b = s01
    n = e1 * e1.dag()
    return e0, e1, s01, s10, X01, Y01, b, n

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

def build_c_ops(T1, T2, Nq):
    gamma = 1.0 / float(T1)
    gamma_phi = 1.0 / float(T2) - 0.5 * gamma
    if gamma_phi < 0.0:
        gamma_phi = 0.0
    _, _, _, _, _, _, b_single, n_single = single_qubit_ops()
    c_ops = []
    for i in range(Nq):
        b_i = lift_one(b_single, i, Nq)
        n_i = lift_one(n_single, i, Nq)
        if gamma > 0.0:
            c_ops.append(np.sqrt(gamma) * b_i)
        if gamma_phi > 0.0:
            c_ops.append(np.sqrt(2.0 * gamma_phi) * n_i)
    return c_ops

def gaussian_from_params(duration_ns, sigma_ns, omega_peak, dt):
    duration = duration_ns * 1e-9
    sigma = sigma_ns * 1e-9
    N = int(np.round(duration / dt)) + 1
    t = np.arange(N, dtype=float) * dt
    t0 = 0.5 * duration
    g = np.exp(-0.5 * ((t - t0) / sigma) ** 2)
    if g.max() > 0:
        g = omega_peak * g / g.max()
    return t, g

def append_block_all(t, I0, Q0, I1, Q1, J, dt, blk_I0, blk_Q0, blk_I1, blk_Q1, blk_J):
    if t.size == 0:
        t_blk = np.arange(blk_I0.size, dtype=float) * dt
    else:
        t_blk = t[-1] + dt + np.arange(blk_I0.size, dtype=float) * dt
    t_out = np.concatenate([t, t_blk]) if t.size else t_blk
    I0_out = np.concatenate([I0, blk_I0]) if I0.size else blk_I0
    Q0_out = np.concatenate([Q0, blk_Q0]) if Q0.size else blk_Q0
    I1_out = np.concatenate([I1, blk_I1]) if I1.size else blk_I1
    Q1_out = np.concatenate([Q1, blk_Q1]) if Q1.size else blk_Q1
    J_out = np.concatenate([J, blk_J]) if J.size else blk_J
    return t_out, I0_out, Q0_out, I1_out, Q1_out, J_out

def build_sequence(p_pi, p_sw, dt):
    phases = {0: 0.0, 1: 0.0}
    t = np.array([])
    I0 = np.array([])
    Q0 = np.array([])
    I1 = np.array([])
    Q1 = np.array([])
    J = np.array([])
    dur_pi_ns = float(p_pi["duration_ns"])
    sig_pi_ns = float(p_pi["sigma_ns"])
    A_pi = float(p_pi["omega_peak_rad_per_s"])
    t_pi, omega_pi = gaussian_from_params(dur_pi_ns, sig_pi_ns, A_pi, dt)
    blk_I0 = 1.0 * omega_pi * np.cos(phases[0])
    blk_Q0 = 1.0 * omega_pi * np.sin(phases[0])
    zeros = np.zeros_like(blk_I0)
    t, I0, Q0, I1, Q1, J = append_block_all(t, I0, Q0, I1, Q1, J, dt, blk_I0, blk_Q0, zeros, zeros, zeros)
    dur_sw_ns = float(p_sw["duration_ns"])
    J_amp = float(p_sw["J_amp_rad_per_s"])
    dur_sqrt_ns = 0.5 * dur_sw_ns
    N_sw = int(np.round((dur_sqrt_ns * 1e-9) / dt)) + 1
    blk_J = np.full(N_sw, J_amp, dtype=float)
    zeros_sw = np.zeros(N_sw, dtype=float)
    t, I0, Q0, I1, Q1, J = append_block_all(t, I0, Q0, I1, Q1, J, dt, zeros_sw, zeros_sw, zeros_sw, zeros_sw, blk_J)
    return t, I0, Q0, I1, Q1, J

def build_H(t, dt, I0, Q0, I1, Q1, J):
    _, _, s01, s10, X01, Y01, _, _ = single_qubit_ops()
    X0 = lift_one(X01, 0, NQ)
    Y0 = lift_one(Y01, 0, NQ)
    X1 = lift_one(X01, 1, NQ)
    Y1 = lift_one(Y01, 1, NQ)
    exch = lift_two(s10, 0, s01, 1, NQ) + lift_two(s01, 0, s10, 1, NQ)
    H = []
    if I0.size:
        H.append([X0, coeff_from_array(0.5 * I0, dt)])
    if Q0.size:
        H.append([Y0, coeff_from_array(0.5 * Q0, dt)])
    if I1.size:
        H.append([X1, coeff_from_array(0.5 * I1, dt)])
    if Q1.size:
        H.append([Y1, coeff_from_array(0.5 * Q1, dt)])
    if J.size:
        H.append([exch, coeff_from_array(J, dt)])
    return H

def proj_ops():
    e0, e1, _, _, _, _, _, _ = single_qubit_ops()
    p00 = qt.tensor([e0, e0]) * qt.tensor([e0, e0]).dag()
    p01 = qt.tensor([e0, e1]) * qt.tensor([e0, e1]).dag()
    p10 = qt.tensor([e1, e0]) * qt.tensor([e1, e0]).dag()
    p11 = qt.tensor([e1, e1]) * qt.tensor([e1, e1]).dag()
    return [p00, p01, p10, p11]

def dm_to_json(rho):
    M = rho.full()
    return {"basis": ["|00>", "|01>", "|10>", "|11>"], "rho_real": M.real.tolist(), "rho_imag": M.imag.tolist()}

def main():
    cal = load_calibration(CALIB_JSON)
    p_pi, p_sw, dt = load_params(PARAMS_JSON)
    t, I0, Q0, I1, Q1, J = build_sequence(p_pi, p_sw, dt)
    H = build_H(t, dt, I0, Q0, I1, Q1, J)
    c_ops = build_c_ops(cal["T1"], cal["T2"], NQ)
    e0, e1, _, _, _, _, _, _ = single_qubit_ops()
    psi0 = qt.tensor([e0, e0])
    opts = {"method": "bdf", "rtol": 1e-6, "atol": 1e-8, "nsteps": 200000, "max_step": dt, "progress_bar": False, "store_states": True}
    res = qt.mesolve(H, psi0, t, c_ops=c_ops, e_ops=proj_ops(), options=opts)
    rho_f = res.states[-1]
    if rho_f.isket:
        rho_f = qt.ket2dm(rho_f)
    with open("final_density_matrix.json", "w") as f:
        json.dump(dm_to_json(rho_f), f, indent=2)
    tt = t * 1e9
    P00, P01, P10, P11 = [np.asarray(sig) for sig in res.expect]
    plt.figure(figsize=(9, 4))
    plt.plot(tt, P00, label="P|00>")
    plt.plot(tt, P01, label="P|01>")
    plt.plot(tt, P10, label="P|10>")
    plt.plot(tt, P11, label="P|11>")
    plt.xlabel("Time [ns]")
    plt.ylabel("Population")
    plt.title("Computational-basis populations")
    plt.legend(ncol=4, fontsize=9)
    plt.tight_layout()
    plt.figure(figsize=(9, 4))
    plt.plot(tt, I0, label="I q0")
    plt.plot(tt, Q0, label="Q q0")
    plt.plot(tt, I1, label="I q1")
    plt.plot(tt, Q1, label="Q q1")
    plt.xlabel("Time [ns]")
    plt.ylabel("Drive")
    plt.title("XY pulses")
    plt.legend(ncol=4, fontsize=9)
    plt.tight_layout()
    plt.figure(figsize=(9, 3))
    plt.plot(tt, J, label="J 0-1")
    plt.xlabel("Time [ns]")
    plt.ylabel("J [rad/s]")
    plt.title("Exchange envelope")
    plt.legend()
    plt.tight_layout()
    print("Final time [ns]:", float(tt[-1]))
    print("Saved density matrix to final_density_matrix.json")
    plt.show()

if __name__ == "__main__":
    main()
