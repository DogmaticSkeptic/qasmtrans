import json
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

CALIB_PATH = "calibration_params.json"
PULSES_PATH = "pulses_envelopes_9q.json"
NQ = 9
NTRAJ = 64
SEED = 1234

def load_calibration(path):
    with open(path, "r") as f:
        C = json.load(f)
    eta = float(C["anharmonicity"].get("eta1_rad_per_s", C["anharmonicity"].get("eta_rad_per_s")))
    A_pi = float(C["single_qubit_drive"]["A_gaussian_I_peak_pi"])
    T1 = float(C["dissipation"]["T1_s"])
    T2 = float(C["dissipation"]["T2_s"])
    b01 = float(C["linearized_couplings_at_on_cancel"]["b01_rad_per_s_per_Wb"])
    b20 = float(C["linearized_couplings_at_on_cancel"]["b20_rad_per_s_per_Wb"])
    Phi0 = 2.067833848e-15
    dphi_over_phi0 = float(C["biases"]["delta_phic_on_minus_off_over_phi0_abs"])
    dphi_abs = dphi_over_phi0 * Phi0
    return {"eta": eta, "A_pi": A_pi, "T1": T1, "T2": T2, "b01": b01, "b20": b20, "Phi0": Phi0, "dphi_abs": dphi_abs}

def load_pulses(path):
    with open(path, "r") as f:
        p = json.load(f)
    dt = float(p["dt_s"])
    N = int(p["num_samples"])
    I = {}
    Q = {}
    for qrec in p["qubits"]:
        k = int(qrec["index"])
        I[k] = np.asarray(qrec["I_envelope"], dtype=float)
        Q[k] = np.asarray(qrec["Q_envelope"], dtype=float)
        if len(I[k]) != N or len(Q[k]) != N:
            raise ValueError("envelope length mismatch")
    S = {}
    edges = []
    for erec in p["edges"]:
        a, b = erec["pair"]
        key = tuple(sorted((int(a), int(b))))
        S[key] = np.asarray(erec["s_envelope"], dtype=float)
        if len(S[key]) != N:
            raise ValueError("s-envelope length mismatch")
        edges.append(key)
    t = np.arange(N, dtype=float) * dt
    return dt, t, I, Q, S, edges

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

def build_full_H(cal, dt, t, I_env, Q_env, S_env, edges, Nq):
    _, _, s01, s10, X01, Y01, _, _ = single_qubit_ops()
    X01_i = [lift_one(X01, i, Nq) for i in range(Nq)]
    Y01_i = [lift_one(Y01, i, Nq) for i in range(Nq)]
    pair_terms = {}
    for (a, b) in edges:
        i, j = a, b
        term01 = lift_two(s10, i, s01, j, Nq) + lift_two(s01, i, s10, j, Nq)
        pair_terms[(i, j)] = term01
    H = []
    A_pi = cal["A_pi"]
    for i in range(Nq):
        Ii = np.asarray(I_env.get(i, np.zeros_like(t)), dtype=float)
        Qi = np.asarray(Q_env.get(i, np.zeros_like(t)), dtype=float)
        Omx = 0.5 * A_pi * Ii
        Omy = 0.5 * A_pi * Qi
        if np.any(Omx != 0.0):
            H.append([X01_i[i], coeff_from_array(Omx, dt)])
        if np.any(Omy != 0.0):
            H.append([Y01_i[i], coeff_from_array(Omy, dt)])
    b01 = cal["b01"]
    Phi0 = cal["Phi0"]
    dphi_abs = cal["dphi_abs"]
    for (i, j) in edges:
        s_arr = np.asarray(S_env[(i, j)], dtype=float)
        J01 = b01 * (Phi0 * dphi_abs) * s_arr
        term01 = pair_terms[(i, j)]
        if np.any(J01 != 0.0):
            H.append([term01, coeff_from_array(J01, dt)])
    return H

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

def plot_pulses(t, I_env, Q_env, S_env, edges, Nq):
    tt = t * 1e9
    plt.figure(figsize=(9, 6))
    for i in range(Nq):
        if i in I_env:
            plt.plot(tt, I_env[i], label=f"I q{i}")
        if i in Q_env:
            plt.plot(tt, Q_env[i], label=f"Q q{i}")
    plt.xlabel("Time [ns]")
    plt.ylabel("Amplitude")
    plt.title("Single qubit I Q envelopes")
    plt.legend(ncol=3, fontsize=9)
    plt.tight_layout()

    plt.figure(figsize=(9, 6))
    for (i, j) in edges:
        key = tuple(sorted((i, j)))
        plt.plot(tt, S_env[key], label=f"s ({i},{j})")
    plt.xlabel("Time [ns]")
    plt.ylabel("Amplitude")
    plt.title("Coupler s envelopes per edge")
    plt.legend(ncol=3, fontsize=9)
    plt.tight_layout()

cal = load_calibration(CALIB_PATH)
dt, t, I_env, Q_env, S_env, edges = load_pulses(PULSES_PATH)
plot_pulses(t, I_env, Q_env, S_env, edges, NQ)
H = build_full_H(cal, dt, t, I_env, Q_env, S_env, edges, NQ)
c_ops = build_c_ops(cal["T1"], cal["T2"], NQ)
e0, e1, _, _, _, _, _, _ = single_qubit_ops()
psi0 = qt.tensor([e0] * NQ)
P1_loc = [lift_one(e1 * e1.dag(), i, NQ) for i in range(NQ)]
np.random.seed(SEED)
opts = {"method": "bdf", "rtol": 1e-5, "atol": 1e-7, "nsteps": 200000, "max_step": dt, "progress_bar": False}
res = qt.mesolve(H, psi0, t, c_ops=c_ops, e_ops=P1_loc, options=opts)
tt = t * 1e9
plt.figure(figsize=(8, 6))
for i, sig in enumerate(res.expect):
    plt.plot(tt, sig, label=f"P1(q{i})")
plt.xlabel("Time [ns]")
plt.ylabel("Population in |1>")
plt.title("9-qubit mesolve: P(|1>) per qubit")
plt.legend(ncol=3, fontsize=9)
plt.tight_layout()
plt.show()
