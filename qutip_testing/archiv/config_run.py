import json
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

Phi0 = 2.067833848e-15
GHz = 2*np.pi*1e9
MHz = 2*np.pi*1e6

def load_calibration(path="calibration_params.json"):
    with open(path, "r") as f:
        cfg = json.load(f)
    dphi_abs = cfg["biases"]["delta_phic_on_minus_off_over_phi0_abs"]
    b01 = cfg["linearized_couplings_at_on_cancel"]["b01_rad_per_s_per_Wb"]
    b20 = cfg["linearized_couplings_at_on_cancel"]["b20_rad_per_s_per_Wb"]
    A_pi = cfg["single_qubit_drive"]["A_gaussian_I_peak_pi"]
    eta1 = cfg["anharmonicity"]["eta1_rad_per_s"]
    eta2 = cfg["anharmonicity"]["eta2_rad_per_s"]
    T1 = cfg["dissipation"]["T1_s"]
    T2 = cfg["dissipation"]["T2_s"]
    return dphi_abs, b01, b20, A_pi, eta1, eta2, T1, T2

def load_envelopes(path="pulses_envelopes.json"):
    with open(path, "r") as f:
        p = json.load(f)
    s = np.asarray(p["s_envelope"], dtype=float)
    I1 = np.asarray(p["I1_envelope"], dtype=float)
    Q1 = np.asarray(p["Q1_envelope"], dtype=float)
    I2 = np.asarray(p["I2_envelope"], dtype=float)
    Q2 = np.asarray(p["Q2_envelope"], dtype=float)
    n = len(s)
    if not (len(I1)==n and len(Q1)==n and len(I2)==n and len(Q2)==n):
        raise ValueError("all envelopes must have the same length")
    t = np.arange(n, dtype=float)*0.5e-9
    return t, s, I1, Q1, I2, Q2

def qutrit_ops():
    e0 = qt.basis(3,0)
    e1 = qt.basis(3,1)
    e2 = qt.basis(3,2)
    s01 = e0*e1.dag()
    s10 = s01.dag()
    s12 = e1*e2.dag()
    s21 = s12.dag()
    b = s01 + np.sqrt(2.0)*s12
    n = b.dag()*b
    P00 = qt.tensor(e0,e0)*qt.tensor(e0,e0).dag()
    P01 = qt.tensor(e0,e1)*qt.tensor(e0,e1).dag()
    P10 = qt.tensor(e1,e0)*qt.tensor(e1,e0).dag()
    P11 = qt.tensor(e1,e1)*qt.tensor(e1,e1).dag()
    P02 = qt.tensor(e0,e2)*qt.tensor(e0,e2).dag()
    P20 = qt.tensor(e2,e0)*qt.tensor(e2,e0).dag()
    return e0,e1,e2,s01,s10,s12,s21,b,n,P00,P01,P10,P11,P02,P20

def make_coeff(tgrid, ygrid):
    ygrid = np.asarray(ygrid, dtype=float)
    def c(t):
        return float(np.interp(t, tgrid, ygrid, left=ygrid[0], right=ygrid[-1]))
    return c

def build_c_ops_qutrit(T1, T2):
    gamma = 1.0/float(T1)
    gamma_phi = max(1.0/float(T2) - 0.5*gamma, 0.0)
    e0,e1,e2,s01,s10,s12,s21,b,n,P00,P01,P10,P11,P02,P20 = qutrit_ops()
    I3 = qt.qeye(3)
    b1 = qt.tensor(b, I3)
    b2 = qt.tensor(I3, b)
    n1 = qt.tensor(n, I3)
    n2 = qt.tensor(I3, n)
    c_ops = []
    if gamma > 0.0:
        c_ops += [np.sqrt(gamma)*b1, np.sqrt(gamma)*b2]
    if gamma_phi > 0.0:
        c_ops += [np.sqrt(2.0*gamma_phi)*n1, np.sqrt(2.0*gamma_phi)*n2]
    return c_ops

def build_H_time_dependent(tgrid, s_env, I1_env, Q1_env, I2_env, Q2_env,
                            dphi_abs, b01, b20, A_pi, eta1, eta2):
    e0,e1,e2,s01,s10,s12,s21,b,n,P00,P01,P10,P11,P02,P20 = qutrit_ops()
    I3 = qt.qeye(3)
    S22_1 = qt.tensor(e2*e2.dag(), I3)
    S22_2 = qt.tensor(I3, e2*e2.dag())
    H0 = -eta1*S22_1 - eta2*S22_2
    S01_1 = qt.tensor(s01, I3)
    S10_1 = S01_1.dag()
    S01_2 = qt.tensor(I3, s01)
    S10_2 = S01_2.dag()
    X01_1 = S01_1 + S10_1
    Y01_1 = -1j*S01_1 + 1j*S10_1
    X01_2 = S01_2 + S10_2
    Y01_2 = -1j*S01_2 + 1j*S10_2
    X01 = qt.tensor(e1,e0)*qt.tensor(e0,e1).dag()
    X02 = qt.tensor(e1,e1)*qt.tensor(e0,e2).dag()
    X20 = qt.tensor(e1,e1)*qt.tensor(e2,e0).dag()
    J01 = b01*(Phi0*dphi_abs)*np.asarray(s_env, dtype=float)
    J20 = b20*(Phi0*dphi_abs)*np.asarray(s_env, dtype=float)
    J02 = J20
    Om1x = 0.5*A_pi*np.asarray(I1_env, dtype=float)
    Om1y = 0.5*A_pi*np.asarray(Q1_env, dtype=float)
    Om2x = 0.5*A_pi*np.asarray(I2_env, dtype=float)
    Om2y = 0.5*A_pi*np.asarray(Q2_env, dtype=float)
    c_J01 = make_coeff(tgrid, J01)
    c_J02 = make_coeff(tgrid, J02)
    c_J20 = make_coeff(tgrid, J20)
    c_1x = make_coeff(tgrid, Om1x)
    c_1y = make_coeff(tgrid, Om1y)
    c_2x = make_coeff(tgrid, Om2x)
    c_2y = make_coeff(tgrid, Om2y)
    H = [H0,
         [X01 + X01.dag(), c_J01],
         [X02 + X02.dag(), c_J02],
         [X20 + X20.dag(), c_J20],
         [X01_1, c_1x], [Y01_1, c_1y],
         [X01_2, c_2x], [Y01_2, c_2y]]
    projs = (P00,P01,P10,P11,P02,P20)
    return H, projs, {"J01_max_MHz": float(np.max(np.abs(J01))/(2*np.pi*1e6)),
                      "J02_max_MHz": float(np.max(np.abs(J02))/(2*np.pi*1e6)),
                      "J20_max_MHz": float(np.max(np.abs(J20))/(2*np.pi*1e6)),
                      "Omega1_max_MHz": float(np.max(np.hypot(Om1x, Om1y))/(2*np.pi*1e6)),
                      "Omega2_max_MHz": float(np.max(np.hypot(Om2x, Om2y))/(2*np.pi*1e6))}

def run_from_json(calib_path="calibration_params.json", pulses_path="pulses_envelopes.json"):
    dphi_abs, b01, b20, A_pi, eta1, eta2, T1, T2 = load_calibration(calib_path)
    t, s_env, I1_env, Q1_env, I2_env, Q2_env = load_envelopes(pulses_path)
    H, projs, diag = build_H_time_dependent(t, s_env, I1_env, Q1_env, I2_env, Q2_env,
                                            dphi_abs, b01, b20, A_pi, eta1, eta2)
    c_ops = build_c_ops_qutrit(T1, T2)
    e1 = qt.basis(3,1)
    e0 = qt.basis(3,0)
    psi0 = qt.tensor(e1, e0)
    opts = qt.Options(method="bdf", rtol=1e-8, atol=1e-10, nsteps=300000, max_step=0.5e-9)
    res = qt.mesolve(H, psi0, t, c_ops=c_ops, e_ops=list(projs), options=opts)
    P00, P01, P10, P11, P02, P20 = res.expect
    t_ns = t*1e9
    print("J01_max [MHz]:", f"{diag['J01_max_MHz']:.3f}")
    print("J02_max [MHz]:", f"{diag['J02_max_MHz']:.3f}")
    print("J20_max [MHz]:", f"{diag['J20_max_MHz']:.3f}")
    print("Omega1_max [MHz]:", f"{diag['Omega1_max_MHz']:.3f}")
    print("Omega2_max [MHz]:", f"{diag['Omega2_max_MHz']:.3f}")
    plt.figure(figsize=(7,3))
    plt.plot(t_ns, s_env, label="s")
    plt.plot(t_ns, I1_env, label="I1")
    plt.plot(t_ns, Q1_env, label="Q1")
    plt.plot(t_ns, I2_env, label="I2")
    plt.plot(t_ns, Q2_env, label="Q2")
    plt.xlabel("Time [ns]")
    plt.ylabel("Envelope")
    plt.title("Control envelopes")
    plt.legend()
    plt.tight_layout()
    plt.figure(figsize=(7,4))
    plt.plot(t_ns, P10, label="|10|")
    plt.plot(t_ns, P01, label="|01|")
    plt.plot(t_ns, P11, label="|11|")
    plt.plot(t_ns, P02, label="|02|")
    plt.plot(t_ns, P20, label="|20|")
    plt.xlabel("Time [ns]")
    plt.ylabel("Population")
    plt.title("Two-qubit populations")
    plt.legend()
    plt.tight_layout()
    plt.show()

run_from_json("calibration_params.json", "pulses_envelopes.json")

