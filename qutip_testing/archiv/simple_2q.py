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
    a01 = cfg["linearized_couplings_at_on_cancel"]["a01_rad_per_s"]
    b01 = cfg["linearized_couplings_at_on_cancel"]["b01_rad_per_s_per_Wb"]
    a20 = cfg["linearized_couplings_at_on_cancel"]["a20_rad_per_s"]
    b20 = cfg["linearized_couplings_at_on_cancel"]["b20_rad_per_s_per_Wb"]
    A_pi = cfg["single_qubit_drive"]["A_gaussian_I_peak_pi"]
    T1 = cfg["dissipation"]["T1_s"]
    T2 = cfg["dissipation"]["T2_s"]
    eta1 = cfg["anharmonicity"]["eta1_rad_per_s"]
    eta2 = cfg["anharmonicity"]["eta2_rad_per_s"]
    return dphi_abs, a01, b01, a20, b20, A_pi, T1, T2, eta1, eta2

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
    if gamma > 0:
        c_ops += [np.sqrt(gamma)*b1, np.sqrt(gamma)*b2]
    if gamma_phi > 0:
        c_ops += [np.sqrt(2.0*gamma_phi)*n1, np.sqrt(2.0*gamma_phi)*n2]
    return c_ops

def g_from_s_time(b_slope, dphi_abs, s_t):
    return b_slope * (Phi0 * dphi_abs) * s_t

def make_interp(tlist, y):
    y = np.asarray(y, dtype=float)
    def coeff(t):
        return float(np.interp(t, tlist, y, left=y[0], right=y[-1]))
    return coeff

def build_H_qutrit_time_dependent(
    tlist,
    dphi_abs, b01, b20,
    eta1, eta2,
    s_t,
    A_pi,
    I1_t, Q1_t,
    I2_t, Q2_t,
    scale_g01=5.0, scale_g02=5.0, scale_g20=5.0
):
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
    s_t = np.asarray(s_t, dtype=float)
    I1_t = np.asarray(I1_t, dtype=float)
    Q1_t = np.asarray(Q1_t, dtype=float)
    I2_t = np.asarray(I2_t, dtype=float)
    Q2_t = np.asarray(Q2_t, dtype=float)
    g01_raw = g_from_s_time(b01, dphi_abs, s_t)
    g20_raw = g_from_s_time(b20, dphi_abs, s_t)
    g02_raw = g20_raw
    g01_t = scale_g01 * g01_raw
    g02_t = scale_g02 * g02_raw
    g20_t = scale_g20 * g20_raw
    Om1x_t = 0.5*A_pi*I1_t
    Om1y_t = 0.5*A_pi*Q1_t
    Om2x_t = 0.5*A_pi*I2_t
    Om2y_t = 0.5*A_pi*Q2_t
    c_g01 = make_interp(tlist, g01_t)
    c_g02 = make_interp(tlist, g02_t)
    c_g20 = make_interp(tlist, g20_t)
    c_1x  = make_interp(tlist, Om1x_t)
    c_1y  = make_interp(tlist, Om1y_t)
    c_2x  = make_interp(tlist, Om2x_t)
    c_2y  = make_interp(tlist, Om2y_t)
    H_td = [H0,
            [X01 + X01.dag(), c_g01],
            [X02 + X02.dag(), c_g02],
            [X20 + X20.dag(), c_g20],
            [X01_1, c_1x], [Y01_1, c_1y],
            [X01_2, c_2x], [Y01_2, c_2y]]
    projs = (P00,P01,P10,P11,P02,P20)
    diag = {
        "g01_max_MHz": float(np.max(np.abs(g01_t))/(2*np.pi*1e6)),
        "g02_max_MHz": float(np.max(np.abs(g02_t))/(2*np.pi*1e6)),
        "g20_max_MHz": float(np.max(np.abs(g20_t))/(2*np.pi*1e6)),
        "Omega1_max_MHz": float(np.max(np.hypot(Om1x_t, Om1y_t))/(2*np.pi*1e6)),
        "Omega2_max_MHz": float(np.max(np.hypot(Om2x_t, Om2y_t))/(2*np.pi*1e6)),
    }
    return H_td, projs, diag

def simulate_td(tlist, H_td, c_ops, projs, psi0=None):
    if psi0 is None:
        e1 = qt.basis(3,1)
        e0 = qt.basis(3,0)
        psi0 = qt.tensor(e1, e0)
    opts = qt.Options(method="bdf", rtol=1e-8, atol=1e-10, nsteps=300000)
    res = qt.mesolve(H_td, psi0, tlist, c_ops=c_ops, e_ops=list(projs), options=opts)
    return res

dphi_abs, a01, b01, a20, b20, A_pi, T1, T2, eta1, eta2 = load_calibration("calibration_params.json")
c_ops = build_c_ops_qutrit(T1, T2)

T = 200e-9
npts = 3001
tlist = np.linspace(0.0, T, npts)

s_amp = 1.0
t0 = 40e-9
tau = 80e-9
s_t = np.where((tlist >= t0) & (tlist <= t0 + tau), s_amp, 0.0)

I1_t = np.ones_like(tlist) * 0.5
Q1_t = np.zeros_like(tlist)
I2_t = np.zeros_like(tlist)
Q2_t = np.zeros_like(tlist)

H_td, projs, diag = build_H_qutrit_time_dependent(
    tlist, dphi_abs, b01, b20, eta1, eta2,
    s_t, A_pi, I1_t, Q1_t, I2_t, Q2_t,
    scale_g01=5.0, scale_g02=5.0, scale_g20=5.0
)

res = simulate_td(tlist, H_td, c_ops, projs)
P00, P01, P10, P11, P02, P20 = res.expect
tns = tlist*1e9

print("g01_max_MHz:", f"{diag['g01_max_MHz']:.3f}")
print("g02_max_MHz:", f"{diag['g02_max_MHz']:.3f}")
print("g20_max_MHz:", f"{diag['g20_max_MHz']:.3f}")
print("Omega1_max_MHz:", f"{diag['Omega1_max_MHz']:.3f}")

plt.figure(figsize=(7,3))
plt.plot(tns, s_t)
plt.xlabel("Time [ns]")
plt.ylabel("s(t)")
plt.tight_layout()

plt.figure(figsize=(7,4))
plt.plot(tns, P10, label="|10|")
plt.plot(tns, P01, label="|01|")
plt.plot(tns, P11, label="|11|")
plt.plot(tns, P02, label="|02|")
plt.plot(tns, P20, label="|20|")
plt.xlabel("Time [ns]")
plt.ylabel("Population")
plt.title("TD exchange 5x + 1Q drive unchanged")
plt.legend()
plt.tight_layout()

plt.show()
