import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import json

# ----------------- constants -----------------
Phi0 = 2.067833848e-15
GHz = 2*np.pi*1e9
MHz = 2*np.pi*1e6

# ----------------- device parameters -----------------
w1   = 3.4e9 * 2*np.pi
eta1 = 0.2e9 * 2*np.pi
EC1  = eta1

EC2      = 0.2e9 * 2*np.pi
EJ2_sum  = 9.0e9  * 2*np.pi
alpha2   = 0.0

ECc      = 0.25e9 * 2*np.pi
EJc_sum  = 22.0e9 * 2*np.pi
alphac   = 0.0

g12      = 2*np.pi*4.7e6
g1c_ref  = 2*np.pi*92.3e6
g2c_ref  = 2*np.pi*92.3e6

# dissipation targets
T1_target = 25e-6
T2_target = 25e-6
gamma     = 1.0 / T1_target
gamma_phi = 1.0 / T2_target - 0.5 * gamma

# ----------------- helpers -----------------
def EJ_of_flux(EJ_sum, alpha, Phi):
    phi = np.pi * Phi / Phi0
    return EJ_sum * np.sqrt(np.cos(phi)**2 + (alpha**2)*(np.sin(phi)**2))

def w01_from_EJ_EC(EJ, EC):
    return np.sqrt(8.0*EJ*EC) - EC

def EJ_from_w01_EC(w01, EC):
    return (w01 + EC)**2 / (8.0*EC)

def xi_from_EJ_EC(EJ, EC):
    EJ_eff = max(EJ, 1e-30)
    return np.sqrt(2.0*EC / EJ_eff)

def w2_of_flux(Phi2):
    EJ2 = EJ_of_flux(EJ2_sum, alpha2, Phi2)
    return w01_from_EJ_EC(EJ2, EC2)

def wc_of_flux(Phic):
    EJc = EJ_of_flux(EJc_sum, alphac, Phic)
    return w01_from_EJ_EC(EJc, ECc)

# === Exact Eq. A16 ===
# g_jc = (E_jc / sqrt(2)) * ((EJ_j/EC_j) * (EJ_c(Φ)/EC_c))^(1/4) * [1 - (ξ_c + ξ_j)/8]
def gjc_A16(Ejc, EJj, ECj, EJcPhi):
    xi_j = xi_from_EJ_EC(EJj,   ECj)
    xi_c = xi_from_EJ_EC(EJcPhi, ECc)
    pref = Ejc / np.sqrt(2.0)
    scale = ((EJj/ECj) * (EJcPhi/ECc))**0.25
    corr  = (1.0 - 0.125*(xi_c + xi_j))
    return pref * scale * corr

# Back out geometric coupling energies E1c, E2c from the reference couplings at Φc=0
def backsolve_Ejc_from_ref(g_ref, EJj, ECj, EJc0):
    xi_j = xi_from_EJ_EC(EJj,  ECj)
    xi_c = xi_from_EJ_EC(EJc0, ECc)
    denom = ((EJj/ECj) * (EJc0/ECc))**0.25 * (1.0 - 0.125*(xi_c + xi_j))
    return g_ref * np.sqrt(2.0) / max(denom, 1e-30)

# Effective couplings
def couplings_static(Phi2, Phic, EJ1, EC1, eta1, eta2, E1c, E2c):
    EJ2 = EJ_of_flux(EJ2_sum, alpha2, Phi2)
    EJcPhi = EJ_of_flux(EJc_sum, alphac, Phic)
    w2 = w01_from_EJ_EC(EJ2, EC2)
    wc = w01_from_EJ_EC(EJcPhi, ECc)
    g1c = gjc_A16(E1c, EJ1, EC1, EJcPhi)
    g2c = gjc_A16(E2c, EJ2, EC2, EJcPhi)

    D1 = max(wc - w1, 1e-30)
    D2 = max(wc - w2, 1e-30)
    g01 = g12 - 0.5*g1c*g2c*(1.0/D1 + 1.0/D2)
    g02 = np.sqrt(2.0)*g12 - 0.5*np.sqrt(2.0)*g1c*g2c*(1.0/D1 + 1.0/max(D2+eta2, 1e-30))
    g20 = np.sqrt(2.0)*g12 - 0.5*np.sqrt(2.0)*g1c*g2c*(1.0/D2 + 1.0/max(D1+eta1, 1e-30))
    return g01, g02, g20, g1c, g2c, w2, wc

# Solvers
def solve_Phi2_on():
    def f(Phi2): return w2_of_flux(Phi2) - w1
    lo = -0.49*Phi0; hi = 0.49*Phi0
    flo, fhi = f(lo), f(hi)
    if flo*fhi > 0:
        xs = np.linspace(-0.49, 0.49, 2001)*Phi0
        vals = [f(x) for x in xs]
        sgn = np.sign(vals)
        idx = np.where(sgn[1:]*sgn[:-1] < 0)[0]
        if len(idx)==0: raise RuntimeError("could not bracket root for Phi2_on")
        lo, hi = xs[idx[0]], xs[idx[0]+1]
    a, b = lo, hi
    for _ in range(100):
        m = 0.5*(a+b)
        if f(a)*f(m) <= 0: b = m
        else: a = m
    return 0.5*(a+b)

def solve_Phic_zero_g01(Phi2_fix, EJ1, EC1, eta1, eta2, E1c, E2c):
    def g(Phi_c): return couplings_static(Phi2_fix, Phi_c, EJ1, EC1, eta1, eta2, E1c, E2c)[0]
    grid = np.linspace(-0.45, 0.45, 2001)*Phi0
    gvals = np.array([g(x) for x in grid])
    s = np.sign(gvals)
    idx = np.where(s[1:]*s[:-1] < 0)[0]
    if len(idx)==0:
        j = np.argmin(np.abs(gvals))
        return grid[j], False
    a, b = grid[idx[0]], grid[idx[0]+1]
    for _ in range(100):
        m = 0.5*(a+b)
        if g(a)*g(m) <= 0: b = m
        else: a = m
    return 0.5*(a+b), True

# iSWAP, DRAG, linearization, c_ops, and T1/T2 routines (unchanged except for passing E1c,E2c)
def simulate_iswap(Phi2, Phic, EJ1, EC1, eta1, eta2, E1c, E2c, c_ops, t_end=None, t_points=2001):
    g01, g02, g20, g1c, g2c, w2, wc = couplings_static(Phi2, Phic, EJ1, EC1, eta1, eta2, E1c, E2c)
    Delta = w2 - w1
    e0 = qt.basis(3,0); e1 = qt.basis(3,1)
    ket01 = qt.tensor(e0, e1); ket10 = qt.tensor(e1, e0)
    P10_01 = ket10 * ket01.dag()
    def c01(t):  return g01 * np.exp(-1j * Delta * t)
    def c01c(t): return np.conj(g01 * np.exp(-1j * Delta * t))
    H = [0*P10_01, [P10_01, c01], [P10_01.dag(), c01c]]
    psi0 = ket10
    if t_end is None:
        gabs = max(abs(g01), 1e-12)
        t_end = 1.2*np.pi/(2*gabs)
        t_end = np.clip(t_end, 20e-9, 300e-9)
    tlist = np.linspace(0.0, t_end, t_points)
    e_ops = [ket10*ket10.dag(), ket01*ket01.dag()]
    res = qt.mesolve(H, psi0, tlist, c_ops, e_ops, options={"nsteps": 300000, "rtol": 1e-9, "atol": 1e-9})
    info = {"g01": g01, "g02": g02, "g20": g20, "Delta": Delta, "w2": w2, "wc": wc, "g1c": g1c, "g2c": g2c, "t_end": t_end}
    return tlist, res.expect, info

# --- replace the whole simulate_single_qubit_drag with this single-transmon version ---
def simulate_single_qubit_drag(eta, alpha_drag=1.0, beta=1.0, sigma=10e-9, n_sigma=6, A_scale=1.0):
    e0 = qt.basis(3,0); e1 = qt.basis(3,1); e2 = qt.basis(3,2)
    s01 = e0*e1.dag(); s12 = e1*e2.dag()

    # DRAG Gaussian envelopes
    t_max = n_sigma*sigma
    tlist = np.linspace(0.0, t_max, 2001)
    A = A_scale * np.pi/(np.sqrt(2.0*np.pi)*sigma)
    t0 = 0.5*t_max
    def Ienv(t):
        x = (t - t0)
        return A*np.exp(-0.5*(x/sigma)**2)
    def Qenv(t):
        dI = -(t - t0)/(sigma**2) * Ienv(t)
        return -alpha_drag/eta * dI
    def Om(t):   return beta*(Ienv(t) + 1j*Qenv(t))
    def Om12(t): return np.sqrt(2.0)*Om(t)*np.exp(1j*eta*t)

    # Single-transmon Hamiltonian with DRAG
    H = [0*s01,
         [s01,        lambda t: 0.5*Om(t)],
         [s01.dag(),  lambda t: 0.5*np.conj(Om(t))],
         [s12,        lambda t: 0.5*Om12(t)],
         [s12.dag(),  lambda t: 0.5*np.conj(Om12(t))]]

    # Single-transmon collapse operators to realize T1, T2
    b = e0*e1.dag() + np.sqrt(2.0)*e1*e2.dag()
    n = b.dag()*b
    c_ops_single = [np.sqrt(gamma)*b, np.sqrt(2.0*gamma_phi)*n]

    rho0 = e0*e0.dag()
    e_ops = [e0*e0.dag(), e1*e1.dag(), e2*e2.dag()]
    res = qt.mesolve(H, rho0, tlist, c_ops_single, e_ops,
                     options={"nsteps": 300000, "rtol": 1e-9, "atol": 1e-9})
    return tlist, res.expect, {"A_base": A, "sigma": sigma, "t0": t0, "alpha_drag": alpha_drag, "beta": beta}

def linearize_g_about(Phi2_fix, Phic0, EJ1, EC1, eta1, eta2, E1c, E2c, dPhi=1e-3):
    d = dPhi*Phi0
    g01_m, _, g20_m, _, _, _, _ = couplings_static(Phi2_fix, Phic0 - d, EJ1, EC1, eta1, eta2, E1c, E2c)
    g01_0, _, g20_0, _, _, _, _ = couplings_static(Phi2_fix, Phic0,       EJ1, EC1, eta1, eta2, E1c, E2c)
    g01_p, _, g20_p, _, _, _, _ = couplings_static(Phi2_fix, Phic0 + d, EJ1, EC1, eta1, eta2, E1c, E2c)
    b01 = (g01_p - g01_m)/(2.0*d)
    b20 = (g20_p - g20_m)/(2.0*d)
    a01 = g01_0; a20 = g20_0
    return a01, b01, a20, b20

def build_c_ops(gamma, gamma_phi):
    e0 = qt.basis(3,0); e1 = qt.basis(3,1); e2 = qt.basis(3,2)
    b = e0*e1.dag() + np.sqrt(2.0)*e1*e2.dag()
    n = b.dag()*b
    I = qt.qeye(3)
    b1 = qt.tensor(b, I); b2 = qt.tensor(I, b)
    n1 = qt.tensor(n, I); n2 = qt.tensor(I, n)
    return [np.sqrt(gamma)*b1, np.sqrt(gamma)*b2,
            np.sqrt(2.0*gamma_phi)*n1, np.sqrt(2.0*gamma_phi)*n2]

def simulate_T1_T2(c_ops, t_T1=150e-6, t_T2=150e-6, npts=4001):
    e0 = qt.basis(3,0); e1 = qt.basis(3,1); I = qt.qeye(3)
    ket10 = qt.tensor(e1, e0); P10 = ket10*ket10.dag()
    H0 = 0*P10
    t1 = np.linspace(0.0, t_T1, npts)
    res1 = qt.mesolve(H0, ket10, t1, c_ops, [P10], options={"nsteps": 200000, "rtol": 1e-9, "atol": 1e-9})
    sx01 = e0*e1.dag() + e1*e0.dag(); Sx1 = qt.tensor(sx01, I)
    psi_plus = (qt.tensor((e0+e1).unit(), e0)).unit()
    t2 = np.linspace(0.0, t_T2, npts)
    res2 = qt.mesolve(H0, psi_plus, t2, c_ops, [Sx1], options={"nsteps": 200000, "rtol": 1e-9, "atol": 1e-9})
    return t1, res1.expect[0], t2, res2.expect[0]

# ----------------- calibration of Ejc from refs at Φc=0 -----------------
EJc0 = EJ_of_flux(EJc_sum, alphac, 0.0)
EJ2_0 = EJ_of_flux(EJ2_sum, alpha2, 0.0)
EJ1   = EJ_from_w01_EC(w1, EC1)

E1c = backsolve_Ejc_from_ref(g1c_ref, EJ1,  EC1,  EJc0)
E2c = backsolve_Ejc_from_ref(g2c_ref, EJ2_0, EC2, EJc0)

# ----------------- bias solves -----------------
Phi2_off = 0.0
Phic_off, ok_off = solve_Phic_zero_g01(Phi2_off, EJ1, EC1, eta1, EC2, E1c, E2c)

Phi2_on = solve_Phi2_on()
Phic_on_cancel, ok_on = solve_Phic_zero_g01(Phi2_on, EJ1, EC1, eta1, EC2, E1c, E2c)

delta_Phic = Phic_on_cancel - Phic_off

# ----------------- compute/linearize -----------------
g_off, g02_off, g20_off, g1c_off, g2c_off, w2_off, wc_off = couplings_static(Phi2_off, Phic_off, EJ1, EC1, eta1, EC2, E1c, E2c)
g_on_at_cancel, g02_on, g20_on_at_cancel, g1c_on, g2c_on, w2_on, wc_on = couplings_static(Phi2_on, Phic_on_cancel, EJ1, EC1, eta1, EC2, E1c, E2c)
a01, b01, a20, b20 = linearize_g_about(Phi2_on, Phic_on_cancel, EJ1, EC1, eta1, EC2, E1c, E2c, dPhi=1e-3)

# ----------------- dynamics (with dissipation) -----------------
c_ops = build_c_ops(gamma, gamma_phi)

t_drag_pi,  pop_drag_pi,  info_pi  = simulate_single_qubit_drag(eta1, alpha_drag=1.0, beta=1.0, sigma=10e-9, n_sigma=6, A_scale=1.0)
t_drag_pi2, pop_drag_pi2, info_pi2 = simulate_single_qubit_drag(eta1, alpha_drag=1.0, beta=1.0, sigma=10e-9, n_sigma=6, A_scale=0.5)

t_on,  ex_on,  info_on  = simulate_iswap(Phi2_on, 0.0, EJ1, EC1, eta1, EC2, E1c, E2c, c_ops, t_end=None)
t_end_on = info_on["t_end"]
t_off, ex_off, info_off = simulate_iswap(Phi2_off, Phic_off, EJ1, EC1, eta1, EC2, E1c, E2c, c_ops, t_end=t_end_on)

t1, P10_t, t2, Sx_t = simulate_T1_T2(c_ops, t_T1=150e-6, t_T2=150e-6, npts=3001)

# ----------------- JSON dump -----------------
params = {
    "biases": {
        "phi2_off_over_phi0_abs": float(abs(Phi2_off/Phi0)),
        "phic_off_over_phi0_abs": float(abs(Phic_off/Phi0)),
        "phi2_on_over_phi0_abs": float(abs(Phi2_on/Phi0)),
        "phic_on_cancel_over_phi0_abs": float(abs(Phic_on_cancel/Phi0)),
        "delta_phic_on_minus_off_over_phi0_abs": float(abs(delta_Phic/Phi0))
    },
    "anharmonicity": {"eta1_rad_per_s": float(eta1), "eta2_rad_per_s": float(EC2)},
    "linearized_couplings_at_on_cancel": {
        "a01_rad_per_s": float(a01), "b01_rad_per_s_per_Wb": float(b01),
        "a20_rad_per_s": float(a20), "b20_rad_per_s_per_Wb": float(b20),
        "reference_phic_on_cancel_Wb": float(Phic_on_cancel)
    },
    "single_qubit_drive": {
        "alpha_drag": 1.0, "beta": 1.0, "sigma_s": 10e-9,
        "A_gaussian_I_peak_pi": float(info_pi["A_base"]),
        "A_gaussian_I_peak_pi_over_2": float(0.5*info_pi["A_base"])
    },
    "dissipation": {
        "T1_s": float(T1_target), "T2_s": float(T2_target),
        "gamma_relax_sinv": float(gamma), "gamma_phi_sinv": float(gamma_phi)
    },
    "A16_parameters": {
        "E1c_rad_per_s": float(E1c), "E2c_rad_per_s": float(E2c),
        "EJc0_rad_per_s": float(EJc0), "EJ2_0_rad_per_s": float(EJ2_0), "EJ1_rad_per_s": float(EJ1)
    }
}
with open("calibration_params.json", "w") as f:
    json.dump(params, f, indent=2)

# ----------------- plots -----------------
Phic_grid = np.linspace(Phic_on_cancel-0.05*Phi0, Phic_on_cancel+0.05*Phi0, 801)
g01_true = []; g01_lin=[]
for ph in Phic_grid:
    g01_val, _, _, _, _, _, _ = couplings_static(Phi2_on, ph, EJ1, EC1, eta1, EC2, E1c, E2c)
    g01_true.append(g01_val/MHz)
    g01_lin.append((a01 + b01*(ph - Phic_on_cancel))/MHz)
g01_true = np.array(g01_true); g01_lin = np.array(g01_lin)

plt.figure(figsize=(7,4))
plt.plot(Phic_grid/Phi0, g01_true, label="g01 true [MHz]")
plt.plot(Phic_grid/Phi0, g01_lin, linestyle="--", label="g01 linearized [MHz]")
plt.axhline(0.0, color="k", lw=0.8)
plt.axvline(Phic_on_cancel/Phi0, color="C1", lw=1.0, linestyle="--", label="Phic on cancel")
plt.xlabel("Phic / Phi0"); plt.ylabel("g01 [MHz]")
plt.title("g01 near ON cancel: true vs linearized"); plt.legend(); plt.tight_layout()

plt.figure(figsize=(7,4))
plt.plot(t_on*1e9, ex_on[0], label="P10 ON"); plt.plot(t_on*1e9, ex_on[1], label="P01 ON")
plt.xlabel("Time [ns]"); plt.ylabel("Population"); plt.title("iSWAP ON with dissipation")
plt.legend(); plt.tight_layout()

plt.figure(figsize=(7,4))
plt.plot(t_off*1e9, ex_off[0], label="P10 OFF"); plt.plot(t_off*1e9, ex_off[1], label="P01 OFF")
plt.xlabel("Time [ns]"); plt.ylabel("Population"); plt.title("iSWAP OFF with dissipation")
plt.legend(); plt.tight_layout()

plt.figure(figsize=(7,4))
plt.plot(t_drag_pi*1e9,  pop_drag_pi[0], label="P0, pi")
plt.plot(t_drag_pi*1e9,  pop_drag_pi[1], label="P1, pi")
plt.plot(t_drag_pi*1e9,  pop_drag_pi[2], label="P2, pi")
plt.plot(t_drag_pi2*1e9, pop_drag_pi2[1], linestyle="--", label="P1, pi/2")
plt.xlabel("Time [ns]"); plt.ylabel("Population"); plt.title("Single-qubit Gaussian DRAG with dissipation")
plt.legend(); plt.tight_layout()

plt.figure(figsize=(7,4))
plt.plot(t1*1e6, P10_t, label="P10(t)"); plt.plot(t1*1e6, np.exp(-t1/T1_target), linestyle="--", label="exp(-t/T1)")
plt.xlabel("Time [us]"); plt.ylabel("Population"); plt.title("T1 relaxation curve"); plt.legend(); plt.tight_layout()

plt.figure(figsize=(7,4))
plt.plot(t2*1e6, Sx_t, label="<sigma_x>"); plt.plot(t2*1e6, np.exp(-t2/T2_target), linestyle="--", label="exp(-t/T2)")
plt.xlabel("Time [us]"); plt.ylabel("Signal"); plt.title("T2 dephasing curve"); plt.legend(); plt.tight_layout()

plt.show()
print("Saved calibration_params.json")
