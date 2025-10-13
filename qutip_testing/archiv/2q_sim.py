import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

# ===================== Parameters (angular units) =====================
Phi0 = 2.067833848e-15

# Qubit 1
w1   = 3.4e9 * 2*np.pi
eta1 = 0.2e9 * 2*np.pi
EC1  = eta1

# Qubit 2 (tunable via Phi2)
EC2     = 0.2e9 * 2*np.pi
EJ2_sum = 9.0e9  * 2*np.pi
alpha2  = 0.0

# Coupler (tunable via Phic)
ECc     = 0.25e9 * 2*np.pi
EJc_sum = 22.0e9  * 2*np.pi
alphac  = 0.0

# Static capacitive couplings (ref values at Phic=0)
g12     = 2*np.pi*4.7e6
g1c_ref = 2*np.pi*92.3e6
g2c_ref = 2*np.pi*92.3e6

# ---- Bias points (constant during the run; no coupler envelope) ----
Phi2_over_Phi0 = -0.143566                   # choose your qubit-2 bias
Phic_over_Phi0 = -0.247879 * 1.3                  # choose your coupler bias
Phi2 = Phi2_over_Phi0 * Phi0
Phic = Phic_over_Phi0 * Phi0

# ---- Single-qubit drive knobs (constant complex amplitudes) ----
# Complex Rabi amplitudes: Omega_k = Omega_k_max * (I_k + i Q_k)
Omega1_max = 2*np.pi*10e6
Omega2_max = 2*np.pi*10e6

I1_amp = 0.4   # set your amplitudes
Q1_amp = 0.0
I2_amp = 0.0
Q2_amp = 0.0

# Drive detunings in the rotating frame used below (often 0)
delta1 = 0.0
delta2 = 0.0

# ---- Simulation window ----
t_end_ns = 200.0
t_points = 2501

# ===================== Helpers =====================
GHz = 2*np.pi*1e9
MHz = 2*np.pi*1e6

def EJ_of_flux(EJ_sum, alpha, Phi):
    # EJ_eff = EJ_sum * sqrt( cos^2(phi) + alpha^2 sin^2(phi) ), phi=πΦ/Φ0
    phi = np.pi * Phi / Phi0
    return EJ_sum * np.sqrt(np.cos(phi)**2 + (alpha**2)*(np.sin(phi)**2))

def w01_from_EJ_EC(EJ, EC):
    # transmon 0-1 angular frequency
    return np.sqrt(8.0*EJ*EC) - EC

def EJ_from_w01_EC(w01, EC):
    return (w01 + EC)**2 / (8.0*EC)

def xi_from_EJ_EC(EJ, EC):
    return 2.0*EC / max(EJ,1e-30)

def w2_of_flux(Phi2):
    EJ2 = EJ_of_flux(EJ2_sum, alpha2, Phi2)
    return w01_from_EJ_EC(EJ2, EC2)

def wc_of_flux(Phic):
    EJc = EJ_of_flux(EJc_sum, alphac, Phic)
    return w01_from_EJ_EC(EJc, ECc)

# g_jc(Φc) from Eq. (A16) normalized to g_ref at Φc=0
def gjc_via_A16_normalized(g_ref, EJj, ECj, Phic):
    EJc_Phi = EJ_of_flux(EJc_sum, alphac, Phic)
    EJc_0   = EJ_of_flux(EJc_sum, alphac, 0.0)
    xi_j     = xi_from_EJ_EC(EJj, ECj)
    xi_c_Phi = xi_from_EJ_EC(EJc_Phi, ECc)
    xi_c_0   = xi_from_EJ_EC(EJc_0,   ECc)
    f_Phi = (EJc_Phi/ECc)**0.25 * (1.0 - 0.125*(xi_c_Phi + xi_j))
    f_0   = (EJc_0  /ECc)**0.25 * (1.0 - 0.125*(xi_c_0   + xi_j))
    return g_ref * (f_Phi / max(f_0,1e-30))

def g01_g02_g20(Phi2, Phic, EJ1, EC1):
    w2 = w2_of_flux(Phi2)
    wc = wc_of_flux(Phic)
    EJ2 = EJ_of_flux(EJ2_sum, alpha2, Phi2)

    g1c = gjc_via_A16_normalized(g1c_ref, EJ1, EC1, Phic)
    g2c = gjc_via_A16_normalized(g2c_ref, EJ2, EC2, Phic)

    D1 = wc - w1
    D2 = wc - w2
    L1 = wc + w1
    L2 = wc + w2
    eta2 = EC2

    # Full A25 with counter-rotating contributions via (wc + wj) terms
    eps = 1e-30
    invs01 = (1.0/max(D1,eps) + 1.0/max(D2,eps) + 1.0/max(L1,eps) + 1.0/max(L2,eps))
    g01 = g12 - 0.5*g1c*g2c*invs01

    # Leakage channels (rotating only; good approx in our regime)
    g02 = (np.sqrt(2.0)*g12
           - 0.5*np.sqrt(2.0)*g1c*g2c*( 1.0/max(D1,eps) + 1.0/max(D2+eta2,eps) ))
    g20 = (np.sqrt(2.0)*g12
           - 0.5*np.sqrt(2.0)*g1c*g2c*( 1.0/max(D2,eps) + 1.0/max(D1+eta1,eps) ))

    return g01, g02, g20, w2, wc, g1c, g2c

# ===================== Build Hamiltonian =====================
# basis |0>,|1>,|2> per qubit
e0 = qt.basis(3,0); e1 = qt.basis(3,1); e2 = qt.basis(3,2)
I3 = qt.qeye(3)

# Two-qubit states used in exchange/leakage
ket01 = qt.tensor(e0,e1)
ket10 = qt.tensor(e1,e0)
ket11 = qt.tensor(e1,e1)
ket02 = qt.tensor(e0,e2)
ket20 = qt.tensor(e2,e0)

# Single-qubit ladder operators in the 3-level truncation
s01_1 = qt.tensor(e0*e1.dag(), I3)
s12_1 = qt.tensor(e1*e2.dag(), I3)
s01_2 = qt.tensor(I3, e0*e1.dag())
s12_2 = qt.tensor(I3, e1*e2.dag())

# Exchange/leakage projectors
P10_01 = ket10 * ket01.dag()
P11_02 = ket11 * ket02.dag()
P20_11 = ket20 * ket11.dag()

# Observables
P10 = ket10*ket10.dag()
P01 = ket01*ket01.dag()
P1_only_q1 = qt.tensor(e1*e1.dag(), I3)
P1_only_q2 = qt.tensor(I3, e1*e1.dag())

# Frequencies, couplings from the chosen DC biases
EJ1 = EJ_from_w01_EC(w1, EC1)
g01, g02, g20, w2, wc, g1c, g2c = g01_g02_g20(Phi2, Phic, EJ1, EC1)
Delta = w2 - w1
eta2  = EC2

# Constant complex Rabi amplitudes you set via I/Q knobs
Omega1 = Omega1_max * (I1_amp + 1j*Q1_amp)
Omega2 = Omega2_max * (I2_amp + 1j*Q2_amp)

# Time-dependent coefficient functions (rotating-frame phases only)
def c01(t, args):   return args["g01"] * np.exp(-1j * args["Delta"] * t)
def c01c(t, args):  return np.conj(args["g01"]) * np.exp(+1j * args["Delta"] * t)

def c02(t, args):   return args["g02"] * np.exp(-1j * (args["Delta"] - args["eta2"]) * t)
def c02c(t, args):  return np.conj(args["g02"]) * np.exp(+1j * (args["Delta"] - args["eta2"]) * t)

def c20(t, args):   return args["g20"] * np.exp(-1j * (args["Delta"] + args["eta1"]) * t)
def c20c(t, args):  return np.conj(args["g20"]) * np.exp(+1j * (args["Delta"] + args["eta1"]) * t)

# Single-qubit drives (constant amplitude, rotating at chosen detunings)
def cd01_1(t, args):  return 0.5*args["Omega1"] * np.exp(-1j * args["delta1"] * t)
def cd01_1c(t, args): return 0.5*np.conj(args["Omega1"]) * np.exp(+1j * args["delta1"] * t)
def cd12_1(t, args):  return 0.5*np.sqrt(2.0)*args["Omega1"] * np.exp(-1j * (args["delta1"] - args["eta1"]) * t)
def cd12_1c(t, args): return 0.5*np.sqrt(2.0)*np.conj(args["Omega1"]) * np.exp(+1j * (args["delta1"] - args["eta1"]) * t)

def cd01_2(t, args):  return 0.5*args["Omega2"] * np.exp(-1j * args["delta2"] * t)
def cd01_2c(t, args): return 0.5*np.conj(args["Omega2"]) * np.exp(+1j * args["delta2"] * t)
def cd12_2(t, args):  return 0.5*np.sqrt(2.0)*args["Omega2"] * np.exp(-1j * (args["delta2"] - args["eta2"]) * t)
def cd12_2c(t, args): return 0.5*np.sqrt(2.0)*np.conj(args["Omega2"]) * np.exp(+1j * (args["delta2"] - args["eta2"]) * t)

# Assemble H(t): exchange + leakage + single-qubit drives
H0 = 0 * P10
H = [H0,
     [P10_01, c01],   [P10_01.dag(), c01c],
     [P11_02, c02],   [P11_02.dag(), c02c],
     [P20_11, c20],   [P20_11.dag(), c20c],
     [s01_1,  cd01_1],[s01_1.dag(),  cd01_1c],
     [s12_1,  cd12_1],[s12_1.dag(),  cd12_1c],
     [s01_2,  cd01_2],[s01_2.dag(),  cd01_2c],
     [s12_2,  cd12_2],[s12_2.dag(),  cd12_2c]
    ]

args = dict(
    g01=g01, g02=g02, g20=g20,
    Delta=Delta, eta1=eta1, eta2=EC2,
    Omega1=Omega1, Omega2=Omega2,
    delta1=delta1, delta2=delta2
)

# ===================== Run a quick sim =====================
psi0 = ket10  # start in |10>
t_end = t_end_ns * 1e-9
tlist = np.linspace(0.0, t_end, t_points)
e_ops = [P10, P01, P1_only_q1, P1_only_q2]

print("=== DC biases & derived ===")
print(f"Phi2/Φ0 = {Phi2_over_Phi0:.6f}   Phic/Φ0 = {Phic_over_Phi0:.6f}")
print(f"w1, w2, wc [GHz]: {w1/GHz:.6f}, {w2/GHz:.6f}, {wc/GHz:.6f}")
print(f"Delta [MHz] = {(Delta)/MHz:.6f}")
print(f"g01,g02,g20 [MHz]: {g01/MHz:.6f}, {g02/MHz:.6f}, {g20/MHz:.6f}")
print(f"g1c,g2c [MHz]: {g1c/MHz:.3f}, {g2c/MHz:.3f}")
print("=== Single-qubit drives ===")
print(f"Omega1 = {Omega1/(2*np.pi*1e6):.3f} MHz   (I1={I1_amp}, Q1={Q1_amp}),  delta1/(2π)={delta1/(2*np.pi*1e6):.3f} MHz")
print(f"Omega2 = {Omega2/(2*np.pi*1e6):.3f} MHz   (I2={I2_amp}, Q2={Q2_amp}),  delta2/(2π)={delta2/(2*np.pi*1e6):.3f} MHz")

res = qt.mesolve(H, psi0, tlist, [], e_ops, args=args,
                 options={"nsteps": 300000, "rtol": 1e-9, "atol": 1e-9})

P10_t = res.expect[0]
P01_t = res.expect[1]
P1q1  = res.expect[2]
P1q2  = res.expect[3]

plt.figure(figsize=(7,4))
plt.plot(tlist*1e9, P10_t, label="P(|10⟩)")
plt.plot(tlist*1e9, P01_t, label="P(|01⟩)")
plt.plot(tlist*1e9, P1q1,  '--', label="P1 (q1)")
plt.plot(tlist*1e9, P1q2,  '--', label="P1 (q2)")
plt.xlabel("Time (ns)")
plt.ylabel("Population")
plt.title("Exchange + leakage + constant I/Q drives")
plt.legend()
plt.tight_layout()
plt.show()
