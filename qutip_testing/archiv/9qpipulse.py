import json
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

CALIB_PATH = "calibration_params.json"
GRID_K     = 3
S_CONST    = 0.0
NPTS_PULSE = 2001         # fine sampling over the pulse window
NPTS_TAIL  = 600          # coarser sampling over the decay window
DECAY_SPAN = 2.0          # simulate to DECAY_SPAN * T1

with open(CALIB_PATH, "r") as f:
    C = json.load(f)

A_pi = float(C["single_qubit_drive"]["A_gaussian_I_peak_pi"])
T1   = float(C["dissipation"]["T1_s"])
T2   = float(C["dissipation"]["T2_s"])
b01  = float(C["linearized_couplings_at_on_cancel"]["b01_rad_per_s_per_Wb"])
Phi0 = 2.067833848e-15
dphi = float(C["biases"]["delta_phic_on_minus_off_over_phi0_abs"]) * Phi0

sx = qt.sigmax()
sy = qt.sigmay()
sz = qt.sigmaz()
sm = qt.sigmam()
I2 = qt.qeye(2)
P1 = (I2 - sz)/2.0

def lift_one(op, i, N):
    ops = [I2]*N
    ops[i] = op
    return qt.tensor(ops)

def lift_two(op_i, i, op_j, j, N):
    ops = [I2]*N
    ops[i] = op_i
    ops[j] = op_j
    return qt.tensor(ops)

def grid_order(k=3):
    return [r*k + c for r in range(k) for c in range(k)]

def edges_for_active(active_nodes, k=3):
    active = set(active_nodes)
    def idx(r,c): return r*k + c
    E = []
    for r in range(k):
        for c in range(k):
            u = idx(r,c)
            if u not in active: continue
            if c+1 < k:
                v = idx(r,c+1)
                if v in active and u < v: E.append((u, v))
            if r+1 < k:
                v = idx(r+1,c)
                if v in active and u < v: E.append((u, v))
    remap = {node:i for i,node in enumerate(active_nodes)}
    return [(remap[i], remap[j]) for (i,j) in E]

# system
active_nodes = grid_order(GRID_K)
N = len(active_nodes)
edges = edges_for_active(active_nodes, GRID_K)

# time grid with fine pulse region
t_pi    = np.pi / A_pi
t_final = DECAY_SPAN * T1
t_pulse = np.linspace(0.0, 5.0*t_pi, NPTS_PULSE)          # fully resolves the π pulse
t_tail  = np.linspace(5.0*t_pi, t_final, NPTS_TAIL, endpoint=True)
tlist   = np.unique(np.concatenate([t_pulse, t_tail]))
dt_min  = np.min(np.diff(tlist))

def u_drive(t):
    return 1.0 if t <= t_pi else 0.0

# Hamiltonian
H = []
for i in range(N):
    H.append([lift_one(sx, i, N), (lambda t, args=None, i=i: 0.5 * A_pi * u_drive(t))])

# optional exchange
J01 = b01 * (Phi0 * dphi) * float(S_CONST)
if abs(J01) > 0:
    H0 = 0
    for (ii, jj) in edges:
        XX = lift_two(sx, ii, sx, jj, N)
        YY = lift_two(sy, ii, sy, jj, N)
        H0 = H0 + 0.5 * J01 * (XX + YY)
    H.insert(0, H0)

# dissipation
gamma = 1.0 / float(T1)
gamma_phi = max(0.0, 1.0/float(T2) - 0.5*gamma)

c_ops = []
for i in range(N):
    c_ops.append(np.sqrt(gamma) * lift_one(sm, i, N))
    if gamma_phi > 0.0:
        c_ops.append(np.sqrt(gamma_phi/2.0) * lift_one(sz, i, N))

# initial state and observables
psi0 = qt.tensor([qt.basis(2,0)] * N)
e_ops = [lift_one(P1, i, N) for i in range(N)]

# solver options: restrict max_step so pulse is resolved
opts = dict(method="bdf", rtol=1e-6, atol=1e-8, nsteps=500000, max_step=t_pi/50.0)
sol = qt.mesolve(H, psi0, tlist, c_ops=c_ops, e_ops=e_ops, options=opts)

# plot
fig, axes = plt.subplots(GRID_K, GRID_K, figsize=(9, 8), sharex=True, sharey=True)
axes = np.array(axes).reshape(GRID_K, GRID_K)

for r in range(GRID_K):
    for c in range(GRID_K):
        idx = r*GRID_K + c
        ax = axes[r, c]
        ax.plot(tlist*1e6, sol.expect[idx], lw=1.5)
        ax.axvline(t_pi*1e6, color='k', ls='--', lw=0.8)
        ax.set_title(f"Qubit ({r},{c})")
        if r == GRID_K-1: ax.set_xlabel("Time [µs]")
        if c == 0:        ax.set_ylabel("P₁(t)")

fig.suptitle("Simultaneous π pulse on all 9 qubits, T₁/T₂ decay", y=0.98)
plt.tight_layout()
plt.show()

print(f"A_pi [rad/s] = {A_pi:.3e}")
print(f"t_pi  [ns]   = {t_pi*1e9:.3f}")
print(f"T1    [us]   = {T1*1e6:.2f},  T2 [us] = {T2*1e6:.2f}")
print(f"J01   [MHz]  = {J01/(2*np.pi*1e6):.3f} (S_CONST={S_CONST})")
print(f"min dt in tlist [ns] = {dt_min*1e9:.3f}  (should be << t_pi)")
print("Expected P1 decay after t>t_pi: P1(t) ≈ exp(-(t-t_pi)/T1)")
