import json
import time
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

# ---------------- user settings ----------------
CALIB_PATH  = "calibration_params.json"
MAX_QUBITS  = 9            # add one qubit at a time up to this many (≤ 9)
GRID_K      = 3            # 3×3 grid topology; we populate the first N sites
NPTS        = 201          # time samples per simulation
S_CONST     = 1.0          # constant exchange envelope s ∈ [0, 1]

# ---------------- load calibration ----------------
def load_calibration(path):
    with open(path, "r") as f:
        C = json.load(f)
    A_pi = float(C["single_qubit_drive"]["A_gaussian_I_peak_pi"])
    T1   = float(C["dissipation"]["T1_s"])
    T2   = float(C["dissipation"]["T2_s"])
    b01  = float(C["linearized_couplings_at_on_cancel"]["b01_rad_per_s_per_Wb"])
    Phi0 = 2.067833848e-15
    dphi = float(C["biases"]["delta_phic_on_minus_off_over_phi0_abs"]) * Phi0
    return {"A_pi": A_pi, "T1": T1, "T2": T2, "b01": b01, "Phi0": Phi0, "dphi": dphi}

CAL = load_calibration(CALIB_PATH)

# ---------------- single-qubit (2-level) operators ----------------
sx = qt.sigmax()
sy = qt.sigmay()
sz = qt.sigmaz()
sm = qt.sigmam()
I2 = qt.qeye(2)

def lift_one(op, site, N):
    ops = [I2]*N
    ops = list(ops)
    ops[site] = op
    return qt.tensor(ops)

def lift_two(op_i, i, op_j, j, N):
    ops = [I2]*N
    ops = list(ops)
    ops[i] = op_i
    ops[j] = op_j
    return qt.tensor(ops)

# ---------------- grid helpers ----------------
def grid_order(k=3):
    # row-major: (0,0)->0, (0,1)->1, (0,2)->2, (1,0)->3, ...
    return [r*k + c for r in range(k) for c in range(k)]

def edges_for_active(active_nodes, k=3):
    # nearest neighbors among currently active nodes; compactly reindexed
    active = set(active_nodes)
    def idx(r,c): return r*k + c
    E = []
    for r in range(k):
        for c in range(k):
            u = idx(r,c)
            if u not in active:
                continue
            if c+1 < k:
                v = idx(r,c+1)
                if v in active and u < v:
                    E.append((u, v))
            if r+1 < k:
                v = idx(r+1,c)
                if v in active and u < v:
                    E.append((u, v))
    remap = {node:i for i,node in enumerate(active_nodes)}
    return [(remap[i], remap[j]) for (i,j) in E]

# ---------------- Hamiltonian assembly (2-level) ----------------
def build_H(N, edges, A_drive, J01):
    """
    Rotating frame, RWA.
    H = sum_i 0.5*A_drive * X_i + sum_(i,j) (J01/2)*(X_i X_j + Y_i Y_j)
    """
    H = 0
    for i in range(N):
        H += 0.5*A_drive * lift_one(sx, i, N)
    for (i,j) in edges:
        XX = lift_two(sx, i, sx, j, N)
        YY = lift_two(sy, i, sy, j, N)
        H += 0.5*J01*(XX + YY)
    return H

# ---------------- collapse operators (T1, T2) ----------------
def build_c_ops(N, T1, T2):
    gamma = 1.0/float(T1)
    gamma_phi = max(0.0, 1.0/float(T2) - 0.5*gamma)
    c_ops = []
    for i in range(N):
        # relaxation
        c_ops.append(np.sqrt(gamma) * lift_one(sm, i, N))
        # pure dephasing: choose L = sqrt(gamma_phi/2) * sz so that coherence decays at rate gamma_phi
        if gamma_phi > 0.0:
            c_ops.append(np.sqrt(gamma_phi/2.0) * lift_one(sz, i, N))
    return c_ops

# ---------------- initial state and observables ----------------
def initial_state(N):
    return qt.tensor([qt.basis(2,0)]*N)

def observables(N):
    # measure P1 = |1><1| on each qubit
    P1 = (I2 - sz)/2.0
    return [lift_one(P1, i, N) for i in range(N)]

# ---------------- incremental benchmark ----------------
def run_incremental(cal, max_qubits=6, k=3, npts=201, s_const=1.0):
    order = grid_order(k)
    A_drive = cal["A_pi"]                      # rad/s
    t_final = np.pi/float(A_drive)             # π pulse length
    tlist = np.linspace(0.0, t_final, npts)

    # exchange strength from linearization
    J01 = cal["b01"] * (cal["Phi0"]*cal["dphi"]) * float(s_const)   # rad/s

    Ns, times = [], []
    for N in range(1, max_qubits+1):
        active = order[:N]
        edges  = edges_for_active(active, k)

        H     = build_H(N, edges, A_drive, J01)
        c_ops = build_c_ops(N, cal["T1"], cal["T2"])
        psi0  = initial_state(N)
        e_ops = observables(N)

        opts = {"method": "bdf", "rtol": 1e-5, "atol": 1e-7, "nsteps": 200000}

        t0 = time.perf_counter()
        _  = qt.mesolve(H, psi0, tlist, c_ops=c_ops, e_ops=e_ops, options=opts)
        t1 = time.perf_counter()

        elapsed = t1 - t0
        print(f"N={N:2d}  edges={edges}  time={elapsed:.3f} s")
        Ns.append(N); times.append(elapsed)

    return Ns, times, t_final, J01, A_drive

Ns, times, t_final, J01, A_drive = run_incremental(CAL, max_qubits=MAX_QUBITS, k=GRID_K, npts=NPTS, s_const=S_CONST)

plt.figure(figsize=(6,4))
plt.plot(Ns, times, marker="o")
plt.xlabel("Number of qubits N (added in 3×3 grid order)")
plt.ylabel("Wall time [s]")
plt.title("mesolve runtime vs N (2-level qubits, π pulse, exchange + T1/T2)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print(f"A_pi [rad/s] = {A_drive:.3e}")
print(f"t_pi  [ns]   = {t_final*1e9:.3f}")
print(f"J01   [MHz]  = {J01/(2*np.pi*1e6):.3f}")
