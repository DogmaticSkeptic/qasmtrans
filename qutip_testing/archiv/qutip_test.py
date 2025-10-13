
"""
qc9_model.py — 9‑qubit (3x3) superconducting-qubit Hamiltonian & QuTiP scaffolding (drive-frame / RWA).

Requirements:
    - QuTiP >= 5.0 (for the new data layer; works with 4.x too if you drop the backend bits)
    - (Optional) qutip-cupy and CuPy if you want GPU

Model overview:
    - 9 qubits in a 3x3 grid; nearest-neighbor couplings (12 edges total).
    - Controls (30 total):
        * 9 microwave complex drives u_k(t) = I_k(t) + i Q_k(t)
          (we pass I_k and Q_k separately as real-valued functions)
        * 9 qubit detunings Δ_k(t) = ω_k(f_k(t)) - ω_d,k    (time-dependent Z terms)
        * 12 tunable exchange couplings J_ij(t) via couplers (XY = (XX + YY))
      Optional: ZZ terms ζ_ij(t).

Hamiltonian (drive frame):
    H(t) = 1/2 Σ_k Δ_k(t) σ_z^k
         + 1/2 Σ_k [ I_k(t) σ_x^k + Q_k(t) σ_y^k ]
         + 1/2 Σ_(i,j in E) J_ij(t) (σ_x^i σ_x^j + σ_y^i σ_y^j)
         + Σ_(i,j in E) ζ_ij(t) σ_z^i σ_z^j             (optional)

Example included:
    - √iSWAP on edge (4,5) (center pair in row-major indexing) using a simple on-resonance exchange pulse
    - X90 on qubit 0
    - Measurement of Z-basis populations

Usage:
    python qc9_model.py            # runs the example schedule (CPU by default)
    python qc9_model.py --gpu      # requests CuPy backend (QuTiP 5 + qutip-cupy required)

Edit the "Pulses" class or pass your own callables to build_H().
"""

from __future__ import annotations
import math
import argparse
import numpy as np

try:
    import qutip as qt
except Exception as exc:
    raise SystemExit("This module requires QuTiP. Install qutip>=5.0 (and qutip-cupy for GPU).") from exc


# ---------- Grid & indexing helpers ----------

def grid_index(r: int, c: int) -> int:
    """Row-major mapping (r,c in {0,1,2}) -> k in {0..8}."""
    return 3*r + c

def edges_3x3() -> list[tuple[int,int]]:
    """Nearest-neighbor edge list for 3x3 grid (12 undirected edges)."""
    E = []
    for r in range(3):
        for c in range(3):
            k = grid_index(r,c)
            if c < 2:   # horizontal neighbor
                E.append((k, grid_index(r,c+1)))
            if r < 2:   # vertical neighbor
                E.append((k, grid_index(r+1,c)))
    return E  # length 12


# ---------- Operator factory ----------

def single_qubit_ops(N: int):
    """Return lists [sigmax_k], [sigmay_k], [sigmaz_k] as Qobj for k=0..N-1."""
    sx_list, sy_list, sz_list = [], [], []
    sx, sy, sz, I2 = qt.sigmax(), qt.sigmay(), qt.sigmaz(), qt.qeye(2)
    for k in range(N):
        ops = [I2]*N
        ops[k] = sx
        sx_list.append(qt.tensor(ops))
        ops[k] = sy
        sy_list.append(qt.tensor(ops))
        ops[k] = sz
        sz_list.append(qt.tensor(ops))
    return sx_list, sy_list, sz_list

def two_qubit_xy_ops(N: int, E: list[tuple[int,int]]):
    """Return list of XY=XX+YY operators (as Qobj) for each edge (i,j)."""
    sx_list, sy_list, sz_list = single_qubit_ops(N)
    XY_ops = []
    for (i,j) in E:
        XX = (sx_list[i] * sx_list[j])
        YY = (sy_list[i] * sy_list[j])
        XY_ops.append(XX + YY)
    return XY_ops

def two_qubit_zz_ops(N: int, E: list[tuple[int,int]]):
    """Return list of ZZ operators for each edge (i,j)."""
    _, _, sz_list = single_qubit_ops(N)
    ZZ_ops = []
    for (i,j) in E:
        ZZ_ops.append(sz_list[i] * sz_list[j])
    return ZZ_ops


# ---------- Control container (30 pulses) ----------

class Pulses:
    """
    Container for 30 time-dependent real-valued callables:
      - I[k](t), Q[k](t) : k=0..8  (9 microwave channels => 18 functions)
      - Delta[k](t)      : k=0..8  (9 detunings)
      - J[e](t)          : e=0..11 (12 exchanges for the 12 edges)
      - Zeta[e](t)       : optional ZZ terms (default 0)

    Each function must accept a float t and return a float.

    Defaults are zeros (no drive, no coupling). Fill as needed.
    """
    def __init__(self, N=9, E=None):
        self.N = N
        self.E = list(range(len(edges_3x3()))) if E is None else list(range(len(E)))
        zero = (lambda t: 0.0)

        self.I = [zero for _ in range(N)]
        self.Q = [zero for _ in range(N)]
        self.Delta = [zero for _ in range(N)]
        self.J = [zero for _ in self.E]
        self.Zeta = [zero for _ in self.E]  # optional

    # convenience setters
    def set_IQ(self, k, I_func, Q_func):
        self.I[k] = I_func
        self.Q[k] = Q_func
    def set_Delta(self, k, func):
        self.Delta[k] = func
    def set_J(self, e_index, func):
        self.J[e_index] = func
    def set_Zeta(self, e_index, func):
        self.Zeta[e_index] = func


# ---------- Hamiltonian assembly (drive-frame / RWA) ----------

def build_H(pulses: Pulses, include_ZZ: bool = False):
    """
    Return QuTiP time-dependent Hamiltonian list for N=9.
    H(t) = 1/2 sum_k Delta_k(t) Z_k
         + 1/2 sum_k [ I_k(t) X_k + Q_k(t) Y_k ]
         + 1/2 sum_e J_e(t) XY_e
         + sum_e Zeta_e(t) ZZ_e            (optional)
    """
    N = pulses.N
    E_pairs = edges_3x3()
    sx_list, sy_list, sz_list = single_qubit_ops(N)
    XY_ops = two_qubit_xy_ops(N, E_pairs)
    ZZ_ops = two_qubit_zz_ops(N, E_pairs) if include_ZZ else None

    H = []  # time-dependent list form

    # 1) Detunings (Z)
    for k in range(N):
        H.append([0.5 * sz_list[k], pulses.Delta[k]])

    # 2) IQ drives
    for k in range(N):
        H.append([0.5 * sx_list[k], pulses.I[k]])
        H.append([0.5 * sy_list[k], pulses.Q[k]])

    # 3) XY exchanges on edges
    for e_idx, (i,j) in enumerate(E_pairs):
        H.append([0.5 * XY_ops[e_idx], pulses.J[e_idx]])

    # 4) Optional ZZ
    if include_ZZ:
        for e_idx, (i,j) in enumerate(E_pairs):
            H.append([ZZ_ops[e_idx], pulses.Zeta[e_idx]])

    return H


# ---------- Example schedule (√iSWAP on center edge + X90 on qubit 0) ----------

def example_pulses(T=100e-9, tlist=None):
    """
    Build a minimal example:
      - √iSWAP on edge (4,5) (center pair in row-major indexing).
      - X90 on qubit 0 via a simple resonant I-quadrature pulse.
      - All other controls set to 0.
    """
    N = 9
    E_pairs = edges_3x3()
    center_edge = None
    # pick an edge with nodes 4 and 5 (middle row, middle and right qubits in row-major [0..8])
    for idx, (i,j) in enumerate(E_pairs):
        if {i,j} == {4,5}:
            center_edge = idx
            break
    if center_edge is None:
        raise RuntimeError("Center edge (4,5) not found—check edges_3x3().")

    pulses = Pulses(N=N, E=E_pairs)

    # Time window [0, T]
    if tlist is None:
        tlist = np.linspace(0.0, T, 501)

    # Exchange rate for √iSWAP: want g_ex * T = π/2  => g_ex = π/(2T)
    g_ex = math.pi / (2*T)

    # Use a cosine ramped square for J(t) on the chosen edge
    def flat_top(t, T=T, Tr=6e-9):
        # cosine ramps of duration Tr at both ends
        if t < 0 or t > T:
            return 0.0
        if t < Tr:
            x = t/Tr
            return 0.5 - 0.5*math.cos(math.pi*x)
        if t > T - Tr:
            x = (T - t)/Tr
            return 0.5 - 0.5*math.cos(math.pi*x)
        return 1.0

    pulses.set_J(center_edge, lambda t, g=g_ex: g*flat_top(t))

    # Simple X90 on qubit 0: Ω * T1q = π/2  (let's use T1q = T/5 here)
    T1q = T/5
    Om = (math.pi/2) / T1q
    def I_drive(t, T1q=T1q, Om=Om):
        if 0 <= t <= T1q:
            return Om   # resonant I-only (phase 0)
        return 0.0
    def Q_drive(t):  # no Q component
        return 0.0

    pulses.set_IQ(0, I_drive, Q_drive)

    # Optionally add tiny detuning drift on all qubits (here zero)
    for k in range(N):
        pulses.set_Delta(k, lambda t: 0.0)

    return pulses, tlist


# ---------- Measurement helpers ----------

def z_projectors(N=9):
    """Return lists [Pi0_k], [Pi1_k] for k=0..N-1 in the computational basis."""
    b0 = qt.basis(2,0); b1 = qt.basis(2,1)
    Pi0_1q = [b0*b0.dag(), b1*b1.dag()]  # index 0 is |0><0|, 1 is |1><1|
    Pi0_list, Pi1_list = [], []
    for k in range(N):
        ops = [qt.qeye(2)]*N
        ops[k] = Pi0_1q[0]
        Pi0_list.append(qt.tensor(ops))
        ops[k] = Pi0_1q[1]
        Pi1_list.append(qt.tensor(ops))
    return Pi0_list, Pi1_list


# ---------- Example simulation ----------

def run_example(cpu_or_gpu: str = "cpu"):
    if cpu_or_gpu.lower() == "gpu":
        try:
            qt.set_data_default("cupy")
            print("[info] Using CuPy backend (GPU).")
        except Exception as exc:
            print("[warn] Could not enable CuPy backend; falling back to CPU.", exc)

    N = 9
    pulses, tlist = example_pulses(T=100e-9)

    psi0 = qt.tensor([qt.basis(2, 0) for _ in range(N)])
    rho0 = psi0 * psi0.dag()

    H = build_H(pulses, include_ZZ=False)

    c_ops = []

    Pi0, Pi1 = z_projectors(N=N)
    e_ops = [Pi1[0], Pi1[4], Pi1[5]]

    opts = {
        "nsteps": 100000,
        "atol": 1e-10,
        "rtol": 1e-8,
        "store_states": False,
        "progress_bar": None
    }

    result = qt.mesolve(H, rho0, tlist, c_ops, e_ops=e_ops, options=opts)

    expect_array = np.asarray(result.expect, dtype=float)
    if expect_array.ndim == 2:
        traces = expect_array.T
    else:
        traces = np.vstack([np.asarray(v, dtype=float).squeeze() for v in result.expect]).T

    return tlist, traces


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpu", action="store_true", help="Use CuPy backend (QuTiP 5 + qutip-cupy required).")
    args = ap.parse_args()
    mode = "gpu" if args.gpu else "cpu"
    tlist, traces = run_example(mode)
    # Simple console preview
    print(f"Simulated {len(tlist)} steps. Example: P(|1>) on qubits [0,4,5] at final time:",
          traces[-1, 0], traces[-1, 1], traces[-1, 2])

