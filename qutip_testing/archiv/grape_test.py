import numpy as np
import matplotlib.pyplot as plt
from qutip import Qobj, basis, tensor, qeye
import qutip.control.pulseoptim as cpo

def sx01():
    e0 = basis(3, 0)
    e1 = basis(3, 1)
    return e0*e1.dag() + e1*e0.dag()

def sy01():
    e0 = basis(3, 0)
    e1 = basis(3, 1)
    return -1j*e0*e1.dag() + 1j*e1*e0.dag()

def exch_10_01():
    e0 = basis(3, 0)
    e1 = basis(3, 1)
    return tensor(e1*e0.dag(), e0*e1.dag()) + tensor(e0*e1.dag(), e1*e0.dag())

def embed_u4_in_u9(U4):
    U = np.eye(9, dtype=complex)
    U[:4, :4] = U4
    return Qobj(U, dims=[[3, 3], [3, 3]])

def random_su4(rng):
    Z = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    Q, R = np.linalg.qr(Z)
    d = np.diag(R)
    ph = d/np.abs(d)
    U = Q*ph
    detU = np.linalg.det(U)
    U = U / detU**0.25
    return U

rng = np.random.default_rng(1)

H_d = 0.0 * tensor(qeye(3), qeye(3))

Hc_I1 = tensor(sx01(), qeye(3))
Hc_Q1 = tensor(sy01(), qeye(3))
Hc_I2 = tensor(qeye(3), sx01())
Hc_Q2 = tensor(qeye(3), sy01())
Hc_g01 = exch_10_01()
H_controls = [Hc_I1, Hc_Q1, Hc_I2, Hc_Q2, Hc_g01]

U_0 = tensor(qeye(3), qeye(3))
U4_targ = random_su4(rng)
U_targ = embed_u4_in_u9(U4_targ)

num_tslots = 200
evo_time = 10.0

res = cpo.optimize_pulse_unitary(
    H_d,
    H_controls,
    U_0,
    U_targ,
    num_tslots=num_tslots,
    evo_time=evo_time,
    amp_lbound=[-5.0, -5.0, -5.0, -5.0, 0.0],
    amp_ubound=[5.0, 5.0, 5.0, 5.0, 5.0],
    init_pulse_type="CRAB_FOURIER",
    init_pulse_params={"num_waves": 3},
    ramping_pulse_type="GAUSSIAN_EDGE",
    alg="CRAB",
    fid_err_targ=1e-4,
    max_iter=800,
    gen_stats=True
)

amps = np.asarray(res.final_amps)
if amps.shape[0] == num_tslots:
    amps = amps.T

t_edges = np.asarray(res.time)
t_mid = 0.5 * (t_edges[:-1] + t_edges[1:])

labels = ["I1", "Q1", "I2", "Q2", "g01"]
plt.figure()
for k in range(5):
    plt.step(t_mid, amps[k], where="mid", label=labels[k])
plt.xlabel("time")
plt.ylabel("amplitude")
plt.title("Optimized controls: I1, Q1, I2, Q2, g01")
plt.legend()
plt.tight_layout()

print("fid_err:", res.fid_err)
plt.show()
