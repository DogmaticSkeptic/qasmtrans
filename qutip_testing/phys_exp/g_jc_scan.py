import numpy as np
import matplotlib.pyplot as plt

Phi0 = 2.067833848e-15

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

Phi2_fixed = 0.14355 * Phi0

def EJ_of_flux(EJ_sum, alpha, Phi):
    phi = np.pi * Phi / Phi0
    return EJ_sum * np.sqrt(np.cos(phi)**2 + (alpha**2)*(np.sin(phi)**2))

def w01_from_EJ_EC(EJ, EC):
    return np.sqrt(8.0*EJ*EC) - EC

def EJ_from_w01_EC(w01, EC):
    return (w01 + EC)**2 / (8.0*EC)

def xi_from_EJ_EC(EJ, EC):
    return 2.0*EC / EJ

def wc_of_flux(Phic):
    EJc = EJ_of_flux(EJc_sum, alphac, Phic)
    return w01_from_EJ_EC(EJc, ECc)

def w2_of_flux(Phi2):
    EJ2 = EJ_of_flux(EJ2_sum, alpha2, Phi2)
    return w01_from_EJ_EC(EJ2, EC2)

def gjc_via_A16_normalized(g_ref, EJj, ECj, Phic):
    EJc_Phi  = EJ_of_flux(EJc_sum, alphac, Phic)
    EJc_0    = EJ_of_flux(EJc_sum, alphac, 0.0)
    xi_j     = xi_from_EJ_EC(EJj, ECj)
    xi_c_Phi = xi_from_EJ_EC(EJc_Phi, ECc)
    xi_c_0   = xi_from_EJ_EC(EJc_0,   ECc)
    f_Phi = (EJc_Phi/ECc)**0.25 * (1.0 - 0.125*(xi_c_Phi + xi_j))
    f_0   = (EJc_0  /ECc)**0.25 * (1.0 - 0.125*(xi_c_0   + xi_j))
    return g_ref * (f_Phi / f_0)

def g01_full(Phi2, Phic, EJ1, EC1):
    w2 = w2_of_flux(Phi2)
    wc = wc_of_flux(Phic)
    EJ2 = EJ_of_flux(EJ2_sum, alpha2, Phi2)
    g1c = gjc_via_A16_normalized(g1c_ref, EJ1, EC1, Phic)
    g2c = gjc_via_A16_normalized(g2c_ref, EJ2, EC2, Phic)
    D1 = wc - w1
    D2 = wc - w2
    L1 = wc + w1
    L2 = wc + w2
    eps = 1e-30
    invs = (1.0/max(D1,eps) + 1.0/max(D2,eps) + 1.0/max(L1,eps) + 1.0/max(L2,eps))
    return g12 - 0.5*g1c*g2c*invs, g1c, g2c, w2, wc, D1, D2, L1, L2

EJ1 = EJ_from_w01_EC(w1, EC1)

N = 2001
Phic_vals = np.linspace(-0.35, 0.35, N) * Phi0
g01_list = []
g1c_list = []
g2c_list = []
wc_list  = []
w2_fixed = w2_of_flux(Phi2_fixed)

for Phic in Phic_vals:
    g01, g1c, g2c, w2, wc, D1, D2, L1, L2 = g01_full(Phi2_fixed, Phic, EJ1, EC1)
    g01_list.append(g01)
    g1c_list.append(g1c)
    g2c_list.append(g2c)
    wc_list.append(wc)

g01_arr = np.array(g01_list)
g1c_arr = np.array(g1c_list)
g2c_arr = np.array(g2c_list)
wc_arr  = np.array(wc_list)

MHz = 2*np.pi*1e6
GHz = 2*np.pi*1e9
idx_min = np.argmin(np.abs(g01_arr))
zero_cross_exists = np.any(np.signbit(g01_arr[1:]*g01_arr[:-1]))

print("=== Scan summary ===")
print("Phi2_fixed / Phi0:", Phi2_fixed / Phi0)
print("w2_fixed / 2π [GHz]:", w2_fixed/GHz)
print("min |g01| [MHz]:", np.abs(g01_arr[idx_min])/MHz)
print("at Phic / Phi0:", Phic_vals[idx_min]/Phi0)
print("zero_cross_exists:", bool(zero_cross_exists))
print("g01(Phic=0) [MHz]:", g01_full(Phi2_fixed, 0.0, EJ1, EC1)[0]/MHz)
print("wc range [GHz]:", np.min(wc_arr)/GHz, np.max(wc_arr)/GHz)

plt.figure(figsize=(7,4))
plt.plot(Phic_vals/Phi0, g01_arr/MHz, lw=1.2, label="g01(Φc) [MHz]")
plt.axhline(0.0, color='k', lw=0.8)
plt.xlabel("Φc / Φ0")
plt.ylabel("g01 [MHz]")
plt.title("Effective exchange coupling g01 vs Φc")
plt.legend()
plt.tight_layout()

plt.figure(figsize=(7,4))
plt.plot(Phic_vals/Phi0, g1c_arr/MHz, lw=1.0, label="g1c(Φc) [MHz]")
plt.plot(Phic_vals/Phi0, g2c_arr/MHz, lw=1.0, label="g2c(Φc) [MHz]", linestyle="--")
plt.xlabel("Φc / Φ0")
plt.ylabel("g_jc [MHz]")
plt.title("Qubit–coupler couplings vs Φc")
plt.legend()
plt.tight_layout()

plt.figure(figsize=(7,4))
plt.plot(Phic_vals/Phi0, wc_arr/GHz, lw=1.2, label="wc(Φc) [GHz]")
plt.axhline(w1/GHz, color='C1', lw=1.0, linestyle='--', label="w1 [GHz]")
plt.axhline(w2_fixed/GHz, color='C2', lw=1.0, linestyle=':', label="w2 fixed [GHz]")
plt.xlabel("Φc / Φ0")
plt.ylabel("Frequency [GHz]")
plt.title("Coupler and qubit frequencies vs Φc")
plt.legend()
plt.tight_layout()

plt.show()
