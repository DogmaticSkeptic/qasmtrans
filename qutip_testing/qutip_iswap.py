from __future__ import annotations
import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

def projectors_and_sx(mode, dims):
    d = dims[mode]
    P01_loc = qt.basis(d,1)*qt.basis(d,1).dag()
    P12_loc = qt.basis(d,2)*qt.basis(d,2).dag() if d >= 3 else 0*qt.qeye(d)
    s01_loc = qt.basis(d,0)*qt.basis(d,1).dag()
    s12_loc = qt.basis(d,1)*qt.basis(d,2).dag() if d >= 3 else 0*qt.qeye(d)
    P01 = qt.tensor([qt.qeye(dims[k]) if k != mode else P01_loc for k in range(len(dims))])
    P12 = qt.tensor([qt.qeye(dims[k]) if k != mode else P12_loc for k in range(len(dims))])
    s01 = qt.tensor([qt.qeye(dims[k]) if k != mode else s01_loc for k in range(len(dims))])
    s12 = qt.tensor([qt.qeye(dims[k]) if k != mode else s12_loc for k in range(len(dims))])
    sx = s01 + s01.dag() + (np.sqrt(2.0)*(s12 + s12.dag()) if d >= 3 else 0*s01)
    return P01, P12, sx

def EJ_of_flux(EJS, EJL, phi_e):
    return np.sqrt(EJS*EJS + EJL*EJL + 2.0*EJS*EJL*np.cos(phi_e))

def omega01_of_flux(EC, EJ_eff):
    xi = 2.0*EC/np.maximum(EJ_eff, 1e-18)
    return np.sqrt(8.0*EJ_eff*EC) - EC*(1.0 + 0.25*xi)

def g_jc_of_flux(Ejc, EJj, ECj, EJc_eff, ECc):
    term = (EJj/ECj)*(EJc_eff/ECc)
    pref = np.sqrt(2.0)*Ejc*np.power(np.maximum(term, 1e-24), 0.25)
    xij = 2.0*ECj/np.maximum(EJj, 1e-18)
    xic = 2.0*ECc/np.maximum(EJc_eff, 1e-18)
    return pref*(1.0 - 0.125*(xij + xic))

def flattop_cosine(T, tau):
    def u(t):
        if t <= 0.0 or t >= T:
            return 0.0
        if t < tau:
            return 0.5*(1.0 - np.cos(np.pi*t/tau))
        if t <= T - tau:
            return 1.0
        x = (t - (T - tau))/tau
        return 0.5*(1.0 + np.cos(np.pi*x))
    return u

def pulses_inverted_qubit(T, tau, phi1_dc, phi2_dc, phic_dc, dphi2_peak, dphic_peak):
    uC = flattop_cosine(T, tau)
    def phi1_t(t):
        return phi1_dc
    def phi2_t(t):
        return phi2_dc - dphi2_peak*uC(t)
    def phic_t(t):
        return phic_dc + dphic_peak*uC(t)
    return phi1_t, phi2_t, phic_t, uC

class FluxDrivenHamiltonian:
    def __init__(self, phys, dims=(3,2,3)):
        self.dims = dims
        self.EC1 = phys["EC1"]
        self.EC2 = phys["EC2"]
        self.ECc = phys["ECc"]
        self.EJS1 = phys["EJS1"]
        self.EJL1 = phys["EJL1"]
        self.EJS2 = phys["EJS2"]
        self.EJL2 = phys["EJL2"]
        self.EJSc = phys["EJSc"]
        self.EJLc = phys["EJLc"]
        self.Ej1c = phys["Ej1c"]
        self.Ej2c = phys["Ej2c"]
        self.g12  = phys["g12"]
        self.P01_1, self.P12_1, self.sx1 = projectors_and_sx(0, dims)
        self.P01_c, self.P12_c, self.sxc = projectors_and_sx(1, dims)
        self.P01_2, self.P12_2, self.sx2 = projectors_and_sx(2, dims)
        self.Hxy_dir = self.g12*self.sx1*self.sx2

    def coeffs_from_flux(self, tlist, phi1_t, phi2_t, phic_t, w1_ref, w2_ref, wc_ref):
        phi1 = np.array([phi1_t(t) for t in tlist])
        phi2 = np.array([phi2_t(t) for t in tlist])
        phic = np.array([phic_t(t) for t in tlist])
        EJ1 = EJ_of_flux(self.EJS1, self.EJL1, phi1)
        EJ2 = EJ_of_flux(self.EJS2, self.EJL2, phi2)
        EJc = EJ_of_flux(self.EJSc, self.EJLc, phic)
        w1 = omega01_of_flux(self.EC1, EJ1)
        w2 = omega01_of_flux(self.EC2, EJ2)
        wc = omega01_of_flux(self.ECc, EJc)
        dw1 = w1 - w1_ref
        dw2 = w2 - w2_ref
        dwc = wc - wc_ref
        dw1_12 = 2.0*dw1 - self.EC1
        dw2_12 = 2.0*dw2 - self.EC2
        dwc_12 = 2.0*dwc - self.ECc
        g1c = g_jc_of_flux(self.Ej1c, EJ1, self.EC1, EJc, self.ECc)
        g2c = g_jc_of_flux(self.Ej2c, EJ2, self.EC2, EJc, self.ECc)
        return dw1, dw2, dwc, dw1_12, dw2_12, dwc_12, g1c, g2c, w1, w2, wc

    def build_list_with_arrays(self, dw1, dw2, dwc, dw1_12, dw2_12, dwc_12, g1c, g2c):
        H = [self.Hxy_dir,
             [self.P01_1, np.asarray(dw1, dtype=float)],
             [self.P12_1, np.asarray(dw1_12, dtype=float)],
             [self.P01_2, np.asarray(dw2, dtype=float)],
             [self.P12_2, np.asarray(dw2_12, dtype=float)],
             [self.P01_c, np.asarray(dwc, dtype=float)],
             [self.P12_c, np.asarray(dwc_12, dtype=float)],
             [self.sx1*self.sxc, np.asarray(g1c, dtype=float)],
             [self.sxc*self.sx2, np.asarray(g2c, dtype=float)]]
        return H

def init_state_10(dims):
    d1, dc, d2 = dims
    return qt.tensor(qt.basis(d1,1), qt.basis(dc,0), qt.basis(d2,0))

def reduce_to_qubits(rho, dims):
    rho_qs = qt.ptrace(rho, [0,2])
    d = dims[0]
    Pq = qt.basis(d,0)*qt.basis(d,0).dag() + qt.basis(d,1)*qt.basis(d,1).dag()
    P = qt.tensor(Pq, Pq)
    return P*rho_qs*P

def run_once_and_plot(phys, dims, T, tau, dt, phi1_dc, phi2_dc, phic_dc, dphi2_peak, dphic_peak):
    Hdrv = FluxDrivenHamiltonian(phys, dims)
    tlist = np.arange(0.0, T + 0.5*dt, dt, dtype=float)
    phi1_t, phi2_t, phic_t, uC = pulses_inverted_qubit(T, tau, phi1_dc, phi2_dc, phic_dc, dphi2_peak, dphic_peak)
    EJ1_ref = EJ_of_flux(phys["EJS1"], phys["EJL1"], phi1_dc)
    EJ2_ref = EJ_of_flux(phys["EJS2"], phys["EJL2"], phi2_dc)
    EJc_ref = EJ_of_flux(phys["EJSc"], phys["EJLc"], phic_dc)
    w1_ref = omega01_of_flux(phys["EC1"], EJ1_ref)
    w2_ref = omega01_of_flux(phys["EC2"], EJ2_ref)
    wc_ref = omega01_of_flux(phys["ECc"], EJc_ref)
    dw1, dw2, dwc, dw1_12, dw2_12, dwc_12, g1c, g2c, w1_arr, w2_arr, wc_arr = Hdrv.coeffs_from_flux(tlist, phi1_t, phi2_t, phic_t, w1_ref, w2_ref, wc_ref)

    GHz = 1.0/(2.0*np.pi*1e9)
    MHz = 1.0/(2.0*np.pi*1e6)
    print("max |Δω2| [MHz] =", np.max(np.abs(dw2))*MHz)
    print("max |Δωc| [MHz] =", np.max(np.abs(dwc))*MHz)
    H = Hdrv.build_list_with_arrays(dw1, dw2, dwc, dw1_12, dw2_12, dwc_12, g1c, g2c)
    psi0 = init_state_10(dims)
    d1, dc, d2 = dims
    ket100 = qt.tensor(qt.basis(d1,1), qt.basis(dc,0), qt.basis(d2,0))
    ket001 = qt.tensor(qt.basis(d1,0), qt.basis(dc,0), qt.basis(d2,1))
    ket101 = qt.tensor(qt.basis(d1,1), qt.basis(dc,0), qt.basis(d2,1))
    P100 = ket100*ket100.dag()
    P001 = ket001*ket001.dag()
    P101 = ket101*ket101.dag()
    Pq = qt.basis(d1,0)*qt.basis(d1,0).dag() + qt.basis(d1,1)*qt.basis(d1,1).dag()
    Pq2 = qt.basis(d2,0)*qt.basis(d2,0).dag() + qt.basis(d2,1)*qt.basis(d2,1).dag()
    P_comp = qt.tensor(Pq, qt.qeye(dc), Pq2)
    e_ops = [P100, P001, P101, P_comp]
    opts = qt.Options(nsteps=120000, rtol=1e-7, atol=1e-9, store_states=True, progress_bar=None)
    res = qt.sesolve(H, psi0, tlist, e_ops=e_ops, options=opts)
    rec = {}
    rec["t"] = tlist
    rec["P10"] = np.array(res.expect[0])
    rec["P01"] = np.array(res.expect[1])
    rec["P11"] = np.array(res.expect[2])
    rec["Pcomp"] = np.array(res.expect[3])
    rec["phi2"] = np.array([phi2_t(t) for t in tlist])
    rec["phic"] = np.array([phic_t(t) for t in tlist])
    rec["w1"] = w1_arr
    rec["w2"] = w2_arr
    rec["wc"] = wc_arr
    rec["g1c"] = g1c
    rec["g2c"] = g2c
    t_ns = 1e9*rec["t"]
    plt.figure()
    plt.plot(t_ns, rec["P10"], label="P_10")
    plt.plot(t_ns, rec["P01"], label="P_01")
    plt.plot(t_ns, rec["P11"], label="P_11")
    plt.plot(t_ns, rec["Pcomp"], label="Tr P_comp")
    plt.xlabel("time ns")
    plt.ylabel("population")
    plt.legend()
    plt.tight_layout()
    plt.figure()
    plt.plot(t_ns, rec["phic"], label="phi_c")
    plt.plot(t_ns, rec["phi2"], label="phi_2")
    plt.xlabel("time ns")
    plt.ylabel("flux rad")
    plt.legend()
    plt.tight_layout()
    GHz = 1.0/(2.0*np.pi*1e9)
    plt.figure()
    plt.plot(t_ns, rec["w1"]*GHz, label="w1_2pi_GHz")
    plt.plot(t_ns, rec["w2"]*GHz, label="w2_2pi_GHz")
    plt.plot(t_ns, rec["wc"]*GHz, label="wc_2pi_GHz")
    plt.xlabel("time ns")
    plt.ylabel("frequency GHz")
    plt.legend()
    plt.tight_layout()
    MHz = 1.0/(2.0*np.pi*1e6)
    plt.figure()
    plt.plot(t_ns, rec["g1c"]*MHz, label="g1c_2pi_MHz")
    plt.plot(t_ns, rec["g2c"]*MHz, label="g2c_2pi_MHz")
    plt.xlabel("time ns")
    plt.ylabel("coupling MHz")
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    twopi = 2.0*np.pi
    ECq = twopi*230.0e6
    ECc = twopi*250.0e6
    EJ1 = (twopi*3.50e9 + ECq)**2/(8.0*ECq)
    EJ2 = (twopi*3.52e9 + ECq)**2/(8.0*ECq)
    EJc = (twopi*5.90e9 + ECc)**2/(8.0*ECc)
    asym_q = 0.1
    asym_c = 0.1
    EJS1 = EJ1*(1.0+asym_q)/2.0
    EJL1 = EJ1*(1.0-asym_q)/2.0
    EJS2 = EJ2*(1.0+asym_q)/2.0
    EJL2 = EJ2*(1.0-asym_q)/2.0
    EJSc = EJc*(1.0+asym_c)/2.0
    EJLc = EJc*(1.0-asym_c)/2.0
    Ej1c = twopi*8.0e6
    Ej2c = twopi*8.0e6
    g12 = twopi*4.0e6
    phys = {"EC1":ECq, "EC2":ECq, "ECc":ECc,
            "EJS1":EJS1, "EJL1":EJL1, "EJS2":EJS2, "EJL2":EJL2,
            "EJSc":EJSc, "EJLc":EJLc, "Ej1c":Ej1c, "Ej2c":Ej2c, "g12":g12}
    dims = (3,2,3)
    T = 60e-9
    tau = 10e-9
    dt = 0.5e-9
    phi1_dc = 0.5
    phi2_dc = 0.5
    phic_dc = 0.5
    dphi2_peak = 0.001*np.pi
    dphic_peak = 0.001*np.pi
    run_once_and_plot(phys, dims, T, tau, dt, phi1_dc, phi2_dc, phic_dc, dphi2_peak, dphic_peak)
