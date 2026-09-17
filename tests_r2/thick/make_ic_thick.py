#!/usr/bin/env python3
"""Build the initial stratification for the OPTICALLY THICK spherical two-stream test.

The state handed to the code is the EXACT radiative-equilibrium diffusion solution for a
power-law density and a constant opacity, so the correct radiative flux on it is known
analytically at every face:

    rho(r)      = rho_in (r_in/r)^n
    d(T^4)/dr   = -3 kappa rho(r) L / (16 pi sigma r^2),     T(top) = Teff
    L           = 4 pi r_out^2 sigma Teff^4
  =>  F_diff(r) = -(16 sigma T^3)/(3 kappa rho) dT/dr = L/(4 pi r^2) = F_req(r).

kappa is chosen from the requested TOTAL radial Rosseland depth
    tau_tot = int_{r_in}^{r_out} kappa rho dr.
The T^4 integral is done NUMERICALLY (trapezoid on a fine grid), so nothing here assumes
a closed form and any n -- including n = 1, where the deep power law degenerates -- works.

Ideal gas, <units> all 1.0:  p = rho k_B T / (mu m_H),  eint = p/(gamma-1).

Writes three columns "r[cm] rho eint" ascending in r, spanning the mesh AND its radial
ghosts with a 5-cell margin at each end, and prints the athinput overrides the run needs.
"""
import argparse
import numpy as np

SIGMA = 5.670374419e-5
KB = 1.380649e-16
MH = 1.67262192369e-24


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=float, required=True, help="density power law index")
    ap.add_argument("--tau", type=float, required=True, help="total radial tau_Ross")
    ap.add_argument("--ratio", type=float, required=True, help="r_out/r_in")
    ap.add_argument("--rin", type=float, default=1.0e12)
    ap.add_argument("--rho-in", type=float, default=1.0e-6)
    ap.add_argument("--teff", type=float, default=5000.0)
    ap.add_argument("--mu", type=float, default=1.0)
    ap.add_argument("--gamma", type=float, default=5.0/3.0)
    ap.add_argument("--nx1", type=int, default=128)
    ap.add_argument("--nfine", type=int, default=20001)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rin, rout = a.rin, a.rin*a.ratio
    dr = (rout - rin)/a.nx1
    rlo, rhi = rin - 5.0*dr, rout + 5.0*dr       # mesh + ghosts, with margin

    # kappa from the requested total optical depth (integral over the ACTIVE domain)
    rq = np.linspace(rin, rout, 400001)
    integ = np.trapezoid(a.rho_in*(rin/rq)**a.n, rq)
    kappa = a.tau/integ

    lum = 4.0*np.pi*rout**2*SIGMA*a.teff**4

    # T^4 by inward integration of the diffusion equation on the extended range
    r = np.linspace(rlo, rhi, a.nfine)
    rho = a.rho_in*(rin/r)**a.n
    g = 3.0*kappa*rho*lum/(16.0*np.pi*SIGMA*r**2)     # = -d(T^4)/dr
    # cumulative integral inward from the TOP OF THE EXTENDED RANGE, where the anchor
    # T = Teff is placed.  Anchoring at r_out instead would drive T^4 negative in the
    # outer ghost margin for the thick/thin-shell cases; the anchor value is irrelevant
    # to the test, because the ODE -- and hence F_diff = F_req -- holds for any of them.
    c = np.concatenate(([0.0], np.cumsum(0.5*(g[1:] + g[:-1])*np.diff(r))))
    t4 = a.teff**4 + (c[-1] - c)
    temp = t4**0.25
    pres = rho*(KB/(a.mu*MH))*temp
    eint = pres/(a.gamma - 1.0)

    with open(a.out, "w") as f:
        f.write("# r[cm] rho[g/cm^3] eint[erg/cm^3]   n=%g tau=%g rout/rin=%g "
                "kappa=%.8e L=%.8e\n" % (a.n, a.tau, a.ratio, kappa, lum))
        for i in range(a.nfine):
            f.write("%.10e %.10e %.10e\n" % (r[i], rho[i], eint[i]))

    print("KAPPA=%.10e" % kappa)
    print("LSTAR=%.10e" % lum)
    print("X1MAX=%.10e" % rout)
    print("FINNER=%.10e" % (lum/(4.0*np.pi*rin**2)))
    print("TIN=%.6e TOUT=%.6e" % (np.interp(rin, r, temp), np.interp(rout, r, temp)))


if __name__ == "__main__":
    main()
