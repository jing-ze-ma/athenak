#!/usr/bin/env python3
"""Write the analytic grey Eddington column as a box_convection `problem/ic_profile`.

A static, grey, plane-parallel atmosphere in radiative equilibrium carrying a flux
F = sigma Teff^4 through a CONSTANT opacity kappa has, exactly,

    T(tau)^4 = (3/4) Teff^4 (tau + 2/3)        (the Eddington closure)
    p(tau)   = g tau / kappa                   (hydrostatic, since dp = g dtau/kappa)

so the whole column follows from tau alone.  Starting a run FROM this profile is what
turns the grey-atmosphere problem into a CONSERVATION test: the state is already the
steady solution the scheme is supposed to hold, so any energy the column gains is the
scheme's own book-keeping error and nothing else.  That is how the tau-blend handover
bug was localized -- see inputs/hydro/box_convection_grey_cons.athinput.

Writes "z rho eint", z increasing, covering the mesh and generously beyond it so the
ghosts are covered too (box_convection.cpp refuses a profile that is not).

    python3 tools/box_convection/mk_eddington_column.py edd_column.txt
"""
import sys

import numpy as np

K_B = 1.380649e-16
M_U = 1.66053906660e-24

MU = 0.6                 # mean molecular weight; matches <units>/mu and problem/mu
GAMMA = 5.0 / 3.0
G0 = 1.0e4               # cm/s^2, problem/g0
KAPPA = 0.4              # cm^2/g, the constant table inputs/hydro/kappa_const.txt
TEFF = 5000.0            # K; F = sigma Teff^4 = <hydro>/rad_flux_inner
TAU_BOT = 1000.0         # the optical depth that sits at x1min
RGAS = K_B / (MU * M_U)


def temperature(tau):
    return TEFF * (0.75 * (tau + 2.0 / 3.0)) ** 0.25


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "edd_column.txt"
    # tau far wider than the mesh at both ends, so the ghosts are covered whatever the
    # resolution; the run only ever interpolates inside it
    lt = np.linspace(np.log(1.0e-5), np.log(1.0e5), 40001)
    tau = np.exp(lt)
    temp = temperature(tau)
    # dz/dln tau = -R T/g, integrated upward, then shifted so z = 0 at TAU_BOT
    dzdl = -(RGAS * temp / G0)
    z = np.concatenate([[0.0], np.cumsum(0.5 * (dzdl[1:] + dzdl[:-1]) * np.diff(lt))])
    z -= np.interp(np.log(TAU_BOT), lt, z)
    pres = G0 * tau / KAPPA
    rho = pres / (RGAS * temp)
    eint = pres / (GAMMA - 1.0)
    top = np.interp(np.log(0.002), lt, z)
    keep = (z > -1.0e8) & (z < top + 2.0e8)
    np.savetxt(out, np.c_[z[keep][::-1], rho[keep][::-1], eint[keep][::-1]],
               fmt="%.10e",
               header="z rho eint -- the analytic grey Eddington column: "
                      "kappa = %g cm^2/g, g = %g cm/s^2, mu = %g, Teff = %g K, "
                      "tau at x1min = %g" % (KAPPA, G0, MU, TEFF, TAU_BOT))
    print("%s: %d nodes, z = %.6e .. %.6e; put x1min = 0 and x1max = %.8e "
          "(tau = 0.002)" % (out, keep.sum(), z[keep].min(), z[keep].max(), top))


if __name__ == "__main__":
    main()
