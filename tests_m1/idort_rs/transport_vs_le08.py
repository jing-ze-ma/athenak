#!/usr/bin/env python3
"""transport_vs_le08.py -- MODEL difference: the IDORT transport reference
(Rad_shock_sol/M<M>/*_4096.dat of the IDORT Zenodo package 15032138: S_N grey transport,
Ferguson+17 type, nondimensional rho0 = T0 = 1, p = rho T / gamma, E = T_r^4) against the
LE08 grey non-equilibrium DIFFUSION solution (le08_shock.py, nominal P0 = 1e-4,
sigma_a = 1e6, kappa = 1), with the same metric definitions as rs_eval.py.  A VET/M1 code
that is accurate converges to the TRANSPORT answer, so these numbers are the floor below
which an rs_eval.py difference against LE08 is not an error.
usage: python transport_vs_le08.py [<Rad_shock_sol dir>]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import le08_shock as L  # noqa: E402

REF = ("/raven/ptmp/jinma/suite_1011/ext/idort_mini/x/IDORT_mini/shock_1d/Rad_shock_sol"
       if len(sys.argv) < 2 else sys.argv[1])
GAM = 5.0 / 3.0


def main():
    for m0, lab in ((1.2, "M1p2"), (2.0, "M2"), (3.0, "M3")):
        r = {k: np.loadtxt(os.path.join(REF, lab, "%s_4096.dat" % k))
             for k in ("x", "rho", "p", "E", "u")}
        x, rho = r["x"], r["rho"]
        tg = GAM * r["p"] / rho
        tr = np.maximum(r["E"], 0.0) ** 0.25
        mach = np.abs(r["u"]) / np.sqrt(tg)
        if mach[0] < 1.0:           # flow from the right: mirror
            x, rho, tg, tr, mach = -x[::-1], rho[::-1], tg[::-1], tr[::-1], mach[::-1]
        k = np.nonzero((mach[:-1] >= 1.0) & (mach[1:] < 1.0))[0][0]
        w = (mach[k] - 1.0) / (mach[k] - mach[k + 1])
        xc = x[k] + w * (x[k + 1] - x[k])
        s = L.Solution(L.LE08(m0))
        mt = s.metrics()
        ref = s.sample(x - xc)
        dx = x[1] - x[0]
        away = np.abs(x - xc) > 3 * dx
        l1 = [np.sum(np.abs(q[away] - rr[away])) / np.sum(np.abs(rr[away]))
              for q, rr in ((rho, ref[0]), (tg, ref[2]), (tr, ref[3]))]
        pre = x < xc - 3 * dx
        l1p = np.sum(np.abs(tg[pre] - ref[2][pre])) / np.sum(np.abs(ref[2][pre] - 1.0))
        lp = {}
        for f in (0.1, 0.01):
            lev = 1.0 + f * (mt["Tp"] - 1.0)
            j = np.nonzero((x < xc) & (tg >= lev))[0][0]
            wl = (lev - tg[j - 1]) / (tg[j] - tg[j - 1])
            lp[f] = (xc - (x[j - 1] + wl * dx)) / mt["Lpre_%g" % f] - 1.0
        trs = np.interp(xc, x, tr) / mt["Trs"] - 1.0
        print("TRANSPORT(IDORT) vs LE08 M%-4g L1 rho %.4f Tg %.4f Tr %.4f pre %.4f | "
              "Lpre0.1 %+6.2f%% Lpre0.01 %+6.2f%% | Trs %+6.3f%% | Tmax %.4f vs %.4f (%+.2f%%)"
              % (m0, l1[0], l1[1], l1[2], l1p, 100 * lp[0.1], 100 * lp[0.01], 100 * trs,
                 tg.max(), mt["Tmax"], 100 * (tg.max() / mt["Tmax"] - 1)))


if __name__ == "__main__":
    main()
