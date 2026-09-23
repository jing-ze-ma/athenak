#!/usr/bin/env python3
"""Is the 'upw' instability of vimp_model.py part 2 an artefact of v0 = 0?  Add a uniform
background drift v0 > 0 so that the upwind cell is really the left one, with its own
first-order upwind advection: -v0 D_up on rho, v, T (the gas) and -(4/3) v0 D_up on E
(the upwinded enthalpy flux a E with a = v0 (1 + 1/3)).  F0 has no v0 term, as in the
code.  Output: max over m = 1..32 (N = 64), tau in 1,1e2,1e4,1e6 of rho(G) - 1."""
import math
import numpy as np
import vimp_model as vm
import radwave_disp as rd

nx = 64
dx = 1.0 / nx
print("v0/c_s  lh    cfl   P : imp rhs code codeimp   (upw; cen for reference)")
for v0f in (0.0, 0.1, 1.0):
    for lh in ("full", "Lh=0"):
        for cfl in (0.15, 0.3):
            for pr in (1.0, 10.0, 100.0):
                res = {}
                for enth in ("upw", "cen"):
                    w = {md: -1.0 for md in vm.MODES}
                    for tau in (1.0, 1e2, 1e4, 1e6):
                        bg = rd.Bg(pr, 1.0)
                        bg.kappa = tau / bg.rho0
                        cs = math.sqrt(bg.gamma * bg.t0)
                        dt = cfl * dx / (cs + v0f * cs)
                        v0 = v0f * cs
                        for m in range(1, nx // 2 + 1):
                            th = 2 * math.pi * m / nx
                            o = vm.DOps(bg, th, dx, enth)
                            du = (1 - np.exp(-1j * th)) / dx
                            for q in range(3):
                                o.Lh[q, q] += -v0 * du
                            o.Le[3, 3] += -(4.0 / 3.0) * v0 * du
                            if lh == "Lh=0":
                                o.Lh = np.zeros_like(o.Lh)
                            for md in vm.MODES:
                                lam = np.linalg.eigvals(vm.smap(o, md, dt))
                                r = max(abs(lam)) - 1
                                w[md] = max(w[md], r)
                    res[enth] = w
                print("%4.1f %-5s %.2f %5g : " % (v0f, lh, cfl, pr)
                      + " ".join("%.1e" % res["upw"][md] for md in vm.MODES)
                      + "   | cen "
                      + " ".join("%.1e" % res["cen"][md] for md in vm.MODES))
