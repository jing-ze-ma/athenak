#!/usr/bin/env python3
"""Contraction factor of the Picard pass that re-evaluates v from the force of the
previous pass (enthalpy implicit via lagged v): v_{m+1} = v(Y(v_m)),
Y(v) = (I - g dt Lr)^{-1} (rhs + g dt Le v).  |factor| < 1 needed; continuum symbol,
k up to Nyquist of N = 64, dt = cfl dx / c_s,gas, g = 1 - 1/sqrt 2 (and g = 1: BE)."""
import math
import sys
import numpy as np
sys.path.insert(0, "/viper/u2/jinma/ATHENAK/athenak/tests_m1/runs_3t_time2")
sys.path.insert(0, "/viper/u2/jinma/ATHENAK/athenak/tests_m1/runs_3r_radwave")
import time2_model as tm  # noqa: E402
import radwave_disp as rd  # noqa: E402

for gname, g in (("hesdirk2 g=0.293", tm.GAM), ("BE g=1", 1.0)):
    print("==", gname)
    for cfl in (0.15, 0.3, 0.8):
        worst = {}
        for pr in (0.1, 1.0, 10.0, 100.0):
            w = 0.0
            for tau in (1e-1, 1e1, 1e2, 1e3, 1e4, 1e6):
                for m in (1, 2, 4, 8, 16, 32):
                    bg = rd.Bg(pr, tau, lam=1.0 / m)
                    bg.kappa = tau / bg.rho0
                    o = tm.Ops(bg)
                    dt = cfl * (1.0 / 64) / math.sqrt(bg.gamma * bg.t0)
                    A = np.eye(5) - g * dt * o.Lr
                    col = g * dt * o.Le[:, 1]
                    f = abs(np.linalg.solve(A, col)[1])
                    w = max(w, f)
            worst[pr] = w
        print("  cfl %.2f  " % cfl + "  ".join("P=%g: %.3g" % (p, worst[p]) for p in worst))
