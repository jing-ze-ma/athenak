#!/usr/bin/env python3
"""rw_edd: the analytic Eddington (D = I/3) acoustic root of the linearised rad_m1 system (radwave_disp.matrix,
lab-frame F with the (4/3) E0 v inertia/drag terms = 'edd'; 'eddc' = comoving-F0 form without -(4/3)E0 dv/dt),
c = 10, gamma 5/3, kappa_P = kappa_E = kappa_F, no scattering, for the 18 gated cells; plus the two limit checks
(tau -> 0: adiabatic gas speed; tau -> inf: mixture speed / sqrt(1 + (4/3) E0/(rho c^2)))."""
import math
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rw_amp1e6"))
import radwave_disp as rd  # noqa: E402
CELLS = [(1 / 300, t) for t in (0.1, 10, 1000)] + [(P, t) for P in (1 / 30, 1 / 3)
         for t in (0.01, 0.1, 1, 10, 100, 1000)] + [(10 / 3, t) for t in (0.1, 10, 1000)]
print("%8s %7s | %14s %13s | %14s %13s" % ("P", "tau", "edd w_re", "w_im", "eddc w_re", "w_im"))
for P, t in CELLS:
    bg = rd.Bg(P, t, c=10.0)
    a = rd.moment_root(bg)[0]; b = rd.moment_root(bg, 0.0, "f0")[0]
    print("%8.4g %7g | %14.8f %13.6e | %14.8f %13.6e" % (P, t, a.real, a.imag, b.real, b.imag))
for P in (1 / 300, 10 / 3):
    bg = rd.Bg(P, 1e-6, c=10.0); w = rd.moment_root(bg)[0]
    print("limit tau->0  P %.4g: c_ph %.8f vs sqrt(gamma T0) %.8f" % (P, w.real / bg.k, bg.speeds()[1]))
    bg = rd.Bg(P, 1e6, c=10.0); w = rd.moment_root(bg)[0]
    print("limit tau->inf P %.4g: c_ph %.8f vs c_mix/sqrt(1+4E0/3c^2) %.8f" % (
        P, w.real / bg.k, bg.mixture_cs() / math.sqrt(1 + 4 * bg.e0 / (3 * bg.c ** 2))))
