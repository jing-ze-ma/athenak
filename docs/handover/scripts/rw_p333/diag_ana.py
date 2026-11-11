#!/usr/bin/env python3
"""rw_p333-1112 closure check: code's VET perturbation vs the vetqs_d model at the fundamental.
usage: diag_ana.py RUNDIR N  -> per dump: t, R = E0 dD11_code/(qT T1_code), J1_code/(chi b j0 T1),
E1/T1 (code) vs the eigenvector, dD11/dD_cont(continuous angles)."""
import os, sys, glob
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rw_amp1e6"))
import radwave_disp as rd  # noqa: E402
import common  # noqa: E402

rdir = sys.argv[1]
bg = rd.Bg(10.0 / 3.0, 10.0, lam=1.0, c=10.0)
chi = bg.rho0 * bg.kappa
b = 4.0 * bg.arad
q = rd.quad_code(4, 8, (1, 0, 0))
j0, j1, j2, j3 = rd.jfun_quad(chi + 0j, bg.k, *q)
qT = chi * b * (j2 - j0 / 3.0)
w, vec = rd.all_roots(bg)["vetqs_d"]


def fh(x, a):
    return complex(((a - a.mean()) * np.exp(-1j * bg.k * x)).mean())


for ph in sorted(glob.glob(rdir + "/tab/rw.hydro_w.*.tab"))[1:]:
    n = ph.split(".")[-2]
    h = common.load_tab(ph)
    m = common.load_tab(rdir + "/tab/rw.m1.%s.tab" % n)
    v = common.load_tab(rdir + "/tab/rw.m1_vet.%s.tab" % n)
    x = np.asarray(h.x1)
    T = (2.0 / 3.0) * h.var("eint") / h.var("dens")
    t1 = fh(x, T)
    e1 = fh(x, m.var("m1_e"))
    r1 = fh(x, h.var("dens"))
    d1 = fh(x, v.var("m1_d11"))
    try:
        s = common.load_tab(rdir + "/tab/rw.m1_vsc.%s.tab" % n)
        jn = [k for k in ("m1_vj", "m1_jfs") if k in str(open(rdir + "/tab/rw.m1_vsc.%s.tab" % n).readlines()[1])]
        jj = fh(x, s.var(jn[0])) / (chi * b * j0 * t1) if jn else np.nan
    except Exception:  # noqa: BLE001
        jj = np.nan
    print("t=%.5f  R=E0 dD/(qT T1)=%.6f%+.6fj  |  J1/Jref=%s  |  (E1/T1)code/eig=%.6f%+.6fj  "
          "(T1/rho1)code/eig=%.6f%+.6fj  D0=%.9f" % (
              h.time, *(lambda z: (z.real, z.imag))(bg.e0 * d1 / (qT * t1)),
              ("%.6f%+.6fj" % (jj.real, jj.imag)) if jj == jj else "-",
              *(lambda z: (z.real, z.imag))((e1 / t1) / (vec[3] / vec[2])),
              *(lambda z: (z.real, z.imag))((t1 / r1) / (vec[2] / vec[0])),
              float(np.mean(v.var("m1_d11")))))
