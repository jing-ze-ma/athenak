#!/usr/bin/env python3
"""rw-drag-1112: 4-mode acoustic omega per arm, N->inf limit (Richardson p=2 on the two finest N), offset vs the
reference root, and eps = the relative change of the momentum-exchange drag G = rho kappa_F F0/c (both the gas
momentum and the F equation) that moves the reference root to that limit (viper RW_P333_TERM_1112: eps = -3.15e-4
for the base arm at P 3.33 tau 10 c 10, both closures).
usage: RW_C=<c> python3 ana_drag.py <runs dir> <cases file>"""
import json
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "rw_amp1e6"))
import run_radwave as rw  # noqa: E402


def root_eps(bg, q, s):
    import radwave_disp as rd
    L = rd.matrix(bg, q)
    for (i, j) in ((4, 4), (4, 1), (1, 4), (1, 1)):
        L[i, j] *= (1.0 + s)
    return rd.pick_acoustic(bg, 1j * np.linalg.eigvals(L))


def eps_fit(c, dw):
    import radwave_disp as rd
    bg, kd = rw.background(c)
    q = rd.vetqs_coef(bg, rd.quad_code(4, 8, kd)) if c["clos"] == "vetf" else 0.0
    w0 = root_eps(bg, q, 0.0)
    d = (root_eps(bg, q, 1e-6) - w0) / 1e-6
    s = (np.conj(d) * dw).real / abs(d) ** 2
    return s, abs(d * s - dw) / max(abs(dw), 1e-300)


runs, cases = sys.argv[1], sys.argv[2]
G = {}
for ln in open(cases).read().splitlines():
    if not ln.strip():
        continue
    c = rw.parse_case(ln)
    p = os.path.join(runs, rw.case_tag(c), "results.json")
    arm = " ".join("%s=%s" % kv for kv in sorted(c["extra"].items()) if kv[0] in ("vimp",) or
                   kv[0] in ("rm_implicit_gas_src", "rm_implicit_vimp", "rm_time_scheme"))
    key = (c["clos"], c["prat"], arm)
    if not os.path.exists(p):
        G.setdefault(key, {})[c["N"]] = None
        continue
    r = json.load(open(p))
    G.setdefault(key, {})[c["N"]] = (complex(r["w_re"], r["w_im"]),
                                     complex(r.get("w4_re", np.nan), r.get("w4_im", np.nan)),
                                     r.get("w4_slow_amp", np.nan), r.get("w4_rms", np.nan), c, r["refs"][r["ref_key"]])
print("c = %s" % rw.C_LIGHT)
for key in sorted(G):
    d = G[key]
    print("== %s P %.4g %s" % (key[0], key[1], key[2]))
    for N in sorted(d):
        if d[N] is None:
            print("  N %4d MISSING (FATAL or not run)" % N)
            continue
        w1, w4, sa, rms, c, ref = d[N]
        print("  N %4d  1-mode %.6f%+.6fj  4-mode %.6f%+.6fj  slow amp %.2e  rms %.1e" % (
            N, w1.real, w1.imag, w4.real, w4.imag, sa, rms))
    Ns = [N for N in sorted(d) if d[N] is not None]
    if len(Ns) >= 2:
        a, b = Ns[-2], Ns[-1]
        wi = d[b][1] + (d[b][1] - d[a][1]) / ((b / a) ** 2 - 1.0)
        w1i = d[b][0] + (d[b][0] - d[a][0]) / ((b / a) ** 2 - 1.0)
        ref = complex(*d[b][5])
        s, res = eps_fit(d[b][4], wi - ref)
        print("  inf: 4-mode %.6f%+.6fj  1-mode %.6f%+.6fj  ref %.6f%+.6fj  dw4 = %+.2e%+.2ej  eps(G) = %+.3e (resid %.2f)" % (
            wi.real, wi.imag, w1i.real, w1i.imag, ref.real, ref.imag, (wi - ref).real, (wi - ref).imag, s, res))
