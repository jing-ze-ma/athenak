#!/usr/bin/env python3
"""Budget and per-cell comparison of two restarts on the SAME grid (round trip) or global budgets of any two
restarts (different grids): mass, U = E - rho Phi (thermal + kinetic), E_tot (hydro E incl. rho Phi), KE, angular
momenta int r m_theta, int r m_phi, E_rad, int F_r.  usage: roundtrip_cmp.py A.rst B.rst"""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/he_giant_1006/remap')
import he_remap_rst_w256 as H
from phi_code import phi_cells


def load(fn):
    r = H.read_rst(fn)
    G = H.gather(r)
    ng, n1 = r.ng, r.bind[1]
    p = [float(H.get_param(r.text, 'mesh', k) or 0.0) for k in H.SKEYS]
    f = H.faces_of(p, n1, ng, r.msize[0], r.msize[3])
    x1v, ph = phi_cells(f, ng, H.get_param(r.text, 'problem', 'he_ic_file'),
                        float(H.get_param(r.text, 'problem', 'he_gm')),
                        int(H.get_param(r.text, 'problem', 'he_nfine')))
    a = slice(ng, ng + n1)
    dV = np.diff(f**3 / 3.0)[a]
    u, w = G['hyd'][..., a], G['m1'][..., a]
    ke = 0.5 * (u[1]**2 + u[2]**2 + u[3]**2) / u[0]
    x = x1v[a]
    q = {'mass': u[0], 'E_tot': u[4], 'U': u[4] - u[0] * ph[a], 'KE': ke, 'L_th': u[2] * x,
         'L_ph': u[3] * x, 'E_rad': w[0], 'F_r': w[1]}
    return r, {k: float((v * dV).sum()) for k, v in q.items()}, u, w, x


ra, ta, ua, wa, xa = load(sys.argv[1])
rb, tb, ub, wb, xb = load(sys.argv[2])
print('t %.6g / %.6g, nx1 %d / %d' % (ra.time, rb.time, ra.bind[1], rb.bind[1]))
for k in ta:
    print('  %-6s %+.6e  %+.6e  rel diff %+.3e' % (k, ta[k], tb[k], (tb[k] - ta[k]) / max(abs(ta[k]), 1e-300)))
if ua.shape == ub.shape:
    for nm, A, B in (('rho', ua[0], ub[0]), ('E', ua[4], ub[4]), ('m_r', ua[1], ub[1]), ('E_rad', wa[0], wb[0])):
        d = np.abs(B - A) / np.maximum(np.abs(A), 1e-300)
        dm = d.max(axis=(0, 1))
        i = int(np.argmax(dm))
        sel = xa < 60 * 6.957e10
        print('  per-cell %-5s max rel %.3e at r %.2f Rsun; median %.3e; max over r < 60 Rsun %.3e' % (
            nm, d.max(), xa[i] / 6.957e10, np.median(d), dm[sel].max()))
