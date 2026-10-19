#!/usr/bin/env python3
"""Option-3 radial grid with the EXISTING polynomial + plateau radial stretch
(mesh/use_grid_stretch_r_poly, f_stretch_r_c1..c8, f_stretch_r_p_amp/_xa/_xb/_w; see
src/coordinates/grid_stretch.hpp StretchRPoly, replicated here exactly).
Target: log spacing dr = f r outside a fine zone of nfine cells of dr_f centred on R_acc,
tanh-blended over wcells cells; the plateau term carries the step, c1..c8 fit the rest.
usage: grid_design.py [nx1] [dr_f] [nfine] [wcells]"""
import sys
import numpy as np

import os
# domain and R_acc: env GD_R0, GD_R1, GD_RACC (default RY Per); GD_ZOFF shifts the fine
# zone centre relative to R_acc (Rsun, < 0 = deeper)
R0 = float(os.environ.get('GD_R0', 2.835492977186198))
R1 = float(os.environ.get('GD_R1', 16.241))
RACC = float(os.environ.get('GD_RACC', 4.06))
ZOFF = float(os.environ.get('GD_ZOFF', 0.0))
nx1 = int(sys.argv[1]) if len(sys.argv) > 1 else 800
drf = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0e-4
nfine = int(sys.argv[3]) if len(sys.argv) > 3 else 200
wcells = float(sys.argv[4]) if len(sys.argv) > 4 else 12.0


def lcs(x):
    ax = np.abs(x)
    return ax + np.log1p(np.exp(-2*ax)) - np.log(2.0)


def G(xa, xb, w, xi):
    return 0.5*w*(lcs((xi-xa)/w) - lcs(-xa/w) - lcs((xi-xb)/w) + lcs(-xb/w))


def plateau_u(xi, p):
    return p[0]*(G(p[1], p[2], p[3], xi) - xi*G(p[1], p[2], p[3], 1.0))


def u_of(xi, c, p):
    xi = np.asarray(xi, dtype=float)
    u = xi.copy()
    xik = xi.copy()
    for k in range(8):
        u = u + c[k]*xik*(1-xi)
        xik = xik*xi
    if p[0] != 0:
        u = u + plateau_u(xi, p)
    return u


xi_f = np.linspace(0, 1, nx1 + 1)
w = wcells/nx1


def target_faces(f, ic):
    """faces from r = R0 with dr_i = drf + (f r_i - drf)(1 - g_i), g = plateau in index"""
    xa, xb = (ic - 0.5*nfine)/nx1, (ic + 0.5*nfine)/nx1
    xc = (np.arange(nx1) + 0.5)/nx1
    g = 0.5*(np.tanh((xc - xa)/w) - np.tanh((xc - xb)/w))
    r = np.empty(nx1 + 1)
    r[0] = R0
    for i in range(nx1):
        r[i+1] = r[i] + drf + (f*r[i] - drf)*(1.0 - g[i])
    return r, xa, xb


# solve f (reach R1) for a given zone centre ic, and ic (zone centre at R_acc) by bisection
def solve_f(ic):
    lo, hi = 1e-4, 1e-1
    for it in range(100):
        f = np.sqrt(lo*hi)
        r, xa, xb = target_faces(f, ic)
        if r[-1] < R1:
            lo = f
        else:
            hi = f
    return f


ilo, ihi = 0.5*nfine + 2*wcells, nx1 - 0.5*nfine - 2*wcells
for it in range(60):
    ic = 0.5*(ilo + ihi)
    f = solve_f(ic)
    r, xa, xb = target_faces(f, ic)
    if np.interp(ic, np.arange(nx1 + 1), r) < RACC + ZOFF:
        ilo = ic
    else:
        ihi = ic
r, xa, xb = target_faces(f, ic)
ut = (r - R0)/(R1 - R0)
# plateau amplitude from the step in du/dxi between the zone and just outside it
out_slope = f*RACC*nx1/(R1 - R0)
in_slope = drf*nx1/(R1 - R0)
a = in_slope - out_slope
p = [a, xa, xb, w]
resid = ut - xi_f - plateau_u(xi_f, p)
A = np.array([xi_f**(k+1)*(1 - xi_f) for k in range(8)]).T
c, *_ = np.linalg.lstsq(A, resid, rcond=None)
rf = R0 + (R1 - R0)*u_of(xi_f, c, p)
dr = np.diff(rf)
ok = np.all(dr > 0)
print(f"nx1 {nx1}, f (dr/r outside) {f:.4e}, zone centre index {ic:.1f}, fold-free {ok}")
if ok:
    ratio = np.maximum(dr[1:]/dr[:-1], dr[:-1]/dr[1:])
    iR = np.argmin(np.abs(rf - RACC))
    fine = dr < 1.1*drf
    print(f"  dr_min {dr.min():.3e}  dr(r_in) {dr[0]:.4e}  dr(r_out) {dr[-1]:.4e}  max "
          f"neighbour ratio {ratio.max():.4f}; max |r - r_target| {np.abs(rf - r).max():.2e}")
    print(f"  cells with dr < 1.1 dr_f: {fine.sum()}, r {rf[:-1][fine].min():.4f}.."
          f"{rf[1:][fine].max():.4f}; face nearest R_acc: index {iR}, r {rf[iR]:.7f}")
    for x in (RACC - 0.5, RACC - 0.3, RACC - 0.15, RACC + 0.1, RACC + 0.15,
              RACC + 0.5, R1 - 0.01):
        print(f"  dr at r {x:.2f}: {np.interp(x, rf[:-1], dr):.4e}")
    print("use_grid_stretch_r_poly = true")
    for k in range(8):
        print(f"f_stretch_r_c{k+1} = {c[k]:.17g}")
    print(f"f_stretch_r_p_amp = {p[0]:.17g}\nf_stretch_r_p_xa = {p[1]:.17g}\n"
          f"f_stretch_r_p_xb = {p[2]:.17g}\nf_stretch_r_p_w = {p[3]:.17g}")
    print(f"r_acc (face) = {rf[iR]:.15g}")
