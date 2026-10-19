"""gridfit2.py (s3grid-1017): extended fine radial zone for S3 (AG Car A 268x32x32 test grid).
Target faces: the S3 active faces up to r_keep, then a smooth (tanh in ln r) change of dr/r to QF, kept to R_FINE, then geometric
growth (ratio <= GROW per cell) to x1max; least-squares fit of StretchRPoly (8 poly + 2 bumps + 2 plateaus, grid_stretch.hpp).
usage: gridfit2.py QF R_FINE [GROW] [r_keep]   (R in units of R_ph = 2.7025e13; QF e.g. 0.003)"""
import sys
import numpy as np
from scipy.optimize import least_squares

RPH = 2.7025e13
g = np.array([[float(v) for v in ln.split()] for ln in open('/raven/ptmp/jinma/s3grid_1017/grid/S3_x1grid.txt')
              if not ln.startswith('#') and len(ln.split()) == 4])
fo = np.concatenate([g[3:271, 1], [g[270, 1] + g[270, 3]]])     # 269 active faces of S3
r0, r1 = fo[0], fo[-1]
QF = float(sys.argv[1]); RF = float(sys.argv[2]) * RPH
GROW = float(sys.argv[3]) if len(sys.argv) > 3 else 1.04
RK = float(sys.argv[4]) * RPH if len(sys.argv) > 4 else 1.08 * RPH
# 1) keep the S3 faces to RK
fk = fo[fo <= RK]
f = list(fk)
q0 = (fk[-1] - fk[-2]) / fk[-2]
# 2) dr/r from q0 to QF over a ln r width of 0.03 (smooth), constant QF to RF
while f[-1] < RF:
    x = np.log(f[-1] / fk[-1]) / 0.03
    q = q0 + (QF - q0) * np.tanh(x)
    f.append(f[-1] * (1 + q))
# 3) geometric growth of dr to r1
dr = f[-1] - f[-2]
while f[-1] < r1:
    dr *= GROW
    f.append(f[-1] + dr)
f = np.array(f)
# stretch the tail (from RF) so that the last face is r1 exactly
it = int(np.searchsorted(f, RF))
f[it:] = f[it] + (f[it:] - f[it]) * (r1 - f[it]) / (f[-1] - f[it])
nx1 = len(f) - 1
xi = np.linspace(0, 1, nx1 + 1)
u = (f - r0) / (r1 - r0)


def lc(x):
    ax = np.abs(x)
    return ax + np.log1p(np.exp(-2 * ax)) - np.log(2)


def G(xa, xb, w, x):
    return 0.5 * w * (lc((x - xa) / w) - lc(-xa / w) - lc((x - xb) / w) + lc(-xb / w))


def model(p, x):
    c = p[:8]
    uu = x.copy(); xk = x.copy()
    for k in range(8):
        uu = uu + c[k] * xk * (1 - x); xk = xk * x
    for b in range(2):
        a, xb, w = p[8 + 3 * b:11 + 3 * b]
        uu = uu + a * w * (np.tanh((x - xb) / w) - (1 - x) * np.tanh(-xb / w) - x * np.tanh((1 - xb) / w))
    for q in range(2):
        a, xa, xb, w = p[14 + 4 * q:18 + 4 * q]
        uu = uu + a * (G(xa, xb, w, x) - x * G(xa, xb, w, 1.0))
    return uu


def dudx(p, x, h=1e-7):
    return (model(p, x + h) - model(p, x - h)) / (2 * h)


xm = 0.5 * (xi[1:] + xi[:-1]); dt_ = np.diff(u) / np.diff(xi)
W = np.ones_like(xm); W[f[1:] < 1.3 * RPH] = 3.0      # weight the kept/fine part


def resid(p):
    d = dudx(p, xm)
    return np.r_[W * np.log(np.maximum(d, 1e-6) / dt_), (model(p, xi) - u) * 30]


# initial guess: the S3 coefficients (their plateaus), bumps off
p_s3 = [-6.847926247688e-01, 3.310337586667e+00, -8.118982563059e+00, -6.896582429047e-01, 1.322995653966e-01,
        1.145439618682e+01, 1.761488973219e+00, -1.681689800188e+01]
lo = np.r_[-80 * np.ones(8), -0.5, 0.0, 0.003, -0.5, 0.0, 0.003, -0.99, 0.0, 0.02, 0.003, -0.99, 0.0, 0.02, 0.003]
hi = np.r_[80 * np.ones(8), 0.5, 1.0, 0.3, 0.5, 1.0, 0.3, 3.0, 0.98, 1.0, 0.3, 3.0, 0.98, 1.0, 0.3]
best = None
xf_hi = float(np.interp(RF, f, xi))
for xa1, xb1 in ((0.39, 0.97), (0.40, xf_hi), (0.35, xf_hi + 0.02)):
    for a2 in (-0.97, -0.9):
        p0 = np.r_[p_s3, 0, 0.5, 0.05, 0, 0.8, 0.05, -0.6, xa1, xb1, 0.027, a2, 0.02, 0.24, 0.085]
        p0 = np.clip(p0, lo + 1e-9, hi - 1e-9)
        res = least_squares(resid, p0, bounds=(lo, hi), max_nfev=4000)
        if best is None or res.cost < best.cost:
            best = res
p = best.x
rf = r0 + (r1 - r0) * model(p, xi)
drf, drt = np.diff(rf), np.diff(f)
xx = np.linspace(0, 1, 4097)
print('# QF %.4f R_FINE %.3f GROW %.3f r_keep %.3f: nx1 %d (S3 268, +%d); min du/dxi %.3e; max |dr/dr_tgt-1| %.3f'
      ' (r<1.75 R_ph: %.3f)' % (QF, RF / RPH, GROW, RK / RPH, nx1, nx1 - 268, dudx(p, xx).min(), np.abs(drf / drt - 1).max(),
                                np.abs(drf / drt - 1)[f[1:] < 1.75 * RPH].max()))
print('nx1 = %d' % nx1)
for k in range(8):
    print('f_stretch_r_c%d = %.12e' % (k + 1, p[k]))
for b in range(2):
    print('f_stretch_r_b%d_amp = %.12e\nf_stretch_r_b%d_x = %.12e\nf_stretch_r_b%d_w = %.12e' % (b + 1, p[8 + 3 * b], b + 1,
                                                                                           p[9 + 3 * b], b + 1, p[10 + 3 * b]))
for q, nm in enumerate(('p', 'p2')):
    print('f_stretch_r_%s_amp = %.12e\nf_stretch_r_%s_xa = %.12e\nf_stretch_r_%s_xb = %.12e\nf_stretch_r_%s_w = %.12e' %
          (nm, p[14 + 4 * q], nm, p[15 + 4 * q], nm, p[16 + 4 * q], nm, p[17 + 4 * q]))
rc = 0.5 * (rf[1:] + rf[:-1])
gc = 0.5 * (g[3:271, 1] + g[3:271, 1] + g[3:271, 3])
print('# r/R_ph  dr/r new (target) | S3')
for a in (0.3, 0.6, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.4, 2.9):
    i = int(np.argmin(abs(rc - a * RPH))); j = int(np.argmin(abs(gc - a * RPH)))
    print('#  %.2f  %.4f%% (%.4f%%) | %.4f%%' % (a, 100 * drf[i] / rc[i], 100 * drt[i] / rc[i], 100 * g[3 + j, 3] / gc[j]))
np.save('/raven/ptmp/jinma/s3grid_1017/grid/faces_q%.4f_r%.2f.npy' % (QF, RF / RPH), rf)
