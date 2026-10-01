#!/usr/bin/env python3
"""Radial grid design for he_star_m1 (poly stretch, mesh/f_stretch_r_c1..8).
usage: grid_he_presn.py STRUCT_NPZ NX1 RTOP NC OUT_NPY
The resolution scale is the TRUE local scale height of the gas state,
  H = min(|dln Pg/dr|^-1, |dln rho/dr|^-1, |dln T/dr|^-1),
NOT Pg/(rho g): with Gamma = kappa F/(c g) = 0.8..1.08 the gas feels g (1 - Gamma) and
Pg/(rho g) underestimates the gas scale height by 5x..400x.  Target width
dx_t = clip(H/NC, dx_lo, dx_hi); the map u(xi) of the target is fitted by least squares
with the 8-coefficient poly stretch; prints cells per H, per total Hp, and a dt proxy
dx/c_gas (c_gas^2 = 5/3 Pg/rho) against the uniform grid."""
import sys
import numpy as np
R = 2.3717e11
GM = 4.18143e26
z = np.load(sys.argv[1])
nx, rtop, NC, out = int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
r, rho, T, pg, tau = [z[k] for k in ('r', 'rho', 'T', 'pg', 'tau')]
A = 7.5657332503e-15
P = pg + A * T**4 / 3
g = GM / r**2
Hp = P / (rho * g)


def hl(q):
    return 1.0 / np.maximum(np.abs(np.gradient(np.log(q), r)), 1e-40)


H = np.minimum(np.minimum(hl(pg), hl(rho)), hl(T))
# smooth H over +-1e9 cm (the scale of a cell) so the target is not noisy
w = int(1e9 / np.median(np.diff(r)))
Hs = np.exp(np.convolve(np.log(H), np.ones(2 * w + 1) / (2 * w + 1), mode='same'))
r0, r1 = 0.5 * R, rtop
m = (r >= r0) & (r <= r1)
rr, Ht = r[m], Hs[m]
dxlo, dxhi = 5.0e7, float(__import__('os').environ.get('DXHI', '4e8'))
dxtop = float(__import__('os').environ.get('DXTOP', '1.5e8'))
dxt = np.clip(Ht / NC, np.where(tau[m] < 1.0, dxtop, dxlo), dxhi)
# FeCZ cap (env FECZ_N cells across 0.636-0.933 R) and total-P Hp rule (same NC)
fz = (rr >= 0.636 * R) & (rr <= 0.933 * R)
nfz = float(__import__('os').environ.get('FECZ_N', '0'))
if nfz > 0:
    dxt[fz] = np.minimum(dxt[fz], 0.297 * R / nfz)
dxt = np.minimum(dxt, np.exp(np.interp(rr, r, np.log(Hp))) / NC)
# smooth stretching: limit the log-slope of dx to MAXR per cell (two sweeps)
mr = float(__import__('os').environ.get('MAXR', '1.04'))
for _ in range(3):
    for k in range(1, len(dxt)):
        lim = dxt[k-1] * mr ** ((rr[k] - rr[k-1]) / dxt[k-1])
        dxt[k] = min(dxt[k], lim)
    for k in range(len(dxt) - 2, -1, -1):
        lim = dxt[k+1] * mr ** ((rr[k+1] - rr[k]) / dxt[k+1])
        dxt[k] = min(dxt[k], lim)
# cumulative cell count of the target, rescaled to nx cells
cnt = np.concatenate([[0], np.cumsum(np.diff(rr) / (0.5 * (dxt[1:] + dxt[:-1])))])
print('target cells %.1f (requested %d)' % (cnt[-1], nx))
xi_of_r = cnt / cnt[-1]
u_of_r = (rr - r0) / (r1 - r0)
xs = np.linspace(0, 1, 2001)
us = np.interp(xs, xi_of_r, u_of_r)
B = np.stack([xs**k * (1 - xs) for k in range(1, 9)], axis=1)
c, *_ = np.linalg.lstsq(B[1:-1], (us - xs)[1:-1], rcond=None)
xe = np.arange(nx + 1) / nx
ue = xe + (np.stack([xe**k * (1 - xe) for k in range(1, 9)], axis=1) @ c)
re_ = r0 + (r1 - r0) * ue
dx = np.diff(re_)
assert (dx > 0).all(), 'fold-over'
rc = 0.5 * (re_[1:] + re_[:-1])


def ip(q):
    return np.exp(np.interp(rc, r, np.log(q)))


Hc, Hpc, tauc = ip(Hs), ip(Hp), np.interp(rc, r, tau)
cs = np.sqrt(5.0 / 3.0 * ip(pg) / ip(rho))
print('c1..c8 =', ' '.join('%.10e' % v for v in c))
rat = np.maximum(dx[1:] / dx[:-1], dx[:-1] / dx[1:])
print('adjacent ratio max at r/R %.4f; ratio > 1.05 in %d faces'
      % (rc[1:][rat.argmax()] / R, (rat > 1.05).sum()))
print('dx: min %.3e (r/R %.4f) max %.3e (r/R %.4f); adjacent ratio max %.3f'
      % (dx.min(), rc[dx.argmin()] / R, dx.max(), rc[dx.argmax()] / R,
         np.max(np.maximum(dx[1:] / dx[:-1], dx[:-1] / dx[1:]))))
for nm, z_ in (('tau>100', tauc > 100), ('tau 1-100', (tauc >= 1) & (tauc <= 100)),
               ('tau<1', tauc < 1)):
    print('  %-9s cells %3d  min cells/H_gas %.1f (r/R %.4f)  min cells/Hp %.1f'
          % (nm, z_.sum(), (Hc / dx)[z_].min(), rc[z_][(Hc / dx)[z_].argmin()] / R,
             (Hpc / dx)[z_].min()))
fzc = (rc >= 0.636 * R) & (rc <= 0.933 * R)
print('  FeCZ cells %d; rad zone 0.5-0.636 cells %d'
      % (fzc.sum(), (rc < 0.636 * R).sum()))
for xx in (0.5, 0.55, 0.6, 0.636, 0.7, 0.8, 0.9, 0.933, 0.96, 0.99, 1.0, 1.02):
    j = np.argmin(abs(rc - xx * R))
    print('    r/R %.3f dx %.2e' % (xx, dx[j]))
dtp = (dx / cs).min()
dxu = (r1 - r0) / nx
print('dt proxy min dx/c_gas: %.3f s at r/R %.4f (uniform grid: %.3f s at r/R %.4f)'
      % (dtp, rc[(dx / cs).argmin()] / R, (dxu / cs).min(), rc[(dxu / cs).argmin()] / R))
np.save(out, c)
