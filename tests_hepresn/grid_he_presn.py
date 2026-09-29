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
ip = lambda q: np.exp(np.interp(rc, r, np.log(q)))
Hc, Hpc, tauc = ip(Hs), ip(Hp), np.interp(rc, r, tau)
cs = np.sqrt(5.0 / 3.0 * ip(pg) / ip(rho))
print('c1..c8 =', ' '.join('%.10e' % v for v in c))
print('dx: min %.3e (r/R %.4f) max %.3e (r/R %.4f); adjacent ratio max %.3f'
      % (dx.min(), rc[dx.argmin()] / R, dx.max(), rc[dx.argmax()] / R,
         np.max(np.maximum(dx[1:] / dx[:-1], dx[:-1] / dx[1:]))))
for nm, z_ in (('tau>100', tauc > 100), ('tau 1-100', (tauc >= 1) & (tauc <= 100)),
               ('tau<1', tauc < 1)):
    print('  %-9s cells %3d  min cells/H_gas %.1f (r/R %.4f)  min cells/Hp %.1f'
          % (nm, z_.sum(), (Hc / dx)[z_].min(), rc[z_][(Hc / dx)[z_].argmin()] / R,
             (Hpc / dx)[z_].min()))
dtp = (dx / cs).min()
dxu = (r1 - r0) / nx
print('dt proxy min dx/c_gas: %.3f s at r/R %.4f (uniform grid: %.3f s at r/R %.4f)'
      % (dtp, rc[(dx / cs).argmin()] / R, (dxu / cs).min(), rc[(dxu / cs).argmin()] / R))
np.save(out, c)
