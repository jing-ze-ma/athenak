"""Fit the polynomial radial stretch (mesh/use_grid_stretch_r_poly, f_stretch_r_c1..c4)
to a stratified column, so that the grid carries a CONSTANT number of cells per pressure
scale height.

The map is
    r = r0 + (r1 - r0) u(xi),   u(xi) = xi + sum_k c_k xi^k (1 - xi),
so u(0) = 0 and u(1) = 1 whatever the coefficients (see src/coordinates/grid_stretch.hpp).
Constant cells per scale height means dr ~ H(r), i.e. uniform spacing in the cumulative
scale-height coordinate N_H(r) = int dr'/H(r'). So: build N_H from a column, invert it at
uniform N to get the target r(xi), and least-squares fit u_target in the basis xi^k(1-xi).

H is the DENSITY scale height, -(d ln rho/dr)^-1, which is what actually has to be
resolved -- a temperature-based H misses the ionization zones where rho falls fastest.

usage: fit_radial_stretch.py COLUMN r0 r1 [--nc 4] [--nx 128]
  COLUMN is the pgen's problem/column_dump: r[cm] p T rho ...
"""
import argparse
import numpy as np


ap = argparse.ArgumentParser()
ap.add_argument('column')
ap.add_argument('r0', type=float)
ap.add_argument('r1', type=float)
ap.add_argument('--nc', type=int, default=4, help='number of coefficients (must be 4)')
ap.add_argument('--nx', type=int, default=128, help='radial cells, for the report')
a = ap.parse_args()

c = np.loadtxt(a.column)
r, rho = c[:, 0], c[:, 3]
m = (r >= a.r0) & (r <= a.r1) & np.isfinite(rho) & (rho > 0)
r, rho = r[m], rho[m]
o = np.argsort(r)
r, rho = r[o], rho[o]

# the density scale height, smoothed a little: the raw d ln rho/dr is noisy on a fine grid
lnrho = np.log(rho)
dl = np.gradient(lnrho, r)
H = np.abs(1.0/np.where(np.abs(dl) > 0, dl, np.nan))
H = np.clip(H, 1e-6*(a.r1 - a.r0), a.r1 - a.r0)
k = max(3, len(r)//64) | 1
H = np.convolve(H, np.ones(k)/k, mode='same')
H[:k] = H[k]
H[-k:] = H[-k-1]

# cumulative scale heights, then invert at uniform N
N = np.concatenate([[0.0], np.cumsum(0.5*(1.0/H[1:] + 1.0/H[:-1])*np.diff(r))])
xi_t = np.linspace(0.0, 1.0, 4001)
r_t = np.interp(xi_t*N[-1], N, r)
u_t = (r_t - a.r0)/(a.r1 - a.r0)

# least squares for u(xi) - xi in the basis xi^k (1 - xi)
A = np.stack([xi_t**kk*(1.0 - xi_t) for kk in range(1, a.nc + 1)], axis=1)
coef, *_ = np.linalg.lstsq(A, u_t - xi_t, rcond=None)


def u_of(xi, cf):
    return xi + sum(cf[kk-1]*xi**kk*(1.0 - xi) for kk in range(1, len(cf) + 1))


# the map must be strictly increasing, or Mesh refuses it
xs = np.linspace(0, 1, 20001)
du = np.gradient(u_of(xs, coef), xs)
print('coefficients:  ' + '  '.join('f_stretch_r_c%d = %+.6f' % (i+1, v)
                                    for i, v in enumerate(coef)))
print('du/dxi over [0,1]: min %.4f  max %.4f  %s'
      % (du.min(), du.max(),
         'OK' if du.min() > 0 else 'NOT MONOTONIC -- Mesh will refuse it'))
print('the domain spans %.1f density scale heights, so N cells gives N/%.1f per scale'
      ' height at best; 10 per scale height wants nx1 = %d'
      % (N[-1], N[-1], int(np.ceil(10*N[-1]/10.0)*10)))

# what it buys, in cells per scale height at nx cells
for name, cf in (('uniform', np.zeros(a.nc)), ('fitted', coef)):
    xf = np.linspace(0, 1, a.nx + 1)
    rf = a.r0 + (a.r1 - a.r0)*u_of(xf, cf)
    rc = 0.5*(rf[1:] + rf[:-1])
    dr = np.diff(rf)
    Hc = np.interp(rc, r, H)
    cph = Hc/dr
    print('%-8s cells per scale height: min %.2f  median %.2f  max %.2f'
          % (name, cph.min(), np.median(cph), cph.max()))
