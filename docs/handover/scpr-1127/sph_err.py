#!/usr/bin/env python3
"""SCPR: J of a homogeneous sphere (radius 0.6, S 1, sigma SIG, default 5) against exact transport.
usage: [SIG=50] [BINPY=<athenak>/vis/python] sph_err.py BINFILE VAR [gl NMU NPHI] [sp]
  VAR: sc_rad_J (PR #798) or vet_vj (ours).  sp: spherical-polar output (x1 = r).
Prints (a) interior cells r < 0.36 (PR's r_mask 0.6 R): rel L2 / Linf of J vs the exact continuous J,
and with 'gl' also vs the same-quadrature J of our GL product set (pole on x1, Cartesian only);
(b) the exterior profile: mean J/J_exact in shells r/R = 1.00-1.05, ..., up to the domain edge."""
import os, sys, numpy as np
sys.path.insert(0, os.environ.get('BINPY', '.'))
import bin_convert as bc
R, S = 0.6, 1.0; SIG = float(os.environ.get('SIG', '5.0'))
a = sys.argv; f = bc.read_binary(a[1]); var = a[2]; sp = 'sp' in a
assert f['n_mbs'] == 1, 'one MeshBlock expected'
if var not in f['mb_data']:  # the PR names its J 'rad_J'
    var = {'sc_rad_J': 'rad_J'}.get(var, var)
if var not in f['mb_data']: sys.exit(f"{var} not in {list(f['mb_data'])}")
J = np.asarray(f['mb_data'][var][0], dtype=float)
def cc(c):  # cell centres of the single MeshBlock (uniform grid)
    n = J.shape[{'1': 2, '2': 1, '3': 0}[c]]; lo, hi = f['x' + c + 'min'], f['x' + c + 'max']
    return lo + (np.arange(n) + 0.5)*(hi - lo)/n
x, y, z = cc('1'), cc('2'), cc('3')
Z, Y, X = np.meshgrid(z, y, x, indexing='ij')
if sp:
    r = X; X, Y, Z = r*np.sin(Y)*np.cos(Z), r*np.sin(Y)*np.sin(Z), r*np.cos(Y)
else:
    r = np.sqrt(X*X + Y*Y + Z*Z)
mu, wq = np.polynomial.legendre.leggauss(800)
def ipath(rr, m):  # optical path (unit sigma) inside the sphere behind a point at rr, direction cosine m
    disc = R*R - rr*rr*(1.0 - m*m)
    ok = disc > 0.0; sq = np.sqrt(np.where(ok, disc, 0.0))
    return np.where(rr < R, rr*m + sq, np.where(ok & (m > 0.0), 2.0*sq, 0.0))
def jex(rv):  # outside: nodes mapped onto the cone [mu_c, 1] that sees the sphere
    rv = np.asarray(rv, dtype=float); out = rv >= R
    mc = np.where(out, np.sqrt(np.clip(1.0 - (R/np.maximum(rv, 1e-30))**2, 0.0, 1.0)), -1.0)
    mm = 0.5*(1.0 - mc)[:, None]*mu[None, :] + 0.5*(1.0 + mc)[:, None]
    return 0.5*S*0.5*(1.0 - mc)*np.sum(wq[None, :]*(1.0 - np.exp(-SIG*ipath(rv[:, None], mm))), axis=1)
def err(p, q): return np.sqrt(np.sum((p - q)**2)/np.sum(q*q)), np.max(np.abs(p - q)/q)
tag = f"{os.path.basename(a[1])} {'sp' if sp else 'cart'} sigma={SIG:g} grid={J.shape[::-1]}"
m = r < 0.6*R
if m.any():
    rv = r[m]; je = jex(rv); e2, ei = err(J[m], je)
    print(f'{tag} interior cells={m.sum()} vs exact: L2 {e2:.3e} Linf {ei:.3e}')
    if 'gl' in a and not sp:
        k = a.index('gl'); nmu, nph = int(a[k+1]), int(a[k+2])
        zq, wz = np.polynomial.legendre.leggauss(nmu); ph = 2*np.pi*np.arange(nph)/nph
        d = np.array([[zz, np.sqrt(1-zz*zz)*np.cos(p), np.sqrt(1-zz*zz)*np.sin(p)] for zz in zq for p in ph])
        w = np.repeat(wz/2.0, nph)/nph; rr = np.maximum(rv, 1e-30)
        mur = (np.outer(X[m], d[:, 0]) + np.outer(Y[m], d[:, 1]) + np.outer(Z[m], d[:, 2]))/rr[:, None]
        jq = np.sum(w[None, :]*S*(1.0 - np.exp(-SIG*ipath(rv[:, None], mur))), axis=1)
        e2q, eiq = err(J[m], jq)
        print(f'  vs same quadrature (GL {nmu}x{nph}): L2 {e2q:.3e} Linf {eiq:.3e}; quadrature-only L2 {err(jq, je)[0]:.3e}')
rmax = r.max() if sp else min(abs(x).max(), abs(y).max(), abs(z).max())
print(f'{tag} exterior profile: r/R  ncell  <J/J_exact>  min  max')
for lo in np.arange(1.0, rmax/R - 0.05 + 1e-9, 0.05):
    mm = (r >= lo*R) & (r < (lo + 0.05)*R)
    if not mm.any(): continue
    q = J[mm]/jex(r[mm])
    print(f'  {lo:.2f}-{lo+0.05:.2f} {mm.sum():7d}  {q.mean():.4f}  {q.min():.4f}  {q.max():.4f}')
