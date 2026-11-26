#!/usr/bin/env python3
"""SCPR: J of the homogeneous sphere (R 0.6, sigma 5, S 1) against exact transport.
usage: sph_err.py BINFILE VAR [gl NMU NPHI]   (VAR: sc_rad_J for PR #798, vet_vj for ours)
Interior cells r < 0.36 (the PR's r_mask 0.6 R). Prints rel L2 and Linf of J vs the exact
continuous J, and, with 'gl', also vs the same-quadrature J of our GL product set (pole on x1)."""
import os, sys, numpy as np
sys.path.insert(0, os.environ.get('BINPY', '.'))  # dir of vis/python/bin_convert.py
import bin_convert as bc
R, SIG, S = 0.6, 5.0, 1.0
f = bc.read_binary(sys.argv[1]); var = sys.argv[2]
assert f['n_mbs'] == 1, 'one MeshBlock expected'
J = np.asarray(f['mb_data'][var][0], dtype=float)
x = np.asarray(f['mb_x1v'][0]); y = np.asarray(f['mb_x2v'][0]); z = np.asarray(f['mb_x3v'][0])
Z, Y, X = np.meshgrid(z, y, x, indexing='ij'); r = np.sqrt(X*X + Y*Y + Z*Z)
m = r < 0.6*R
def chord(rr, mu):  # path to the surface looking back along -n from a point at radius rr
    return rr*mu + np.sqrt(R*R - rr*rr*(1.0 - mu*mu))
mu, wq = np.polynomial.legendre.leggauss(400)
rv = r[m]
jex = 0.5*S*np.sum(wq[None, :]*(1.0 - np.exp(-SIG*chord(rv[:, None], mu[None, :]))), axis=1)
def err(a, b): return np.sqrt(np.sum((a - b)**2)/np.sum(b*b)), np.max(np.abs(a - b)/b)
e2, ei = err(J[m], jex)
print(f'{sys.argv[1]} N={J.shape[2]} cells={m.sum()} vs exact: L2 {e2:.3e} Linf {ei:.3e}')
if 'gl' in sys.argv:
    k = sys.argv.index('gl'); nmu, nph = int(sys.argv[k+1]), int(sys.argv[k+2])
    zq, wz = np.polynomial.legendre.leggauss(nmu); ph = 2*np.pi*np.arange(nph)/nph
    d = np.array([[zz, np.sqrt(1-zz*zz)*np.cos(p), np.sqrt(1-zz*zz)*np.sin(p)] for zz in zq for p in ph])
    w = np.repeat(wz/2.0, nph)/nph
    px, py, pz = X[m], Y[m], Z[m]; rr = np.maximum(rv, 1e-30)
    mur = (np.outer(px, d[:, 0]) + np.outer(py, d[:, 1]) + np.outer(pz, d[:, 2]))/rr[:, None]
    jq = np.sum(w[None, :]*S*(1.0 - np.exp(-SIG*chord(rv[:, None], mur))), axis=1)
    e2q, eiq = err(J[m], jq)
    print(f'  vs same quadrature (GL {nmu}x{nph}): L2 {e2q:.3e} Linf {eiq:.3e};  quadrature-only error L2 {err(jq, jex)[0]:.3e}')
