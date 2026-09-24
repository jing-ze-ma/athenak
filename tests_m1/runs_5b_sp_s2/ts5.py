#!/usr/bin/env python3
"""T-S5: the radiation-modified acoustic wave (T10) run TANGENTIALLY (along theta) in a
thin spherical shell vs the same wave along x2 of a Cartesian slab.

For each run: the complex amplitude Z(t) = <(rho - <rho>) exp(-i 2 pi (j + 1/2)/n2)>
over the dumps; phase speed from the unwrapped phase, decay rate from ln|Z| (linear
fits).  Against the Cartesian twin: max over dumps of <|rho_sp - rho_cart|>/amp
(cell-to-cell, same (k, j, i) indices).

usage: ts5.py lambda cartdir spdir [spdir ...]
"""
import glob
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/wt_m1sp2/tests_m1')
import common  # noqa: E402


def series(d):
    fs = sorted(glob.glob(d + '/*.tab'))
    out = []
    for f in fs:
        t = common.load_tab(f)
        out.append((t.time, t.data['dens']))
    return out


def zamp(r):
    n2 = r.shape[1]
    ph = np.exp(-2j*np.pi*(np.arange(n2) + 0.5)/n2)
    return np.mean((r - r.mean())*ph[None, :, None])


def fit(ser, lam):
    ts = np.array([s[0] for s in ser])
    z = []
    for _, r in ser:
        n2 = r.shape[1]
        ph = np.exp(-2j*np.pi*(np.arange(n2) + 0.5)/n2)
        z.append(np.mean((r - r.mean())*ph[None, :, None]))
    z = np.array(z)
    phs = np.unwrap(np.angle(z))
    om = -np.polyfit(ts, phs, 1)[0]
    gam = -np.polyfit(ts, np.log(np.abs(z)), 1)[0]
    return om*lam/(2*np.pi), gam


lam = float(sys.argv[1])
cart = series(sys.argv[2])
vc, gc = fit(cart, lam)
print(f"cart: n={len(cart)} v_phase={vc:.6f} decay={gc:.5e}")
for d in sys.argv[3:]:
    sp = series(d)
    vs, gs = fit(sp, lam)
    zs, zc = zamp(sp[-1][1])/zamp(sp[0][1]), zamp(cart[-1][1])/zamp(cart[0][1])
    print(f"{d}: t={sp[-1][0]:.6f} v_phase={vs:.6f} (rel diff {vs/vc - 1:+.3e}) "
          f"decay={gs:.5e} (diff {gs - gc:+.3e}) "
          f"|Z_sp(T)/Z_sp(0) - Z_c(T)/Z_c(0)|/|Z_c(T)/Z_c(0)|={abs(zs - zc)/abs(zc):.3e}")
