#!/usr/bin/env python3
"""rph.py BIN_CONVERT_DIR KAPPA_TABLE LTOP hydro_w.bin [...]   (RSGCONT 2026-11-24)

Photospheric radius of a cubed-sphere he_star dump (hydro_w + m1 with the same number; the m1 file is the
hydro_w name with hydro_w -> m1).  Old-scheme and clean-scheme dumps are read the same way.

  BIN_CONVERT_DIR  directory with AthenaK's vis/python/bin_convert.py
  KAPPA_TABLE      the run's Rosseland table (<problem>/he_opac_table; RSG: rosseland_mesakap_4c0e8b64d26a.txt)
  LTOP             L_top [erg/s] for T_eff; or the run's user.hst (L_top = column 6 at the dump time, nearest row);
                   0 = 4 pi r^2 <F1> at the top cell

Columns are weighted by their solid angle on the gnomonic equiangular cubed sphere
(x = tan xi, y = tan eta, xi = x2 pi/4: dOmega = (1+x^2)(1+y^2)/(1+x^2+y^2)^1.5 dxi deta).
kappa_R(T_rad, rho), T_rad = (E/a)^(1/4), bilinear in (log T, log rho) as rsgjump/tools/red.py.
Printed per dump:
  R23_med / R23_mean / p10 / p90   per-column tau_R = 2/3 radius (tau summed from the top face inward)
  R1_mean                          per-column tau_R = 1 radius, solid-angle mean
  Rsh23 / Rsh1                     tau = 2/3 / 1 of the SHELL-averaged opacity chi(r) = <kappa rho>_Omega
  R_Teff                           outermost r with a <T_rad^4>_Omega = L / (4 pi sigma r^2)  (pi c r^2 <E> = L)
  T_eff(R23_mean)                  (L / (4 pi sigma R23_mean^2))^(1/4)
"""
import sys
import numpy as np

sys.path.insert(0, sys.argv[1])
import bin_convert as bc  # noqa: E402

A_RAD, C_L = 7.5657e-15, 2.99792458e10
SIG = A_RAD*C_L/4.0
import os  # noqa: E402
tab = sys.argv[2]
HST = np.loadtxt(sys.argv[3]) if os.path.isfile(sys.argv[3]) else None
ltop_arg = 0.0 if HST is not None else float(sys.argv[3])

L = open(tab).readlines()
g = [ln for ln in L if ln.startswith('#')][4].split()[1:]
nT, nD = int(g[0]), int(g[1])
lT0, dlT, lD0, dlD = map(float, g[2:6])
K = np.array([float(ln) for ln in L if not ln.startswith('#')]).reshape(nT, nD)


def lkap(lT, lD):
    x = np.clip((lT - lT0)/dlT, 0, nT - 1.001)
    y = np.clip((lD - lD0)/dlD, 0, nD - 1.001)
    i, j = x.astype(int), y.astype(int)
    fx, fy = x - i, y - j
    return (K[i, j]*(1-fx)*(1-fy) + K[i+1, j]*fx*(1-fy) + K[i, j+1]*(1-fx)*fy
            + K[i+1, j+1]*fx*fy)


def wq(v, w, q):
    """solid-angle weighted quantile"""
    o = np.argsort(v)
    c = np.cumsum(w[o])
    return v[o][np.searchsorted(c, q*c[-1])]


def cross(x1f, tau, chi, lev):
    """per row: radius where tau (from the top, at each cell's inner face) reaches lev, linear inside the cell"""
    n = tau.shape[-1]
    ip = n - 1 - np.argmax(tau[..., ::-1] >= lev, axis=-1)   # outermost cell whose inner-face tau >= lev
    rows = np.arange(tau.shape[0])
    tout = np.where(ip + 1 < n, tau[rows, np.minimum(ip + 1, n - 1)], 0.0)   # tau at that cell's outer face
    return x1f[ip + 1] - (lev - tout)/chi[rows, ip]


for fh in sys.argv[4:]:
    h = bc.read_binary(fh)
    m = bc.read_binary(fh.replace('hydro_w', 'm1'))
    x1f = np.array(h['mb_x1f'][0])
    x1v = np.array(h['mb_x1v'][0])
    n1 = len(x1v)
    dr = np.diff(x1f)
    nb = h['n_mbs']

    def cat(M, v):
        return np.concatenate([np.array(M['mb_data'][v][b]).reshape(-1, n1) for b in range(nb)])

    w = []
    for b in range(nb):
        x2f, x3f = np.array(h['mb_x2f'][b]), np.array(h['mb_x3f'][b])
        x2v, x3v = np.array(h['mb_x2v'][b]), np.array(h['mb_x3v'][b])
        X, Y = np.tan(np.pi/4*x2v)[None, :], np.tan(np.pi/4*x3v)[:, None]
        d = (np.pi/4)**2*np.diff(x2f)[None, :]*np.diff(x3f)[:, None]
        w.append(((1 + X*X)*(1 + Y*Y)/(1 + X*X + Y*Y)**1.5*d).ravel())   # (k, j) order = reshape(-1, n1)
    w = np.concatenate(w)
    rho = cat(h, 'dens')
    E = cat(m, 'm1_e')
    F1 = cat(m, 'm1_f1')
    Tr = (np.maximum(E, 1e-30)/A_RAD)**0.25
    chi = 10**lkap(np.log10(Tr), np.log10(np.maximum(rho, 1e-30)))*rho
    tau = np.cumsum((chi*dr)[:, ::-1], axis=1)[:, ::-1]
    r23 = cross(x1f, tau, chi, 2/3.)
    r1 = cross(x1f, tau, chi, 1.0)
    W = w.sum()
    chis = (w[:, None]*chi).sum(0)/W
    taus = np.cumsum((chis*dr)[::-1])[::-1]
    rs23 = cross(x1f, taus[None, :], chis[None, :], 2/3.)[0]
    rs1 = cross(x1f, taus[None, :], chis[None, :], 1.0)[0]
    Es = (w[:, None]*E).sum(0)/W
    Fs = (w[:, None]*F1).sum(0)/W
    if HST is not None:
        ltop_arg = HST[np.argmin(abs(HST[:, 0] - h['time'])), 5]
    lt = ltop_arg if ltop_arg > 0 else 4*np.pi*x1v[-1]**2*Fs[-1]
    q = np.pi*C_L*x1v**2*Es - lt          # > 0 inside R_Teff
    i = n1 - 1 - np.argmax(q[::-1] >= 0)
    rte = x1v[i] + (x1v[i+1] - x1v[i])*q[i]/(q[i] - q[i+1]) if i + 1 < n1 else x1v[-1]
    r23m = (w*r23).sum()/W
    print('%s t %.6e cyc %d L %.5e | R23 med %.5e mean %.5e p10 %.5e p90 %.5e | R1 mean %.5e | '
          'Rsh23 %.5e Rsh1 %.5e | R_Teff %.5e | T_eff(R23mean) %.1f K | sum(w)/4pi %.6f'
          % (fh.split('/')[-1], h['time'], h['cycle'], lt, wq(r23, w, 0.5), r23m, wq(r23, w, 0.1),
             wq(r23, w, 0.9), (w*r1).sum()/W, rs23, rs1, rte, (lt/(4*np.pi*SIG*r23m**2))**0.25,
             W/(4*np.pi)), flush=True)
