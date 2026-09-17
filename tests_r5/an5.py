#!/usr/bin/env python3
"""Mass budget, Eddington ratio and shell profiles for a he4_presn arm.

Usage: an5.py <armdir> [label]

rt_profile.bin slots (red_giant.cpp:494-506):
    0 rho  1 v1  2 rho v1  3 v1^2  4 v2^2+v3^2  5 wtemp = p/rho  6 eint  7 v1(eint+p)
"""
import sys
import os
import glob
import numpy as np

TURN = 4.705e3
RSTAR = 2.3717e11
LSTAR = 2.3066e38
MSTAR = 6.2651e33
GRAV = 6.67430e-8
CLIGHT = 2.99792458e10
VMLT = 1.45e7
FRACS = [0.60, 0.70, 0.80, 0.90, 0.97]
OPAC = ('/viper/u2/jinma/ATHENAK/bench/hestar_fecz/'
        'rosseland_he_x0.0_z0.02.txt')
COL = None          # the run's own initial column dump, for the T[K] calibration


def read_prof(fn):
    ts, qs, x1v = [], [], None
    with open(fn, "rb") as f:
        while True:
            b = f.read(8)
            if len(b) < 8:
                break
            t = np.frombuffer(b, "<f8")[0]
            nx1, nvar = np.frombuffer(f.read(8), "<i4")
            x = np.frombuffer(f.read(8 * nx1), "<f8")
            q = np.frombuffer(f.read(8 * nvar * nx1), "<f8").reshape(nvar, nx1)
            if x1v is None:
                x1v = x.copy()
            ts.append(t)
            qs.append(q.copy())
    return np.array(ts), x1v, np.array(qs)


def read_opac(fn):
    """regular (log10 T, log10 rho) grid, T slowest; header '# nT nD lT0 dlT lD0 dlD'."""
    nt = None
    vals = []
    for ln in open(fn):
        if ln.lstrip().startswith('#'):
            p = ln.lstrip('# ').split()
            if len(p) == 6 and nt is None:
                try:
                    nt, nd = int(p[0]), int(p[1])
                    lt0, dlt, ld0, dld = [float(x) for x in p[2:]]
                except ValueError:
                    nt = None
            continue
        if ln.strip():
            vals.append(float(ln))
    tab = np.array(vals[:nt * nd]).reshape(nt, nd)
    return (lt0 + dlt * np.arange(nt)), (ld0 + dld * np.arange(nd)), tab


def kappa_of(ltg, ldg, tab, rho, T):
    lt = np.clip(np.log10(T), ltg[0], ltg[-1])
    ld = np.clip(np.log10(rho), ldg[0], ldg[-1])
    it = min(max(int((lt - ltg[0]) / (ltg[1] - ltg[0])), 0), len(ltg) - 2)
    jd = min(max(int((ld - ldg[0]) / (ldg[1] - ldg[0])), 0), len(ldg) - 2)
    ft = (lt - ltg[it]) / (ltg[1] - ltg[0])
    fd = (ld - ldg[jd]) / (ldg[1] - ldg[0])
    return 10.0**((1 - ft) * (1 - fd) * tab[it, jd] + (1 - ft) * fd * tab[it, jd + 1]
                  + ft * (1 - fd) * tab[it + 1, jd] + ft * fd * tab[it + 1, jd + 1])


def main():
    d = sys.argv[1]
    lab = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(d.rstrip('/'))
    pt, x1v, Q = read_prof(os.path.join(d, "rt_profile.bin"))
    nx1 = len(x1v)
    idx = [int(np.argmin(np.abs(x1v - f * RSTAR))) for f in FRACS]
    ltg, ldg, tab = read_opac(OPAC)
    # T[K] = wtemp * mu m_u / k.  mu is not dumped, so calibrate the factor per radius
    # against the run's own initial column (T[K] there is exact) and hold it fixed; the
    # shell means then give T[K] to the accuracy of d(mu)/mu, a few per cent.
    cf = glob.glob(os.path.join(d, "column_*.txt"))
    C = np.loadtxt(cf[0])
    fac = np.interp(x1v, C[:, 0], C[:, 2]) / Q[0, 5, :]
    edd = LSTAR / (4.0 * np.pi * CLIGHT * GRAV * MSTAR)   # Gamma = kappa * edd

    # cell volumes from the face positions (x1f reconstructed as cell-centre midpoints)
    xf = np.empty(nx1 + 1)
    xf[1:-1] = 0.5 * (x1v[:-1] + x1v[1:])
    xf[0] = x1v[0] - (xf[1] - x1v[0])
    xf[-1] = x1v[-1] + (x1v[-1] - xf[-2])
    vol = 4.0 / 3.0 * np.pi * (xf[1:]**3 - xf[:-1]**3)

    hst = sorted(glob.glob(os.path.join(d, "*.hydro.hst")))
    H = np.concatenate([np.loadtxt(h, comments="#", ndmin=2) for h in hst])
    H = H[np.argsort(H[:, 0])]

    print("\n===== ARM %s : %d profile records, t_end = %.4g s = %.3f turnover"
          % (lab, len(pt), pt[-1], pt[-1] / TURN))
    print("hst total mass: %.6e -> %.6e g  (drift %+.3e relative, %+.3e g)"
          % (H[0, 2], H[-1, 2], H[-1, 2] / H[0, 2] - 1.0, H[-1, 2] - H[0, 2]))

    print("\n-- TABLE 1  shell mass flux Mdot(r) = 4 pi r^2 <rho v1> [g/s], + = outward")
    print("  t/turn " + "".join("   r/R=%.2f" % f for f in FRACS)
          + "  r/R=%.3f(top)   M(domain)/M0-1" % (x1v[-1] / RSTAR))
    for k in range(11):
        j = int(np.argmin(np.abs(pt - k * 0.1 * TURN)))
        row = "  %6.2f " % (pt[j] / TURN)
        for i in idx + [nx1 - 1]:
            row += " %+9.3e" % (4.0 * np.pi * x1v[i]**2 * Q[j, 2, i])
        mtot = float((vol * Q[j, 0, :]).sum())
        m0 = float((vol * Q[0, 0, :]).sum())
        row += "      %+.3e" % (mtot / m0 - 1.0)
        print(row)

    print("\n-- TABLE 2  shell-integrated mass in three zones [g] and its drift")
    zon = [("0.50-0.68 R (below FeCZ)", 0.50, 0.68),
           ("0.68-0.92 R (FeCZ, Gamma>1)", 0.68, 0.92),
           ("0.92-0.96 R (the swept shell)", 0.92, 0.96),
           ("0.96-1.01 R (top atmosphere)", 0.96, 1.015)]
    print("  zone                            " + "".join("   t=%.2f" % a
          for a in (0.0, 0.25, 0.5, 0.75, 1.0)) + "     d(1 turn)")
    js = [int(np.argmin(np.abs(pt - a * TURN))) for a in (0.0, .25, .5, .75, 1.0)]
    for nm, a, b in zon:
        s = (x1v >= a * RSTAR) & (x1v < b * RSTAR)
        m = [float((vol[s] * Q[j, 0, s]).sum()) for j in js]
        print("  %-31s" % nm + "".join(" %9.3e" % v for v in m)
              + "   %+7.1f %%" % (100.0 * (m[-1] / m[0] - 1.0)))

    print("\n-- TABLE 3  shell means and the Eddington ratio Gamma = kappa F/(c g)")
    print("   r/R     r[cm]    rho(0)    rho(.5)   rho(1)    v1(.5)     v1(1)"
          "      T0[K]     G(0)   G(.5)   G(1)")
    for i in range(0, nx1, 3):
        g = []
        for j in (js[0], js[2], js[4]):
            g.append(kappa_of(ltg, ldg, tab, Q[j, 0, i], fac[i] * Q[j, 5, i]) * edd)
        print("  %.3f %9.4e %9.3e %9.3e %9.3e %+9.2e %+9.2e %9.3e  %6.3f %6.3f %6.3f"
              % (x1v[i] / RSTAR, x1v[i], Q[js[0], 0, i], Q[js[2], 0, i],
                 Q[js[4], 0, i], Q[js[2], 1, i], Q[js[4], 1, i],
                 fac[i] * Q[js[0], 5, i], g[0], g[1], g[2]))

    print("\n-- TABLE 4  v1 against c_s, v_esc and v_MLT")
    print("   r/R    v1(.5)     v1(1)      c_s(1)     v_esc      v1/c_s  v1/v_esc"
          "  v1/v_MLT  vr_rms/v_MLT")
    for f, i in zip(FRACS, idx):
        ves = np.sqrt(2.0 * GRAV * MSTAR / x1v[i])
        cs = np.sqrt(5.0 / 3.0 * (2.0 / 3.0) * Q[js[4], 6, i] / Q[js[4], 0, i])
        v = Q[js[4], 1, i]
        print("  %.2f  %+9.3e %+9.3e %9.3e %9.3e  %8.4f %8.4f %8.4f %8.4f"
              % (f, Q[js[2], 1, i], v, cs, ves, v / cs, v / ves, v / VMLT,
                 np.sqrt(max(Q[js[4], 3, i], 0.0)) / VMLT))


if __name__ == "__main__":
    main()
