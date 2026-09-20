#!/usr/bin/env python3
"""Arm tables for tests_r6 (multi-turnover): mass flux, shell drift, and the
Eddington/flux hand-over budget from the per-leg mltfaces dumps.

Usage: an6.py <armdir> [label]

rt_profile.bin slots (red_giant.cpp:494-506):
    0 rho  1 v1  2 rho v1  3 v1^2  4 v2^2+v3^2  5 wtemp = p/rho  6 eint  7 v1(eint+p)
mltfaces columns (red_giant.cpp:3517 ff):
    0 i  1 r  2 T_f  3 p_f  4 grad  5 grad_ad  6 x  7 F_mlt  8 F_rad  9 F_used
   10 F_req 11 F_raddiff 12 w_blend 13 F_rad_col0 14 F_conv_res 15 grad_rad 16 D 17 F_2s
"""
import sys
import os
import glob
import re
import numpy as np

TURN = 4.705e3
RSTAR = 2.3717e11
LSTAR = 2.3066e38
MSTAR = 6.2651e33
GRAV = 6.67430e-8
CLIGHT = 2.99792458e10
VMLT = 1.45e7
FRACS = [0.60, 0.70, 0.80, 0.90, 0.97]
GFRACS = [0.70, 0.80, 0.90]
OPAC = ('/viper/u2/jinma/ATHENAK/bench/hestar_fecz/'
        'rosseland_he_x0.0_z0.02.txt')


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
    o = np.argsort(pt)
    pt, Q = pt[o], Q[o]
    nx1 = len(x1v)
    idx = [int(np.argmin(np.abs(x1v - f * RSTAR))) for f in FRACS]
    ltg, ldg, tab = read_opac(OPAC)
    cf = sorted(glob.glob(os.path.join(d, "column_*.txt")))
    C = np.loadtxt(cf[0])
    fac = np.interp(x1v, C[:, 0], C[:, 2]) / Q[0, 5, :]
    edd = LSTAR / (4.0 * np.pi * CLIGHT * GRAV * MSTAR)

    xf = np.empty(nx1 + 1)
    xf[1:-1] = 0.5 * (x1v[:-1] + x1v[1:])
    xf[0] = x1v[0] - (xf[1] - x1v[0])
    xf[-1] = x1v[-1] + (x1v[-1] - xf[-2])
    vol = 4.0 / 3.0 * np.pi * (xf[1:]**3 - xf[:-1]**3)

    hst = sorted(glob.glob(os.path.join(d, "*.hydro.hst")))
    H = np.concatenate([np.loadtxt(h, comments="#", ndmin=2) for h in hst])
    H = H[np.argsort(H[:, 0])]

    tend = pt[-1]
    nq = int(np.floor(tend / (0.25 * TURN))) + 1
    tq = [k * 0.25 * TURN for k in range(nq)]
    if tend - tq[-1] > 1.0:
        tq.append(tend)

    print("\n===== ARM %s : %d profile records, t_end = %.5g s = %.3f turnover"
          % (lab, len(pt), tend, tend / TURN))
    print("hst total mass: %.6e -> %.6e g  (drift %+.3e relative)"
          % (H[0, 2], H[-1, 2], H[-1, 2] / H[0, 2] - 1.0))

    print("\n-- TABLE 1  Mdot(r) = 4 pi r^2 <rho v1> [g/s], + = outward ; "
          "and the domain mass")
    print("  t/turn " + "".join("   r/R=%.2f" % f for f in FRACS)
          + "  r/R=%.3f(top)   M(prof)/M0-1" % (x1v[-1] / RSTAR))
    for tt in tq:
        j = int(np.argmin(np.abs(pt - tt)))
        row = "  %6.2f " % (pt[j] / TURN)
        for i in idx + [nx1 - 1]:
            row += " %+9.3e" % (4.0 * np.pi * x1v[i]**2 * Q[j, 2, i])
        mtot = float((vol * Q[j, 0, :]).sum())
        m0 = float((vol * Q[0, 0, :]).sum())
        row += "      %+.3e" % (mtot / m0 - 1.0)
        print(row)

    print("\n-- TABLE 2  shell-integrated mass per zone [g], and the drift PER TURNOVER")
    zon = [("0.50-0.68 R (below FeCZ)", 0.50, 0.68),
           ("0.68-0.92 R (FeCZ, Gamma>1)", 0.68, 0.92),
           ("0.92-0.96 R (the swept shell)", 0.92, 0.96),
           ("0.96-1.01 R (top atmosphere)", 0.96, 1.015)]
    nturn = int(np.floor(tend / TURN))
    tw = [k * TURN for k in range(nturn + 1)]
    if tend - tw[-1] > 1.0:
        tw.append(tend)
    jw = [int(np.argmin(np.abs(pt - a))) for a in tw]
    print("  zone                          "
          + "".join(" t=%5.2f  " % (pt[j] / TURN) for j in jw))
    for nm, a, b in zon:
        s = (x1v >= a * RSTAR) & (x1v < b * RSTAR)
        m = [float((vol[s] * Q[j, 0, s]).sum()) for j in jw]
        print("  %-29s" % nm + "".join(" %9.3e" % v for v in m))
        print("  %-29s" % "   d/turn [%]"
              + "".join("     --   " if k == 0 else
                        " %+7.2f  " % (100.0 * (m[k] / m[k - 1] - 1.0))
                        for k in range(len(m))))

    print("\n-- TABLE 3  the three diagnostic radii: rho drift, v1, Gamma = kappa L/"
          "(4 pi c G M)")
    for nm, f in (("base", 0.635), ("kappa peak", 0.76), ("top", 0.97)):
        i = int(np.argmin(np.abs(x1v - f * RSTAR)))
        print("  %-10s r/R=%.3f i=%d" % (nm, x1v[i] / RSTAR, i))
        print("    t/turn   rho/rho0-1      v1      vr_rms/vMLT   T[K]    Gamma")
        for tt in tq:
            j = int(np.argmin(np.abs(pt - tt)))
            g = kappa_of(ltg, ldg, tab, Q[j, 0, i], fac[i] * Q[j, 5, i]) * edd
            print("   %7.2f   %+9.3f  %+9.2e   %8.3f   %8.3e  %6.3f"
                  % (pt[j] / TURN, Q[j, 0, i] / Q[0, 0, i] - 1.0, Q[j, 1, i],
                     np.sqrt(max(Q[j, 3, i], 0.0)) / VMLT,
                     fac[i] * Q[j, 5, i], g))

    print("\n-- TABLE 4  THE HAND-OVER BUDGET from the per-leg mltfaces dumps.")
    print("   Gamma_rad = Gamma x F_2s/F_req is the Eddington ratio the RADIATION alone")
    print("   carries; F_mlt/F_req is the sub-grid closure, F_res/F_req the RESOLVED")
    print("   convective flux (mltfaces column F_conv_res).")
    mfs = sorted(glob.glob(os.path.join(d, "mltfaces_*.txt")))
    print("   t/turn " + "".join("   r/R=%.2f: G_rad F_2s/Fq F_mlt/Fq F_res/Fq"
                                 % f for f in GFRACS))
    for mf in mfs:
        M = np.loadtxt(mf, comments="#", ndmin=2)
        k = int(re.search(r"mltfaces_(\d+)", mf).group(1))
        tt = k * 0.5 * TURN
        j = int(np.argmin(np.abs(pt - tt)))
        if abs(pt[j] - tt) > 0.5 * TURN:
            continue
        row = "   %6.2f " % (pt[j] / TURN)
        for f in GFRACS:
            im = int(np.argmin(np.abs(M[:, 1] - f * RSTAR)))
            i = int(np.argmin(np.abs(x1v - M[im, 1])))
            g = kappa_of(ltg, ldg, tab, Q[j, 0, i], fac[i] * Q[j, 5, i]) * edd
            fq = M[im, 10]
            row += "        %6.3f  %7.3f %7.3f  %+7.3f" % (
                g * M[im, 17] / fq, M[im, 17] / fq, M[im, 7] / fq, M[im, 14] / fq)
        print(row)


main()
