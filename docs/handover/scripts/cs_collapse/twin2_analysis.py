#!/usr/bin/env python3
"""ctl2_s01_twin2: the RT flux components (rt_Ft, rt_Fb, rt_Qs, rt_Em) in the collapsing
cubed-sphere vertex column.

Three arms, all restarted from the same rot-5.700 seed:
  exp_g3  rt_semi_implicit = false, diag_gid = 3  -- the block the exp collapse is in
  exp_g6  rt_semi_implicit = false, diag_gid = 6  -- the twin's original block
  si_g6   rt_semi_implicit = true,  diag_gid = 6  -- the runaway column

Everything here is descriptive.  No modelling, no fitting.
"""

import glob
import os
import sys

import numpy as np

sys.path.insert(0, "/viper/u2/jinma/ATHENAK/athenak/docs/handover/scripts")
from read_cyclediag import load           # noqa: E402

ROOT = "/viper/u2/jinma/ATHENAK/bench/cs_ens/ctl2_s01_twin2"
OUT = os.path.join(ROOT, "analysis")
SIGMA = 5.6704e-5

# arm -> (k, j) of the column of interest
ARMS = {"exp_g3": (16, 17), "exp_g6": (17, 2), "si_g6": (17, 2)}

# --- the radial grid (mesh/use_grid_stretch_r_poly, ctl2 coefficients) ----------------
R0, R1 = 9.44e9, 2.0556e10
CPOLY = (-0.068392, -2.191487, 2.464818, -1.366698)
NX1, NGHOST = 128, 2


def stretch_r_poly(r):
    xi = (np.asarray(r) - R0)/(R1 - R0)
    u = xi.copy()
    xik = xi.copy()
    for c in CPOLY:
        u += c*xik*(1.0 - xi)
        xik = xik*xi
    return R0 + (R1 - R0)*u


_du = (R1 - R0)/NX1
FACE = stretch_r_poly(R0 + np.arange(-NGHOST, NX1 + NGHOST + 1)*_du)
DX1 = np.diff(FACE)                       # DX1[i] = width of array cell i

CELLV = ("rt_T", "w_dens", "w_velx", "w_eint", "u_ener", "rt_src", "rt_de", "rt_clip",
         "rt_Ft", "rt_Fb", "rt_Qs", "rt_Em", "cond_f1", "cond_kappa", "cond_keff",
         "cond_dtcell")
FACEV = ("rad_w", "rad_tauf")


def read_arm(arm):
    """Whole column (k, j), all i, every cycle.  Cached as an npz."""
    cache = os.path.join(OUT, "%s_col.npz" % arm)
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return {k: z[k] for k in z.files}
    kk, jj = ARMS[arm]
    files = sorted(glob.glob(os.path.join(ROOT, arm, "cyclediag", "*.dat")))
    rec = {v: [] for v in CELLV + FACEV}
    for k in ("cycle", "time", "dt", "icut", "nclip_blk", "maxT_blk"):
        rec[k] = []
    skipped = []
    for fn in files:
        try:
            d = load(fn)
        except (ValueError, IndexError) as e:
            skipped.append("%s: %s" % (os.path.basename(fn), str(e)))
            continue
        rec["cycle"].append(d["cycle"])
        rec["time"].append(d["time"])
        rec["dt"].append(d["dt"])
        rec["icut"].append(d["rt_icut"][kk, jj])
        sl = (slice(d["ks"], d["ke"]+1), slice(d["js"], d["je"]+1),
              slice(d["is"], d["ie"]+1))
        rec["nclip_blk"].append(int(np.count_nonzero(d["rt_clip"][sl] > 0.5)))
        rec["maxT_blk"].append(float(d["rt_T"][sl].max()))
        for v in CELLV + FACEV:
            rec[v].append(d[v][kk, jj].copy())
        rec["_hdr"] = (d["is"], d["ie"], d["n1"])
    hdr = rec.pop("_hdr", (2, 129, 132))
    out = {k: np.array(v) for k, v in rec.items()}
    out["is"] = np.array(hdr[0])
    out["ie"] = np.array(hdr[1])
    out["skipped"] = np.array(skipped)
    np.savez_compressed(cache, **out)
    return out


def root_cell(rec, lag=20, fac=1.5):
    """First cycle at which any active cell's T exceeds fac x its value lag cycles
    earlier.  Returns (index, cycle, i) or (None, None, None)."""
    T = rec["rt_T"]
    i0, i1 = int(rec["is"]), int(rec["ie"])
    for n in range(lag, T.shape[0]):
        r = T[n, i0:i1+1]/np.maximum(T[n-lag, i0:i1+1], 1.0e-30)
        if np.any(r > fac):
            return n, int(rec["cycle"][n]), i0 + int(np.argmax(r))
    return None, None, None


def trigger_table(rec, arm, lag=20, fac=1.5):
    """Every active i that ever trips the 1.5x-in-20-cycles test, by first trip."""
    T = rec["rt_T"]
    i0, i1 = int(rec["is"]), int(rec["ie"])
    first = {}
    for n in range(lag, T.shape[0]):
        r = T[n, i0:i1+1]/np.maximum(T[n-lag, i0:i1+1], 1.0e-30)
        for w in np.nonzero(r > fac)[0]:
            i = i0 + int(w)
            if i not in first:
                first[i] = (n, int(rec["cycle"][n]), float(r[w]))
    print("\n%s: every active i that ever trips T > %.1f x T(%d cycles earlier), "
          "ordered by first trip" % (arm, fac, lag))
    print("%6s %8s %8s %10s %13s %13s" % ("i", "index", "cycle", "ratio",
                                          "T(n-lag)", "T(n)"))
    for i, (n, c, r) in sorted(first.items(), key=lambda kv: kv[1][0]):
        print("%6d %8d %8d %10.4f %13.5e %13.5e"
              % (i, n, c, r, T[n-lag, i], T[n, i]))
    if not first:
        print("  (none)")


def maxT_table(rec, arm, every):
    print("\n" + "=" * 72)
    print("%s: per-cycle max T in the column (k,j) = (%d,%d), active i only"
          % (arm, ARMS[arm][0], ARMS[arm][1]))
    print("=" * 72)
    print("%8s %13s %12s %6s %12s" % ("cycle", "time", "maxT[K]", "i", "dt"))
    i0, i1 = int(rec["is"]), int(rec["ie"])
    T = rec["rt_T"][:, i0:i1+1]
    am = T.argmax(axis=1) + i0
    n = T.shape[0]
    for m in range(n):
        if m % every and m < n - 10:
            continue
        print("%8d %13.7g %12.5e %6d %12.5e"
              % (rec["cycle"][m], rec["time"][m], T[m].max(), am[m], rec["dt"][m]))


def cond_div(rec, i):
    """index-space conduction divergence D = -(cond_f1[i+1] - cond_f1[i])."""
    return -(rec["cond_f1"][:, i+1] - rec["cond_f1"][:, i])


def root_table(rec, arm, ir, n_root, nback=40):
    print("\n" + "=" * 235)
    print("%s ROOT CELL i = %d, column (k,j) = (%d,%d), from %d cycles before the root "
          "cycle to the end of the run" % (arm, ir, ARMS[arm][0], ARMS[arm][1], nback))
    print("delta_u = u_ener(n+1) - u_ener(n).  D is the INDEX-SPACE conduction "
          "divergence -(cond_f1[i+1]-cond_f1[i]); D/dx1 is per volume.")
    print("dx1 = %.6e cm.  sT4 = sigma*T^4 with sigma = %g." % (DX1[ir], SIGMA))
    print("=" * 235)
    hdr = ("%8s %10s %11s %11s %12s %11s %12s %12s %12s %12s %12s %11s %11s %11s "
           "%9s %9s %10s %10s %12s %10s %6s")
    print(hdr % ("cycle", "dt", "T[K]", "rho", "v_r", "e_int", "delta_u", "rt_src",
                 "rt_de", "rt_Ft", "rt_Fb", "rt_Qs", "rt_Em", "sT4",
                 "radw_i", "radw_i1", "tauf_i", "tauf_i1", "D/dx1", "dtcell",
                 "clip"))
    n = len(rec["cycle"])
    lo = max(0, n_root - nback)
    u = rec["u_ener"][:, ir]
    D = cond_div(rec, ir)
    for m in range(lo, n):
        du = (u[m+1] - u[m]) if m + 1 < n else np.nan
        T = rec["rt_T"][m, ir]
        print(hdr % ("%d" % rec["cycle"][m], "%.3e" % rec["dt"][m],
                     "%.4e" % T, "%.4e" % rec["w_dens"][m, ir],
                     "%.5e" % rec["w_velx"][m, ir], "%.4e" % rec["w_eint"][m, ir],
                     "%.5e" % du, "%.5e" % rec["rt_src"][m, ir],
                     "%.5e" % rec["rt_de"][m, ir], "%.5e" % rec["rt_Ft"][m, ir],
                     "%.5e" % rec["rt_Fb"][m, ir], "%.4e" % rec["rt_Qs"][m, ir],
                     "%.4e" % rec["rt_Em"][m, ir], "%.4e" % (SIGMA*T**4),
                     "%.3e" % rec["rad_w"][m, ir], "%.3e" % rec["rad_w"][m, ir+1],
                     "%.4e" % rec["rad_tauf"][m, ir],
                     "%.4e" % rec["rad_tauf"][m, ir+1],
                     "%.5e" % (D[m]/DX1[ir]), "%.4e" % rec["cond_dtcell"][m, ir],
                     "%.0f" % rec["rt_clip"][m, ir]))


def flux_sanity(rec, arm, n_root, back=5, only=None):
    i0, i1 = int(rec["is"]), int(rec["ie"])
    ms = only if only is not None else (n_root - back, n_root)
    for m in ms:
        if m < 0 or m >= len(rec["cycle"]):
            continue
        T = rec["rt_T"][m]
        Tmax = T[i0:i1+1].max()
        cap = 2.0*SIGMA*Tmax**4
        Qsmax = rec["rt_Qs"][m, i0:i1+1].max()
        print("\n" + "=" * 132)
        print("%s FLUX SANITY, cycle %d (n_root %s %d), column (k,j) = (%d,%d)"
              % (arm, rec["cycle"][m], "-" if m < n_root else "=", abs(n_root - m),
                 ARMS[arm][0], ARMS[arm][1]))
        print("column max T = %.5e K -> sigma*Tmax^4 = %.5e, flag threshold 2x = %.5e"
              % (Tmax, SIGMA*Tmax**4, cap))
        print("column max rt_Qs = %.5e   icut = %d   dt = %.5e"
              % (Qsmax, rec["icut"][m], rec["dt"][m]))
        print("rt_Fb[i] is the net flux on the BOTTOM face of cell i; the 'hotter "
              "adjacent cell' is max(T[i-1], T[i]).")
        print("=" * 132)
        print("%5s %12s %12s %14s %14s %12s %12s %12s %8s"
              % ("i", "T[K]", "rho", "rt_Fb", "rt_Em", "rad_w", "rad_tauf",
                 "|Fb|/sT4adj", "FLAG"))
        for i in range(i0, i1+1):
            Fb = rec["rt_Fb"][m, i]
            Tadj = max(T[i-1], T[i]) if i > 0 else T[i]
            ref = SIGMA*Tadj**4
            print("%5d %12.5e %12.5e %14.6e %14.6e %12.5e %12.5e %12.4e %8s"
                  % (i, T[i], rec["w_dens"][m, i], Fb, rec["rt_Em"][m, i],
                     rec["rad_w"][m, i], rec["rad_tauf"][m, i],
                     abs(Fb)/ref if ref > 0 else np.nan,
                     "***" if abs(Fb) > cap else ""))
        nflag = int(np.count_nonzero(np.abs(rec["rt_Fb"][m, i0:i1+2]) > cap))
        print("faces with |net flux| > 2 sigma Tmax^4: %d" % nflag)


def sign_locality(rec, arm, ir, n_root, nafter=6):
    print("\n" + "=" * 108)
    print("%s SIGN / LOCALITY: the %d cycles after the root cycle, i = %d .. %d"
          % (arm, nafter, ir-3, ir+3))
    print("-(Ft-Fb)/dx1 is recomputed here from the dumped Ft, Fb and the poly-stretch "
          "dx1; rt_src is the code's own value.")
    print("=" * 108)
    n = len(rec["cycle"])
    for m in range(n_root + 1, min(n_root + 1 + nafter, n)):
        print("\ncycle %d   dt = %.5e   time = %.8g"
              % (rec["cycle"][m], rec["dt"][m], rec["time"][m]))
        print("%5s %12s %14s %14s %14s %14s"
              % ("i", "T[K]", "rt_src", "-(Ft-Fb)/dx", "rt_Qs", "rt_Em"))
        for i in range(max(0, ir-3), min(rec["rt_T"].shape[1], ir+4)):
            div = -(rec["rt_Ft"][m, i] - rec["rt_Fb"][m, i])/DX1[i]
            print("%5d %12.5e %14.6e %14.6e %14.6e %14.6e"
                  % (i, rec["rt_T"][m, i], rec["rt_src"][m, i], div,
                     rec["rt_Qs"][m, i], rec["rt_Em"][m, i]))


def plots(recs, roots):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paths = []

    # 1. rt_Fb and T profiles at 6 cycles around the root cycle, both arms
    for arm in ("exp_g3", "si_g6"):
        rec = recs[arm]
        n_root, _, ir = roots[arm]
        if n_root is None:
            continue
        i0, i1 = int(rec["is"]), int(rec["ie"])
        ii = np.arange(i0, i1+1)
        ms = [m for m in range(n_root-3, n_root+3) if 0 <= m < len(rec["cycle"])]
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
        for m in ms:
            lab = "cycle %d" % rec["cycle"][m]
            ax[0].plot(ii, rec["rt_Fb"][m, i0:i1+1], label=lab)
            ax[1].semilogy(ii, rec["rt_T"][m, i0:i1+1], label=lab)
        for a, t in zip(ax, ("rt_Fb (bottom-face net flux)", "T [K]")):
            a.axvline(ir, color="k", ls=":", lw=1)
            a.set_xlabel("i (array index)")
            a.set_title("%s  %s" % (arm, t))
            a.legend(fontsize=7)
        ax[0].set_yscale("symlog", linthresh=1e3)
        fig.tight_layout()
        p = os.path.join(OUT, "profile_Fb_T_%s.png" % arm)
        fig.savefig(p, dpi=130)
        plt.close(fig)
        paths.append(p)

    # 1b. the same profiles at the END of the run
    for arm in ("exp_g3", "si_g6"):
        rec = recs[arm]
        i0, i1 = int(rec["is"]), int(rec["ie"])
        ii = np.arange(i0, i1+1)
        nlast = len(rec["cycle"]) - 1
        ms = [m for m in range(nlast-5, nlast+1) if m >= 0]
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
        for m in ms:
            lab = "cycle %d" % rec["cycle"][m]
            ax[0].plot(ii, rec["rt_Fb"][m, i0:i1+1], label=lab)
            ax[1].semilogy(ii, rec["rt_T"][m, i0:i1+1], label=lab)
        for a, t in zip(ax, ("rt_Fb (bottom-face net flux)", "T [K]")):
            a.set_xlabel("i (array index)")
            a.set_title("%s  %s  (last 6 cycles)" % (arm, t))
            a.legend(fontsize=7)
        ax[0].set_yscale("symlog", linthresh=1e3)
        fig.tight_layout()
        p = os.path.join(OUT, "profile_Fb_T_end_%s.png" % arm)
        fig.savefig(p, dpi=130)
        plt.close(fig)
        paths.append(p)

    # 2. root-cell budget vs cycle
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
    for a, arm in zip(ax, ("exp_g3", "si_g6")):
        rec = recs[arm]
        n_root, _, ir = roots[arm]
        if n_root is None:
            continue
        u = rec["u_ener"][:, ir]
        du = np.diff(u)
        cterm = (cond_div(rec, ir)/DX1[ir])*rec["dt"]
        a.plot(rec["cycle"][:-1], du, label="delta_u")
        a.plot(rec["cycle"], rec["rt_de"][:, ir], label="rt_de")
        a.plot(rec["cycle"], cterm, label="cond*dt")
        a.axvline(rec["cycle"][n_root], color="k", ls=":", lw=1)
        a.set_yscale("symlog", linthresh=1e-14)
        a.set_xlabel("cycle")
        a.set_title("%s, root cell i = %d" % (arm, ir))
        a.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(OUT, "root_budget.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths.append(p)
    return paths


def main():
    recs = {}
    roots = {}
    for arm in ARMS:
        rec = read_arm(arm)
        recs[arm] = rec
        print("%s: %d cyclediag files, cycles %d .. %d, t %.8g .. %.8g"
              % (arm, len(rec["cycle"]), rec["cycle"][0], rec["cycle"][-1],
                 rec["time"][0], rec["time"][-1]))
        for s in rec["skipped"]:
            print("  SKIPPED %s" % s)

    print("\n" + "#" * 72)
    print("# 1. ROOT CELL")
    print("#" * 72)
    maxT_table(recs["exp_g6"], "exp_g6", every=25)
    for arm in ("exp_g3", "si_g6"):
        rec = recs[arm]
        maxT_table(rec, arm, every=25)
        trigger_table(rec, arm)
        n, c, i = root_cell(rec)
        roots[arm] = (n, c, i)
        if n is None:
            print("\n%s: NO cycle has any active cell above 1.5x its T 20 cycles "
                  "earlier." % arm)
            continue
        r = (rec["rt_T"][n, i]/rec["rt_T"][n-20, i])
        print("\n%s ROOT: first cycle with T > 1.5 x T(20 cycles earlier) is cycle %d "
              "(index %d), i = %d, ratio %.4f (T %.5e -> %.5e)"
              % (arm, c, n, i, r, rec["rt_T"][n-20, i], rec["rt_T"][n, i]))
        root_table(rec, arm, i, n)

    print("\n" + "#" * 72)
    print("# 2. FLUX SANITY")
    print("#" * 72)
    for arm in ("exp_g3", "si_g6"):
        n, c, i = roots.get(arm, (None, None, None))
        if n is None:
            continue
        flux_sanity(recs[arm], arm, n)

    print("\n" + "#" * 72)
    print("# 2b. FLUX SANITY at the LAST CYCLE of each arm (supplementary)")
    print("#" * 72)
    for arm in ("exp_g3", "exp_g6", "si_g6"):
        rec = recs[arm]
        nlast = len(rec["cycle"]) - 1
        flux_sanity(rec, arm, nlast, only=(nlast,))

    print("\n" + "#" * 72)
    print("# 3. SIGN / LOCALITY (si_g6)")
    print("#" * 72)
    n, c, i = roots.get("si_g6", (None, None, None))
    if n is not None:
        sign_locality(recs["si_g6"], "si_g6", i, n)

    print("\n" + "#" * 72)
    print("# 3b. SIGN / LOCALITY at the END of each arm (supplementary), around the")
    print("#     hottest active cell of the column on the last cycle")
    print("#" * 72)
    for arm in ("exp_g3", "si_g6"):
        rec = recs[arm]
        nlast = len(rec["cycle"]) - 1
        i0, i1 = int(rec["is"]), int(rec["ie"])
        ih = i0 + int(np.argmax(rec["rt_T"][nlast, i0:i1+1]))
        print("\n%s: hottest active cell on the last cycle (%d) is i = %d, "
              "T = %.5e K" % (arm, rec["cycle"][nlast], ih, rec["rt_T"][nlast, ih]))
        sign_locality(rec, arm, ih, nlast - 6, nafter=6)

    print("\n" + "#" * 72)
    print("# 4. PLOTS")
    print("#" * 72)
    for p in plots(recs, roots):
        print("PLOT %s" % p)


if __name__ == "__main__":
    main()
