#!/usr/bin/env python3
"""rs_eval.py -- ACCURACY of an IDORT radiative-shock run (rsx_m<M>, or the old rs_m<M>
decks) against the Lowrie & Edwards (2008) semi-analytic solution of le08_shock.py
(deck-implied P0/sigma_a; 'nominal' for the old t7-based rs_m decks is NOT used: the
exact physical problem of the deck is the reference in both cases).

Metrics on the LAST tab dump (and on the first, = IC fidelity), x2/x3 slice:
  vshock   shock speed over the last output interval / upstream speed v0
  drift    shock position (adiabatic sonic point M = v/a = 1, sub-cell) minus the deck's
           m1_shock_xs, in cells (the shock is held in its own frame: must stay put)
  L1 rho/Tg/Tr  sum|q - q_ref| / sum|q_ref| over all cells except +-3 cells of the shock,
           reference shifted onto the measured shock (as shockcmp.py / t7)
  L1pre    sum|Tg - Tg_ref| / sum|Tg_ref - 1| over the precursor (x < x_shock - 3 dx):
           precursor SHAPE error relative to the precursor heating
  Lpre(f)  precursor length: distance from the shock to where Tg = 1 + f (Tp_ref - 1),
           f = 0.1 and 0.01 (sub-cell interpolation), vs the exact value
  Trs      T_rad at the shock (continuous across it) vs exact
  Tmax     max Tg vs the exact spike maximum AND vs the exact profile averaged over the
           code's cells (the spike is sub-cell for M0 >= 5 at any affordable resolution)
  down     mean rho, Tg, Tr over the downstream half (x_shock .. x1max - 3 dx) vs the
           exact profile there (= far-field state for M0 >= 5)
  mflux    max |rho v / (rho0 v0) - 1| over the interior, > 10 cells from the shock
           (mflux_near: 3..10 cells from it = post-shock ringing)
usage: python rs_eval.py <rundir> [<rundir> ...]   (rundir has in.athinput + tab/)
       env RS_CACHE = directory for the cached semi-analytic solutions (default refs/)
"""
import glob
import os
import pickle
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import le08_shock as L  # noqa: E402

GAM = 5.0 / 3.0
TUN = 8.2499208e7
ARC = 1.6332399e-46
CACHE = os.environ.get("RS_CACHE", "/raven/ptmp/jinma/idortrs_1017/refs")
# IDORT S_N transport reference (Zenodo 15032138), nominal LE08 parameters, M0 = 1.2/2/3
TREF = os.environ.get("RS_TREF", "/raven/ptmp/jinma/suite_1011/ext/idort_mini/x/IDORT_mini/"
                      "shock_1d/Rad_shock_sol")
TLAB = {1.2: "M1p2", 2.0: "M2", 3.0: "M3"}


def transport(m0):
    """(x - x_sonic, rho, Tg, Tr) of the IDORT transport solution, or None."""
    lab = TLAB.get(m0)
    if lab is None or not os.path.isdir(os.path.join(TREF, lab)):
        return None
    r = {k: np.loadtxt(os.path.join(TREF, lab, "%s_4096.dat" % k))
         for k in ("x", "rho", "p", "E", "u")}
    x, rho = r["x"], r["rho"]
    tg = GAM * r["p"] / rho
    tr = np.maximum(r["E"], 0.0) ** 0.25
    mach = np.abs(r["u"]) / np.sqrt(tg)
    k = np.nonzero((mach[:-1] >= 1.0) & (mach[1:] < 1.0))[0][0]
    xc = x[k] + (mach[k] - 1.0) / (mach[k] - mach[k + 1]) * (x[k + 1] - x[k])
    return x - xc, rho, tg, tr


def solution(m0):
    f = os.path.join(CACHE, "sol_m%g.pkl" % m0)
    if os.path.exists(f):
        return pickle.load(open(f, "rb"))
    p0, sa, ka = L.implied_params()
    s = L.Solution(L.LE08(m0, p0=p0, sigma_a=sa, kappa=ka))
    s.raw = None
    try:
        pickle.dump(s, open(f, "wb"))
    except OSError:
        pass
    return s


def deckval(txt, key):
    mm = re.search(r"^\s*%s\s*=\s*(\S+)" % re.escape(key), txt, re.M)
    return mm.group(1) if mm else None


def rdtab(f):
    lines = open(f).readlines()
    t = float(lines[0].split("time=")[1].split()[0])
    hdr = lines[1].lstrip("#").split()
    a = np.loadtxt(f)
    if len(hdr) != a.shape[1]:   # slice outputs list the sliced coordinates in the header
        hdr = [c for c in hdr if c not in ("k", "x3v", "j", "x2v")]
    return t, {h: a[:, i] for i, h in enumerate(hdr)}


def profile(d, which):
    fs = sorted(glob.glob(os.path.join(d, "tab", "*.m1.*.tab")))
    f = fs[which]
    t, m = rdtab(f)
    _, h = rdtab(f.replace(".m1.", ".hydro_w."))
    return t, h["x1v"], h["dens"], h["velx"], h["eint"], m["m1_e"]


def cellavg(sol, xe):
    """exact profile averaged over cells with edges xe (nondim), 64 sub-points."""
    q = np.linspace(0.0, 1.0, 65)[:-1] + 1.0 / 128.0
    xx = (xe[:-1, None] + (xe[1:] - xe[:-1])[:, None] * q[None, :]).ravel()
    s = sol.sample(xx)
    return s.reshape(5, xe.size - 1, q.size).mean(axis=2)


def evaluate(d, which=-1):
    dk = d if os.path.exists(os.path.join(d, "in.athinput")) else os.path.dirname(d)
    txt = open(os.path.join(dk, "in.athinput")).read()
    try:
        args = open(os.path.join(d, "args.txt")).read()
    except OSError:
        args = ""
    vl = float(deckval(txt, "m1_shock_v_l"))
    rhol = float(deckval(txt, "m1_shock_rho_l"))
    p0, sa, ka = L.implied_params()
    u = L.Units(L.LE08(2.0, p0=p0, sigma_a=sa, kappa=ka))
    m0 = round(vl / u.v, 4)
    sol = solution(m0)
    mt = sol.metrics()
    x1min = float(deckval(txt, "x1min"))
    x1max = float(deckval(txt, "x1max"))
    nx1 = int(deckval(txt, "nx1"))
    xs_deck = float(deckval(txt, "m1_shock_xs"))
    t, x, rho, v, eint, er = profile(d, which)
    nx1 = np.unique(x).size     # the deck's nx1 may be overridden on the command line
    xn = (x - xs_deck) / u.L                    # nondim, deck shock at 0
    dxn = (x1max - x1min) / nx1 / u.L
    rhon = rho / rhol
    tg = (GAM - 1.0) * eint / (rho * TUN) / L.T0
    tr = (np.maximum(er, 0.0) / ARC) ** 0.25 / TUN / L.T0
    vn = v / u.v
    mach = vn / np.sqrt(tg)
    # shock = first downward crossing of M = 1
    k = np.nonzero((mach[:-1] >= 1.0) & (mach[1:] < 1.0))[0]
    if k.size == 0:
        return {"dir": d, "t": t, "error": "no sonic crossing"}
    k = k[0]
    w = (mach[k] - 1.0) / (mach[k] - mach[k + 1])
    xc = xn[k] + w * (xn[k + 1] - xn[k])
    out = {"dir": d, "t": t, "m0": m0, "nx1": nx1, "kind": mt["kind"],
           "drift": xc / dxn, "dx": dxn}
    # shock SPEED over the last output interval, in units of the upstream speed: the
    # relocation itself is the transient from the diffusion IC to the code's (transport)
    # steady state; a steady scheme must end with |v_shock| / v0 << 1
    nf = len(glob.glob(os.path.join(d, "tab", "*.m1.*.tab")))
    if which == -1 and nf >= 3:
        t2, x2, rho2, v2, e2, _ = profile(d, -2)
        tg2 = (GAM - 1.0) * e2 / (rho2 * TUN) / L.T0
        m2 = (v2 / u.v) / np.sqrt(tg2)
        x2n = (x2 - xs_deck) / u.L
        k2 = np.nonzero((m2[:-1] >= 1.0) & (m2[1:] < 1.0))[0]
        if k2.size:
            k2 = k2[0]
            xc2 = x2n[k2] + (m2[k2] - 1.0) / (m2[k2] - m2[k2 + 1]) * (x2n[k2 + 1] - x2n[k2])
            out["vshock"] = float((xc - xc2) * u.L / ((t - t2) * vl))
    xe = np.concatenate([[xn[0] - 0.5 * dxn], 0.5 * (xn[1:] + xn[:-1]), [xn[-1] + 0.5 * dxn]])
    ref = cellavg(sol, xe - xc)
    away = np.abs(xn - xc) > 3.0 * dxn
    for nm, q, r in (("rho", rhon, ref[0]), ("Tg", tg, ref[2]), ("Tr", tr, ref[3])):
        out["L1_" + nm] = float(np.sum(np.abs(q[away] - r[away])) / np.sum(np.abs(r[away])))
    pre = xn < xc - 3.0 * dxn
    out["L1pre"] = float(np.sum(np.abs(tg[pre] - ref[2][pre])) /
                         max(np.sum(np.abs(ref[2][pre] - 1.0)), 1e-300))
    for f in (0.1, 0.01):
        lev = 1.0 + f * (mt["Tp"] - 1.0)
        up = np.nonzero((xn < xc) & (tg >= lev))[0]
        if up.size and up[0] > 0:
            j = up[0]
            wl = (lev - tg[j - 1]) / (tg[j] - tg[j - 1])
            xl = xn[j - 1] + wl * (xn[j] - xn[j - 1])
            out["Lpre_%g" % f] = float(xc - xl)
        else:
            out["Lpre_%g" % f] = float("nan")
        out["Lpre_%g_ref" % f] = mt["Lpre_%g" % f]
        out["eLpre_%g" % f] = out["Lpre_%g" % f] / mt["Lpre_%g" % f] - 1.0
    out["Trs"] = float(np.interp(xc, xn, tr))
    out["Trs_ref"] = mt["Trs"]
    out["eTrs"] = out["Trs"] / mt["Trs"] - 1.0
    out["Tmax"] = float(tg.max())
    out["Tmax_ref"] = mt["Tmax"]
    out["Tmax_refcell"] = float(ref[2].max())
    out["eTmax"] = out["Tmax"] / mt["Tmax"] - 1.0
    out["eTmax_cell"] = out["Tmax"] / out["Tmax_refcell"] - 1.0
    out["spike_fwhm_ref_cells"] = mt["spike_fwhm"] / dxn
    dn = (xn > xc + 0.5 * (xn[-1] - xc)) & (xn < xn[-1] - 3.0 * dxn)
    for nm, q, r in (("rho", rhon, ref[0]), ("Tg", tg, ref[2]), ("Tr", tr, ref[3])):
        out["edown_" + nm] = float(q[dn].mean() / r[dn].mean() - 1.0)
    tt = transport(m0)
    if tt is not None:
        # the same L1 / Trs / Lpre against the TRANSPORT solution (on the code's cells)
        for nm, q, j in (("rho", rhon, 1), ("Tg", tg, 2), ("Tr", tr, 3)):
            rr = np.interp(xn[away] - xc, tt[0], tt[j])
            out["T_L1_" + nm] = float(np.sum(np.abs(q[away] - rr)) / np.sum(np.abs(rr)))
        out["T_eTrs"] = out["Trs"] / float(np.interp(0.0, tt[0], tt[3])) - 1.0
        up = tt[0] < 0.0
        for f in (0.1, 0.01):
            lev = 1.0 + f * (mt["Tp"] - 1.0)
            j = np.nonzero(up & (tt[2] >= lev))[0][0]
            xl = tt[0][j - 1] + (lev - tt[2][j - 1]) / (tt[2][j] - tt[2][j - 1]) * \
                (tt[0][j] - tt[0][j - 1])
            out["T_eLpre_%g" % f] = out["Lpre_%g" % f] / (-xl) - 1.0
        out["T_Tmax"] = float(tt[2].max())
    ii = np.arange(xn.size)
    far = np.abs(xn - xc) > 10.0 * dxn
    inner = (ii > 3) & (ii < xn.size - 4)
    mf = np.abs(rhon * vn / m0 - 1.0)
    out["mflux"] = float(np.max(mf[inner & far]))
    # post-shock oscillation (odd-even ringing of a slowly moving captured shock)
    out["mflux_near"] = float(np.max(mf[inner & ~far & away]))
    return out


def fmt(o):
    if "error" in o:
        return "%s ERROR %s" % (o["dir"], o["error"])
    return ("%-34s t=%.3e M%-5g %-13s n%-5d drift %+6.2f c vs %+.1e | L1 rho %.4f Tg %.4f "
            "Tr %.4f pre %.4f | "
            "Lpre0.1 %+6.2f%% Lpre0.01 %+6.2f%% | Trs %+6.3f%% | Tmax %.4f ref %.4f (cell %.4f; "
            "fwhm %.2g c) | down rho %+.4f%% Tg %+.4f%% Tr %+.4f%% | mflux %.1e near %.1e"
            % (o["dir"][-34:], o["t"], o["m0"], o["kind"], o["nx1"], o["drift"],
               o.get("vshock", float("nan")), o["L1_rho"],
               o["L1_Tg"], o["L1_Tr"], o["L1pre"], 100 * o["eLpre_0.1"], 100 * o["eLpre_0.01"],
               100 * o["eTrs"], o["Tmax"], o["Tmax_ref"], o["Tmax_refcell"],
               o["spike_fwhm_ref_cells"], 100 * o["edown_rho"], 100 * o["edown_Tg"],
               100 * o["edown_Tr"], o["mflux"], o["mflux_near"])
            + ("" if "T_L1_Tg" not in o else
               " || vs TRANSPORT: L1 rho %.4f Tg %.4f Tr %.4f | Lpre0.1 %+.2f%% 0.01 %+.2f%% | "
               "Trs %+.3f%% | Tmax ref %.4f" % (o["T_L1_rho"], o["T_L1_Tg"], o["T_L1_Tr"],
                                             100 * o["T_eLpre_0.1"], 100 * o["T_eLpre_0.01"],
                                             100 * o["T_eTrs"], o["T_Tmax"])))


if __name__ == "__main__":
    import json
    res = []
    for d in sys.argv[1:]:
        for which, tag in ((0, "t0 "), (-1, "end")):
            try:
                o = evaluate(d.rstrip("/"), which)
            except Exception as e:  # noqa: BLE001
                o = {"dir": d, "error": repr(e)}
            o["tag"] = tag
            res.append(o)
            print(tag, fmt(o))
            sys.stdout.flush()
    jf = os.environ.get("RS_JSON")
    if jf:
        json.dump(res, open(jf, "w"), indent=1)
