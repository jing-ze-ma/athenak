#!/usr/bin/env python3
"""prodchk.py (s3grid-1017): READ-ONLY surface/background check of an AG Car (he_star_m1) run, the BGDIAG metrics.

usage: python3 prodchk.py --run DIR [--stdout 'GLOB'] [--bins 'GLOB'] [--dstride 5e6] [--win 2e7] [--rph 2.7025e13]
                          [--dmax 2e-16] [--out PREFIX]
  DIR        run dir with <basename>.user.hst and <basename>.log (event log); bins default DIR/bin/*.hydro_w.*.bin
  --stdout   glob of the run's stdout logs (all links; used for BE share, DT_RETRY, NC, cycle->time map)
  --dstride  use one hydro_w dump per DSTRIDE seconds (default 5e6)
Writes PREFIX_hst.txt (per window), PREFIX_dump.txt (per dump), PREFIX_sum.txt; prints all.
Metrics (definitions as NOTE rol2u_safety_1015 BGDIAG):
 hst  : L_top/L_in, Mdot_top, P_sponge/L_top, budget residual Bres = Etot - int(L_in - L_top + P_esrc - P_sponge) dt
        (window slope erg/s), M_tot.
 log  : m1_pos_un (E-floor energy created, cumulative per link; resets at restart handled), m1_pos_flr, vceil_de, eos_tfloor,
        eos_efloor, eos_vceil (per-row counts).
 stdout: steps, share of steps with 'hesdirk2 stage 1 NOT ADMISSIBLE' (-> backward Euler), DT_RETRY cuts, NON-CONVERGED.
 dumps: r_surf = outermost r/R_ph with sin-weighted shell <rho> > 1e-12; steps = (first cell with rho < 1e-12 outside
        0.9 R_ph) - (last cell before it with rho > 1e-10); cliff fractions of columns with steps <= 1 and <= 2 (BGDIAG's
        '2-dex drop within <= 1 cell' = <= 2 steps, i.e. at most one intermediate cell); median cells across that
        2-dex drop; dr/r and g dr/(e/rho) at r_surf; M(>1.2,1.5,2 R_ph); cells |v_r| > 500 km/s; top 4 cells: max rho,
        fraction with rho < he_sponge_dmax (sponge active only there).
"""
import argparse
import glob
import gzip
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
GM = 6.674e-8 * 55.0 * 1.98847e33   # AG Car A point mass 55 Msun (only for the g dr diagnostic)


def opn(f):
    return gzip.open(f, 'rt', errors='replace') if f.endswith('.gz') else open(f, errors='replace')


def read_hst(f):
    hdr = None
    rows = []
    for ln in opn(f):
        if ln.startswith('#  ['):
            hdr = re.findall(r'\[\d+\]=(\S+)', ln)
            continue
        if ln.startswith('#') or hdr is None:
            continue
        v = ln.split()
        if len(v) == len(hdr):
            try:
                rows.append([float(x) for x in v])
            except ValueError:
                pass
    x = np.array(rows)
    # restarts can rewrite a stretch of time: keep a strictly increasing time series (last write wins)
    t = x[:, 0]
    keep = np.ones(len(t), bool)
    tmax = np.inf
    for k in range(len(t) - 1, -1, -1):
        if t[k] >= tmax:
            keep[k] = False
        else:
            tmax = t[k]
    x = x[keep]
    return {n: x[:, k] for k, n in enumerate(hdr)}


def read_evlog(f):
    names = None
    rows = []
    for ln in opn(f):
        if ln.startswith('#  cycle') or ln.startswith('# cycle'):
            names = ln[1:].split()
            continue
        if ln.startswith('#') or names is None:
            continue
        v = ln.split()
        if len(v) == len(names):
            rows.append([float(x) for x in v])
    x = np.array(rows)
    return {n: x[:, k] for k, n in enumerate(names)}


def incr_cum(c):
    """increments of a per-link cumulative counter that resets to ~0 at each restart"""
    d = np.diff(c, prepend=0.0)
    r = d < 0
    d[r] = c[r]
    return d


def read_stdout(globs):
    cyc_t = {}
    adm = set()
    retr = []
    nc = []
    fat = 0
    fs = sorted(set(sum([glob.glob(g) for g in globs], [])))
    for f in fs:
        for ln in opn(f):
            if ln.startswith('elapsed='):
                m = re.search(r'cycle=(\d+) time=(\S+)', ln)
                if m:
                    cyc_t[int(m.group(1))] = float(m.group(2))
            elif 'NOT ADMISSIBLE at cycle' in ln and 'stage 1' in ln:
                adm.add(int(re.search(r'at cycle (\d+)', ln).group(1)))
            elif ln.startswith('DT_RETRY cycle='):
                m = re.search(r'cycle=(\d+) t=(\S+).*reason=(\S+)', ln)
                if m:
                    retr.append((int(m.group(1)), float(m.group(2)), m.group(3)))
            elif 'NON-CONVERGED after' in ln:
                m = re.search(r'cycle[= ](\d+)', ln)
                nc.append(int(m.group(1)) if m else -1)
            elif 'FATAL' in ln:
                fat += 1
    return fs, cyc_t, adm, retr, nc, fat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True)
    ap.add_argument('--stdout', action='append', default=[])
    ap.add_argument('--bins', default=None)
    ap.add_argument('--hst', default=None)
    ap.add_argument('--evlog', default=None)
    ap.add_argument('--dstride', type=float, default=5e6)
    ap.add_argument('--win', type=float, default=2e7)
    ap.add_argument('--rph', type=float, default=2.7025e13)
    ap.add_argument('--dmax', type=float, default=2e-16)
    ap.add_argument('--out', default='prodchk')
    ap.add_argument('--nodumps', action='store_true')
    a = ap.parse_args()
    out = {k: open('%s_%s.txt' % (a.out, k), 'w') for k in ('hst', 'dump', 'sum')}

    def P(k, s):
        print(s)
        out[k].write(s + '\n')
        out[k].flush()

    hf = a.hst or (glob.glob(os.path.join(a.run, '*.user.hst')) + [None])[0]
    lf = a.evlog or (glob.glob(os.path.join(a.run, '*.log')) + [None])[0]
    H = read_hst(hf)
    t = H['time']
    P('sum', 'run %s\n  user.hst %s: %d rows, t %.4e -> %.4e' % (a.run, hf, len(t), t[0], t[-1]))
    src = H['L_in'] - H['L_top'] + H.get('P_esrc', 0 * t) - H.get('P_sponge', 0 * t)
    S = np.concatenate([[0.0], np.cumsum(0.5 * (src[1:] + src[:-1]) * np.diff(t))])
    B = (H['Etot'] - H['Etot'][0]) - S
    # stdout
    fs, cyc_t, adm, retr, nc, fat = read_stdout(a.stdout) if a.stdout else ([], {}, set(), [], [], 0)
    cy = np.array(sorted(cyc_t))
    ct = np.array([cyc_t[c] for c in cy]) if len(cy) else np.array([])
    P('sum', '  stdout files %d, step lines %d, NOT-ADMISSIBLE cycles %d, DT_RETRY lines %d, NC lines %d, FATAL lines %d'
      % (len(fs), len(cy), len(adm), len(retr), len(nc), fat))
    # event log
    E = None
    if lf and os.path.exists(lf):
        E = read_evlog(lf)
        ecy = E['cycle']
        if len(cy):
            et = np.interp(ecy, cy, ct, left=np.nan, right=np.nan)
        else:   # fallback: event-log rows every output dt (1e5 s) from the first hst time
            et = t[0] + 1e5 * np.arange(1, len(ecy) + 1)
        dun = incr_cum(E['m1_pos_un']) if 'm1_pos_un' in E else None
        dfl = incr_cum(E['m1_pos_flr']) if 'm1_pos_flr' in E else None
        P('sum', '  event log %s: %d rows (cycle %d -> %d)' % (lf, len(ecy), ecy[0], ecy[-1]))
    P('hst', '# window [t0,t1] | L_top/L_in  Mdot_top g/s  P_spg/L_top | Bres slope erg/s (/L_top) | logged m1_pos_un rate erg/s'
      ' | m1_pos_flr  vceil_de  tfloor efloor vceil (per window) | steps  BE share  DT_RETRY  NC')
    ed = np.arange(np.floor(t[0] / a.win) * a.win, t[-1] + a.win, a.win)
    for t0, t1 in zip(ed[:-1], ed[1:]):
        w = (t >= t0) & (t < t1)
        if w.sum() < 3:
            continue
        i0, i1 = np.where(w)[0][[0, -1]]
        dt_ = t[i1] - t[i0]
        bs = (B[i1] - B[i0]) / dt_ if dt_ > 0 else np.nan
        lt = H['L_top'][w].mean()
        s = '%.3e %.3e | %.3f %9.2e %8.1e | %+9.2e (%+.1e)' % (t0, t1, lt / H['L_in'][w].mean(), H['Mdot_top'][w].mean(),
                                                               H.get('P_sponge', 0 * t)[w].mean() / lt, bs, bs / lt)
        if E is not None:
            we = (et >= t0) & (et < t1)
            if we.any():
                s += ' | %9.2e | %9.3e %9.3e %7d %7d %6d' % (
                    dun[we].sum() / (t1 - t0) if dun is not None else np.nan, dfl[we].sum() if dfl is not None else np.nan,
                    E['vceil_de'][we].sum(), E['eos_tfloor'][we].sum(), E['eos_efloor'][we].sum(), E['eos_vceil'][we].sum())
            else:
                s += ' | (no log rows)'
        if len(cy):
            wc = (ct >= t0) & (ct < t1)
            n = wc.sum()
            cs = set(cy[wc].tolist())
            nb = len(cs & adm)
            nr = sum(1 for c, tt, r in retr if t0 <= tt < t1 and 'recovered' not in r)
            nn = sum(1 for c in nc if c in cs)
            s += ' | %7d %.3f %6d %5d' % (n, nb / max(n, 1), nr, nn)
        P('hst', s)
    P('sum', '  whole run: L_top/L_in first/last 2e6 s: %.3f / %.3f; Bres total %.3e erg (rate %.2e erg/s); Mdot_top last 2e6 s %.2e'
      % (H['L_top'][t < t[0] + 2e6].mean() / H['L_in'][t < t[0] + 2e6].mean(), H['L_top'][t > t[-1] - 2e6].mean() /
         H['L_in'][t > t[-1] - 2e6].mean(), B[-1], B[-1] / (t[-1] - t[0]), H['Mdot_top'][t > t[-1] - 2e6].mean()))
    if E is not None and dun is not None:
        P('sum', '  logged m1_pos_un total %.3e erg -> Bres / logged = %.3f' % (np.nansum(dun), B[-1] / max(np.nansum(dun), 1e-30)))
    if a.nodumps:
        return
    import bin_convert as bc
    fb = sorted(glob.glob(a.bins or os.path.join(a.run, 'bin', '*.hydro_w.*.bin')))
    sel, tl = [], -1e30
    for f in fb:
        n = int(re.search(r'\.(\d+)\.bin$', f).group(1))
        sel.append((n, f))
    # stride by dump number using the first two dumps' times
    P('dump', '#      t   r_surf | frac(2dex in 1 step) frac(<=2 steps, BGDIAG def)  med steps/2dex  ncol | dr/r@rs  g*dr/(e/rho)@rs (median col) | M>1.2 M>1.5 M>2 (g)'
      ' | n|vr|>500 max vr | top4: max rho  frac<dmax | <v_r> 1.05 1.10 km/s')
    geo = None
    for n, f in sel:
        h = bc.read_binary_as_athdf(f)
        tt = float(h['Time'])
        if tt < tl + a.dstride * (1 - 1e-6):
            continue
        tl = tt
        rho = np.asarray(h['dens'], float)
        if geo is None:
            r = np.asarray(h['x1v'], float)
            xf = np.asarray(h['x1f'], float)
            th = np.asarray(h['x2v'], float)
            x2f = np.asarray(h['x2f'], float)
            x3f = np.asarray(h['x3f'], float)
            dV = ((xf[1:] ** 3 - xf[:-1] ** 3) / 3)[None, None, :] * (np.cos(x2f[:-1]) - np.cos(x2f[1:]))[None, :, None] * \
                np.diff(x3f)[:, None, None]
            wsin = np.sin(th)[None, :, None] / np.sin(th).sum() / rho.shape[0]
            x = r / a.rph
            dr = np.diff(xf)
            i09 = int(np.searchsorted(x, 0.9))
            geo = 1
        e = np.asarray(h['eint'], float)
        v = np.asarray(h['velx'], float)
        rm = (rho * wsin).sum(axis=(0, 1))
        j = np.where(rm > 1e-12)[0]
        isf = j.max() if len(j) else -1
        rs = x[isf]
        # per column: last cell with rho > 1e-10 (from 0.9 R_ph outward, before the first cell < 1e-12)
        R = rho[:, :, i09:]
        lo = R < 1e-12
        ib = np.where(lo.any(axis=2), lo.argmax(axis=2), -1)          # first cell < 1e-12
        hi = R > 1e-10
        nk, nj, ni = R.shape
        idx = np.arange(ni)[None, None, :]
        ia = np.where(hi & (idx < np.where(ib >= 0, ib, ni)[:, :, None]), idx, -1).max(axis=2)
        ok = (ib >= 0) & (ia >= 0)
        ncell = (ib - ia)[ok]
        cl = (ncell <= 1).mean() if ok.any() else np.nan
        cl2 = (ncell <= 2).mean() if ok.any() else np.nan
        # g dr / (e/rho) at the column's surface cell (first < 1e-12 minus 1)
        isurf = (ia + i09)[ok]
        kk, jj = np.where(ok)
        gdr = GM / r[isurf] ** 2 * dr[isurf] / (e[kk, jj, isurf] / rho[kk, jj, isurf])
        M = [(rho * dV)[:, :, x > q].sum() for q in (1.2, 1.5, 2.0)]
        wv = x > 1.2
        nf = int((np.abs(v[:, :, wv]) > 5e7).sum())
        top = rho[:, :, -4:]
        vm = (v * wsin).sum(axis=(0, 1))
        P('dump', '%.4e %6.3f | %.3f %.3f %6.1f %6d | %.4f %8.3f | %.2e %.2e %.2e | %6d %7.1f | %.2e %.3f | %+.2f %+.2f' % (
            tt, rs, cl, cl2, np.median(ncell) if ok.any() else np.nan, ok.sum(), dr[isf] / r[isf], np.median(gdr) if len(gdr) else np.nan,
            M[0], M[1], M[2], nf, v[:, :, wv].max() / 1e5, top.max(), (top < a.dmax).mean(),
            vm[np.argmin(abs(x - 1.05))] / 1e5, vm[np.argmin(abs(x - 1.10))] / 1e5))


if __name__ == '__main__':
    main()
