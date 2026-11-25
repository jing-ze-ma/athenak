#!/usr/bin/env python3
"""zzzzh4 / zzzzh5 (READ ONLY).  Per-row energy and mass closure of the production hst files.
E: R = Etot(c) - Etot(c-1) - dt (L_in - L_top + P_esrc - P_sponge - E_hflx);  M: RM = M_tot(c) - M_tot(c-1) + dt M_hflx.
dt = time(c) - time(c-1); rates at row c (the step that ended there), also tested at row c-1 and the mean.  Baseline = running median (51 rows)
inside a chunk; outliers |R - baseline| > 10 x MAD (MAD = median |R - baseline| over the star's consecutive rows)."""
import re, sys, glob, os, numpy as np
from numpy.lib.stride_tricks import sliding_window_view as sw
J = '/resnick/groups/carnegie_poc/jingze/'
STARS = [('bsg', 'prod_c1_1124/bsg', 'bsg3d', []), ('t9l62', 'prod_c4_1124/t9l62', 't9l62', [3713]), ('m35mms', 'prod_c1_1124/m35mms', 'm35mms', [16340]),
         ('m35zams', 'prod_c1_1124/m35zams', 'm35zams', [50816]), ('t9mlt', 'prod_mlt_1125/t9mlt', 't9l62', [])]
OUT = J + 'eres_1124/out/'; os.makedirs(OUT, exist_ok=True)
def rhst(f):
    """-> list of chunks (names, array); a header line starts a new chunk"""
    ch = []; n = None; rows = []
    def flush():
        if rows: ch.append((n, np.array(rows)))
    for ln in open(f, errors='replace'):
        if ln.startswith('#'):
            if '[1]=' in ln:
                flush(); rows = []; n = [m.split('=', 1)[1] for m in re.findall(r'\[\d+\]=\S+', ln)]
        elif n and len(ln.split()) == len(n):
            rows.append([float(x) for x in ln.split()])
    flush(); return ch
def rmed(x, w=51):
    h = w // 2; p = np.pad(x, (h, h), mode='reflect') if len(x) > h else np.pad(x, (h, h), mode='edge'); return np.median(sw(p, w), axis=1)
def logs(d):
    """per-cycle counters from log.* (plain read); cumulative floor_audit / fofc counters differenced inside each log file (process)"""
    C = {}; files = sorted(glob.glob(J + d + '/log.[0-9]*'), key=lambda f: int(f.split('.')[-1])); cyct = []
    def put(c, k, v): C.setdefault(c, {})[k] = C.get(c, {}).get(k, 0) + v if k in ('fb_pos_gas', 'fcap_v', 'fcap_T', 'c2p_dfloor', 'floor_pay', 'auto_be', 'e_floor') else v
    for li, f in enumerate(files):
        prev = {}
        for ln in open(f, errors='replace'):
            if ln.startswith('CLEAN FALLBACK cycle'):
                m = re.match(r'CLEAN FALLBACK cycle (\d+): pos_gas (\S+) .*? e_floor (\S+) .*?auto_be (\S+) floorcap v (\S+) T (\S+) c2p_dfloor (\S+) \(mass (\S+)\) floor_pay (\S+) \(unpaid', ln)
                if m:
                    c = int(m.group(1)); [put(c, k, float(v.rstrip(','))) for k, v in zip(('fb_pos_gas', 'e_floor', 'auto_be', 'fcap_v', 'fcap_T', 'c2p_dfloor', 'c2p_dfloor_mass', 'floor_pay'), m.groups()[1:])]
            elif 'floor_audit cycle' in ln:
                m = re.search(r'floor_audit cycle (\d+) \[(update|refresh)\] cell-c2p (\S+) .*?floored (\S+) .*?skip\(dfloor/vceil\) (\S+)', ln)
                if m:
                    c = int(m.group(1)); cur = [float(m.group(i)) for i in (3, 4, 5)]; p = prev.get(m.group(2))
                    if p and cur[0] >= p[0]: put(c, 'floored_' + m.group(2), cur[1] - p[1]); put(c, 'skip_' + m.group(2), cur[2] - p[2])
                    prev[m.group(2)] = cur
                m = re.search(r'floor_audit cycle (\d+) pay: paid (\S+) unpaid (\S+) .*?fofc_retest: retested (\S+) unflagged (\S+)', ln)
                if m:
                    c = int(m.group(1)); cur = [float(m.group(i)) for i in (2, 4, 5)]; p = prev.get('pay')
                    if p and cur[0] >= p[0]: put(c, 'paid', cur[0] - p[0]); put(c, 'fofc_retested', cur[1] - p[1]); put(c, 'fofc_unflagged', cur[2] - p[2])
                    prev['pay'] = cur
            elif 'fofc_report' in ln:
                m = re.search(r'fofc_report t=\S+ cycle=(\d+) stages=(\d+) fired cell-stages: at-floor=(\S+) other=(\d+)', ln)
                if m: put(int(m.group(1)), 'fofc_rep_at', float(m.group(3).rstrip(';'))); put(int(m.group(1)), 'fofc_rep_other', float(m.group(4)))
            elif 'cycle=' in ln and 'time=' in ln and 'elapsed' in ln:
                m = re.search(r'cycle=(\d+) time=([0-9.eE+]+)', ln)
                if m: cyct.append((li, int(m.group(1)), float(m.group(2))))
    return C, cyct
P = []; pr = lambda *a: (P.append(' '.join(str(x) for x in a)), print(*a, flush=True))
for star, d, base, switch in STARS:
    ch = rhst(J + d + '/' + base + '.user.hst'); C, cyct = logs(d)
    pr('\n##### %s (%s): chunks %d; switch cycles named by the coordinator: %s; log cycle prints %d; counters for %d cycles' % (star, d, len(ch), switch or 'none', len(cyct), len(C)))
    ref = [(k, c, t) for k, c, t in cyct]
    # cycle of row 0 of each chunk: chunk 0 row 0 = cycle 0; later chunks: anchor on the first log line of that chunk
    start = [0] * len(ch)
    for k in range(len(ch)):
        t = ch[k][1][:, 0]; q = [(c, tt) for kk, c, tt in ref if kk == k]
        if q:
            c, tt = q[min(2, len(q) - 1)]; start[k] = c - int(np.argmin(abs(t - tt)))
    RE, RM, RM2, dtv, meta, cyc, Etot, M, rn, Ral = [], [], [], [], [], [], [], [], [], []
    for k, (n, a) in enumerate(ch):
        g = lambda nm: a[:, n.index(nm)]
        t = g('time'); E = g('Etot'); m = g('M_tot'); dt_col = g('dt'); dE = np.diff(E); dt = np.diff(t)
        rate = g('L_in') - g('L_top') + g('P_esrc') - g('P_sponge') - g('E_hflx')
        mh = g('M_hflx'); r1 = dE - dt * rate[1:]; r0 = dE - dt * rate[:-1]; rm = dE - dt * 0.5 * (rate[1:] + rate[:-1])
        var = {'row c': r1, 'row c-1': r0, 'mean': rm}
        cons = (dt > 0) & (np.abs(dt / dt_col[1:] - 1) < 1e-6)                                  # dt column = the step that ended at the row
        consp = (dt > 0) & (np.abs(dt / dt_col[:-1] - 1) < 1e-6)
        meta.append((k, len(a), int((dt <= 0).sum()), int(cons.sum()), int(consp.sum()), float(dt[dt > 0].min()) if (dt > 0).any() else 0))
        C0 = np.arange(1, len(a)) + start[k]
        RE.append(r1); RM.append(np.diff(m) + dt * mh[1:]); RM2.append(np.diff(m) + dt * (g('Mdot_top')[1:] - g('Mdot_bot')[1:])); cyc.append(C0); Etot.append(E[1:]); M.append(m[1:]); rn.append((cons, consp, {kk: vv for kk, vv in var.items()}, dt)); Ral.append((t[1:]))
        pr('  chunk %d: rows %d, t %.6e - %.6e, first cycle %d; row steps with dt <= 0: %d; rows with dt(time) = dt column of the same row: %d, of the previous row: %d of %d' % (k, len(a), t[0], t[-1], start[k], (dt <= 0).sum(), cons.sum(), consp.sum(), len(dt)))
    pr('  Etot ~ %.4e erg, eps-scale 1.1e-16 x Etot = %.2e erg; printed %%24.16e = 17 significant digits (round-trips a double), so the print is not the limit; the limit is the double sum over the cells' % (np.median(np.concatenate(Etot)), 1.1e-16 * np.median(np.concatenate(Etot))))
    # choose rate alignment by MAD on the interior rows (chunk 0)
    for k, (n, a) in enumerate(ch):
        pass
    best = None
    for al in ('row c', 'row c-1', 'mean'):
        rr = []
        for k in range(len(ch)):
            v = rn[k][2][al]; E_ = Etot[k]; q = rn[k][0] | rn[k][1]
            x = (v / E_)[q]; 
            if len(x) > 60: rr.append(np.abs(x - rmed(x)))
        mad = np.median(np.concatenate(rr)) if rr else np.nan
        pr('  rate at %-8s: MAD of (R - running median)/Etot = %.3e' % (al, mad)); 
        if best is None or mad < best[0]: best = (mad, al)
    al = 'row c'; pr('  -> rate taken at: row c for every star (source: dt and L_top columns are the step that ended at the row); the MAD of 0 at row c-1 for t9l62 is a running-median artefact on a smooth monotone series, not a better closure')
    X = []; XM = []
    for k in range(len(ch)):
        cons, consp, var, dt = rn[k]; q = cons | consp; v = var[al][q]; E_ = Etot[k][q]; x = v / E_
        xm = (RM[k] / M[k])[q]; xm2 = (RM2[k] / M[k])[q]; c = cyc[k][q]; t = Ral[k][q]
        b = rmed(x); bm = rmed(xm); X.append((c, t, x, x - b, np.median(abs(x - b)))); XM.append((c, t, xm, xm - bm, xm2))
    allx = np.concatenate([x[3] for x in X]); allm = np.concatenate([x[3] for x in XM]); allc = np.concatenate([x[0] for x in X])
    mad = np.median(abs(allx)); madm = np.median(abs(allm))
    out = np.where(abs(allx) > 10 * mad)[0]; om = np.where(abs(allm) > 10 * madm)[0]
    pr('  ENERGY rows used %d (consecutive), MAD %.3e of Etot; median(R/Etot) %.3e; outliers (10 MAD) %d; of them > 1e-10 of Etot: %d; largest |R - med|/Etot %.3e' % (len(allx), mad, np.median(np.concatenate([x[2] for x in X])), len(out), int((abs(allx[out]) > 1e-10).sum()), abs(allx).max()))
    pr('  MASS  rows used %d, MAD of (dM + dt M_hflx)/M_tot %.3e; median %.3e; outliers (10 MAD) %d; of them > 1e-14: %d, > 1e-12: %d; largest %.3e | variant dM + dt(Mdot_top - Mdot_bot): MAD %.3e, largest %.3e' % (len(allm), madm, np.median(np.concatenate([x[2] for x in XM])), len(om), int((abs(allm[om]) > 1e-14).sum()), int((abs(allm[om]) > 1e-12).sum()), abs(allm).max(), np.median(abs(np.concatenate([x[4] for x in XM]) - rmed(np.concatenate([x[4] for x in XM])))), abs(np.concatenate([x[4] for x in XM])).max()))
    keys = sorted({k for v in C.values() for k in v})
    def row(i, arr):
        c = int(allc[i]); s = 'c%-7d t %.6e %s=%+.3e' % (c, np.concatenate([x[1] for x in X])[i], arr, (allx if arr == 'dE/E' else allm)[i])
        for lag in (-1, 0, 1):
            u = C.get(c + lag, {}); s += ' |lag%+d ' % lag + ' '.join('%s=%g' % (k, u[k]) for k in ('fofc_retested', 'fofc_rep_other', 'c2p_dfloor', 'fcap_v', 'fb_pos_gas', 'floor_pay', 'paid', 'floored_update', 'skip_update', 'auto_be', 'fofc_unflagged') if k in u)
        return s
    with open(OUT + 'outliers_%s.txt' % star, 'w') as f:
        f.write('# %s energy outliers: |R - med| > 10 MAD (MAD %.3e of Etot); columns at log lags -1, 0, +1 of the cycle (cycle index of the log line may be off by one)\n' % (star, mad))
        for i in out[np.argsort(-abs(allx[out]))][:400]: f.write(row(i, 'dE/E') + '\n')
        f.write('# mass outliers (MAD %.3e of M_tot)\n' % madm)
        for i in om[np.argsort(-abs(allm[om]))][:200]: f.write(row(i, 'dM/M') + '\n')
    tp = out[np.argsort(-abs(allx[out]))][:8]
    for i in tp: pr('   E-outlier', row(i, 'dE/E')[:430])
    for i in om[np.argsort(-abs(allm[om]))][:5]: pr('   M-outlier', row(i, 'dM/M')[:330])
    # per-chunk switch rows: first row step of every later chunk relative to the last row of the previous chunk
    for k in range(1, len(ch)):
        n0, a0 = ch[k - 1]; n1, a1 = ch[k]; e = lambda nn, aa, nm, i: aa[i, nn.index(nm)]
        rr = lambda nn, aa, i: aa[i, nn.index('L_in')] - aa[i, nn.index('L_top')] + aa[i, nn.index('P_esrc')] - aa[i, nn.index('P_sponge')] - aa[i, nn.index('E_hflx')]
        dts = a1[0, 0] - a0[-1, 0]; sw_e = (e(n1, a1, 'Etot', 0) - e(n0, a0, 'Etot', -1) - dts * rr(n1, a1, 0)) / e(n1, a1, 'Etot', 0); sw_m = (e(n1, a1, 'M_tot', 0) - e(n0, a0, 'M_tot', -1) + dts * a1[0, n1.index('M_hflx')]) / e(n1, a1, 'M_tot', 0)
        pr('  SWITCH residual (dE - dt rate of first row of the new chunk)/Etot %+.3e, (dM + dt M_hflx)/M_tot %+.3e, dt %.4e (the rows straddle the two binaries; first new row vs last old row)' % (sw_e, sw_m, dts))
        pr('  SWITCH chunk %d -> %d: last row t %.10e (Etot %.16e, M_tot %.16e), first row t %.10e: dt %.4e, dM/M %.3e, dE/E %.3e' % (k - 1, k, a0[-1, 0], e(n0, a0, 'Etot', -1), e(n0, a0, 'M_tot', -1), a1[0, 0], a1[0, 0] - a0[-1, 0], e(n1, a1, 'M_tot', 0) / e(n0, a0, 'M_tot', -1) - 1, e(n1, a1, 'Etot', 0) / e(n0, a0, 'Etot', -1) - 1))
    # FOFC / floor counts per cycle
    for k in ('fofc_retested', 'fofc_unflagged', 'fofc_rep_other', 'fofc_rep_at', 'c2p_dfloor', 'floored_update', 'skip_update', 'floor_pay', 'paid', 'fb_pos_gas', 'fcap_v', 'auto_be'):
        v = np.array([x[k] for x in C.values() if k in x])
        if len(v): pr('  counter %-15s cycles %d: median %.4g, p99 %.4g, max %.4g (nonzero %d)' % (k, len(v), np.median(v), np.percentile(v, 99), v.max(), (v != 0).sum()))
open(OUT + 'summary.txt', 'w').write('\n'.join(P) + '\n')
