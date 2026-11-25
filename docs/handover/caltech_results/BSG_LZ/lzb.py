#!/usr/bin/env python3
"""zzzze9 (READ ONLY): L_z residual r[c] = Lz_tot(c) - Lz_tot(c-1) + 0.5 (hflx(c-1) + hflx(c)) dt against the printed counters, per cycle, and
against the ledger sampled in the bin headers (parameter dump), per dump window.  Row index = cycle (fresh per-cycle runs)."""
import re, sys, glob, os, numpy as np
J = '/resnick/groups/carnegie_poc/jingze/'
RUN = [('c1', 'prod_c1_1124/bsg', 'bsg3d', ['log.0'], [(6050, 7000), (8500, 9000)]), ('T2', 'prod_rc4v_1124/bsg', 'bsg3d', ['log.0', 'log.1'], [(6061, 6900)])]
def rhst(f):
    n = None; rows = []
    for ln in open(f, errors='replace'):
        if ln.startswith('#'):
            if '[1]=' in ln: n = [m.split('=', 1)[1] for m in re.findall(r'\[\d+\]=\S+', ln)]
        elif n and len(ln.split()) == len(n): rows.append([float(x) for x in ln.split()])
    return n, np.array(rows)
def spear(x, y):
    rk = lambda q: np.argsort(np.argsort(q)).astype(float)
    return np.corrcoef(rk(x), rk(y))[0, 1] if np.std(x) > 0 and np.std(y) > 0 else np.nan
for lab, d, base, logs, wins in RUN:
    n, a = rhst(J + d + '/' + base + '.user.hst'); t = a[:, 0]; lz = a[:, n.index('Lz_tot')]; lf = a[:, n.index('Lz_hflx')]
    R = np.zeros(len(t)); R[1:] = np.diff(lz) + 0.5 * (lf[1:] + lf[:-1]) * np.diff(t)
    C = {}                                                              # cycle -> dict of columns
    def put(c, k, v): C.setdefault(c, {})[k] = v
    prev = {}
    for lg in logs:
        prev = {}
        for ln in open(J + d + '/' + lg, errors='replace'):
            if ln.startswith('CLEAN FALLBACK cycle'):
                m = re.match(r'CLEAN FALLBACK cycle (\d+): pos_gas (\S+) .*? e_floor (\S+) .*?ffct_left (\S+) auto_be (\S+) floorcap v (\S+) T (\S+) c2p_dfloor (\S+) \(mass (\S+)\) floor_pay (\S+) \(unpaid (\S+)\)', ln)
                if m:
                    c = int(m.group(1)); [put(c, k, float(v.rstrip(','))) for k, v in zip(('fb_pos_gas', 'fb_e_floor', 'fb_ffct_left', 'auto_be', 'fcap_v', 'fcap_T', 'c2p_dfloor', 'c2p_dfloor_mass', 'floor_pay', 'floor_unpaid'), m.groups()[1:])]
            elif 'floor_audit cycle' in ln:
                m = re.search(r'floor_audit cycle (\d+) \[(update|refresh)\] cell-c2p (\S+) .*?floored (\S+) .*?skip\(dfloor/vceil\) (\S+) skip_de (\S+) skip_dm (\S+)', ln)
                if m:
                    c = int(m.group(1)); cur = [float(x) for x in m.groups()[2:]]; p = prev.get(m.group(2))
                    if p and cur[0] >= p[0]:
                        for k, (u, v) in zip(('floored', 'skip', 'skip_de', 'skip_dm'), [(cur[i], p[i]) for i in (1, 2, 3, 4)]): put(c, '%s_%s' % (k, m.group(2)), u - v)
                    prev[m.group(2)] = cur
                m = re.search(r'floor_audit cycle (\d+) pay: paid (\S+) unpaid (\S+) given_to_rad (\S+) cells (\S+) \| fofc_retest: retested (\S+) unflagged (\S+)', ln)
                if m:
                    c = int(m.group(1)); cur = [float(x) for x in m.groups()[1:]]; p = prev.get('pay')
                    if p and cur[4] >= p[4]:
                        for k, i in (('pay_paid', 0), ('pay_unpaid', 1), ('given_to_rad', 2), ('pay_cells', 3), ('fofc_retested', 4), ('fofc_unflagged', 5)): put(c, k, cur[i] - p[i])
                    prev['pay'] = cur
            elif ln.startswith('CLEAN cycle'):
                m = re.match(r'CLEAN cycle (\d+) t2s (\d) it (\d+) \(sum (\d+)\) pos (\d+) mark (\d+) all (\d+) neg (\d+) r1 (\d+) r2 (\d+) .*?guard (\d+) solves (\d+)', ln)
                if m: 
                    c = int(m.group(1)); s = m.group(2)
                    for k, v in zip(('it', 'sum', 'pos', 'mark', 'all', 'neg', 'r1', 'r2', 'guard', 'solves'), m.groups()[2:]): put(c, 't2s%s_%s' % (s, k), float(v))
            elif 'floorcap cycle' in ln and 'dL_z' in ln:
                m = re.search(r'floorcap cycle (\d+) .*?dL_z (\S+)', ln); 
                if m: put(int(m.group(1)), 'fcap_dLz_cum', float(m.group(2)))
    keys = sorted({k for v in C.values() for k in v} - {'fcap_dLz_cum'})
    print('\n# %s: hst rows %d (cycle = row), columns parsed per cycle: %d keys: %s' % (lab, len(t), len(keys), ' '.join(keys)))
    for lo, hi in wins:
        cc = np.array([c for c in range(lo, min(hi, len(t) - 2)) if c in C]); 
        print('\n## %s c%d-%d: %d cycles with log lines; residual sum %.4e, rms %.4e, max |r| %.4e' % (lab, lo, hi, len(cc), R[cc].sum(), np.sqrt(np.mean(R[cc] ** 2)), abs(R[cc]).max()))
        print('%-18s %9s %9s %9s %9s  (Pearson with r, with |r|; Spearman with r, with |r|; best of lags -1,0,+1 by |Pearson with |r||, lag; nonzero cycles)' % ('counter', 'P(r)', 'P(|r|)', 'S(r)', 'S(|r|)'))
        out = []
        for k in keys:
            best = None
            for lag in (-1, 0, 1):
                q = np.array([C.get(c + lag, {}).get(k, np.nan) for c in cc]); ok = np.isfinite(q)
                if ok.sum() < 50 or np.std(q[ok]) == 0: continue
                x = q[ok]; y = R[cc][ok]; v = (np.corrcoef(x, y)[0, 1], np.corrcoef(x, abs(y))[0, 1], spear(x, y), spear(x, abs(y)))
                if best is None or abs(v[1]) > abs(best[0][1]): best = (v, lag, int((x != 0).sum()), ok.sum())
            if best: out.append((abs(best[0][1]), k, best))
            else:
                q = np.array([C.get(c, {}).get(k, np.nan) for c in cc]); out.append((-1, k, (None, 0, int(np.nansum(q != 0)), int(np.isfinite(q).sum()))))
        for s, k, b in sorted(out, reverse=True):
            if b[0] is None: print('%-18s constant or < 50 values (nonzero %d of %d)' % (k, b[2], b[3]))
            else: print('%-18s %+9.3f %+9.3f %+9.3f %+9.3f  lag %+d, nonzero %d of %d' % (k, *b[0], b[1], b[2], b[3]))
        top = cc[np.argsort(-abs(R[cc]))][:10]
        show = [k for k in ('auto_be', 'fb_pos_gas', 'fb_e_floor', 'fb_ffct_left', 'fcap_v', 'fcap_T', 'c2p_dfloor', 'floor_pay', 'floored_update', 'skip_update', 'skip_refresh', 'given_to_rad', 'pay_cells', 'fofc_retested', 'fofc_unflagged', 't2s2_it', 't2s3_it', 't2s2_neg', 't2s3_neg', 't2s3_guard', 't2s3_mark') if k in keys]
        print('\n10 largest |residual| cycles: cycle t residual dt dLz_tot | hst: Picard | ' + ' '.join(show))
        iP = n.index('Picard')
        for c in top:
            print('c%-6d t %.5e r %+.4e dt %.3f dLz %+.4e | Picard %.0f | ' % (c, t[c], R[c], t[c] - t[c - 1], lz[c] - lz[c - 1], a[c, iP]) + ' '.join('%s=%g' % (k, C.get(c, {}).get(k, np.nan)) for k in show))
        fl = sorted((c, v['fcap_dLz_cum']) for c, v in C.items() if 'fcap_dLz_cum' in v and lo - 30 <= c <= hi + 30)
        print('floorcap cumulative dL_z in the window (25-cycle lines): first c%s %.4e, last c%s %.4e' % (fl[0][0], fl[0][1], fl[-1][0], fl[-1][1]) if fl else 'no floorcap dL_z lines in the window')
    # ledger from the bin headers
    L = {}
    for f in sorted(glob.glob(J + d + '/bin/*.bin')):
        try:
            with open(f, 'rb') as fh: h = fh.read(140000)
        except OSError: continue
        mc = re.search(rb'cycle=(\d+)', h[:400]); ml = re.search(rb'clean_lz_ledger\s*=\s*([-0-9.eE+]+)', h)
        if mc and ml: L[int(mc.group(1))] = float(ml.group(1))
    cs = sorted(L); print('\n## %s ledger from %d bin headers (cycles %d-%d); window = between consecutive dump cycles' % (lab, len(cs), cs[0], cs[-1]))
    print('cycles           d(ledger)     sum residual   dLz_tot       int hflx dt    -d(ledger)/residual')
    X, Y = [], []
    for u, v in zip(cs[:-1], cs[1:]):
        if v < len(t):
            r = R[u + 1:v + 1].sum(); dl = L[v] - L[u]; X.append(dl); Y.append(r)
            print('%6d-%-6d  %+.5e  %+.5e  %+.5e  %+.4e   %+.3e' % (u, v, dl, r, lz[v] - lz[u], (0.5 * (lf[u + 1:v + 1] + lf[u:v]) * np.diff(t)[u:v]).sum(), -dl / r if r else np.nan))
    X, Y = np.array(X), np.array(Y); cm = np.array(cs[1:len(X) + 1])
    for nm, q in (('all windows', np.ones(len(X), bool)), ('windows ending after c6050', cm > 6050), ('before c6050', cm <= 6050)):
        if q.sum() > 3: print('%s: N %d  Pearson(d ledger, residual) %+.4f, Spearman %+.4f; |d ledger| median %.2e, |residual| median %.2e' % (nm, q.sum(), np.corrcoef(X[q], Y[q])[0, 1], spear(X[q], Y[q]), np.median(abs(X[q])), np.median(abs(Y[q]))))
