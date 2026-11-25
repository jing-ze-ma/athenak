#!/usr/bin/env python3
"""zzzzc5: t9l62 histories (READ ONLY).  t9hist.py OUTDIR : runs listed below; user.hst columns by header name; rows a restart rewrote dropped."""
import os, re, sys, glob, numpy as np
J = '/resnick/groups/carnegie_poc/jingze/'
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
QUIET_DT = 96.85
RUNS = """prod_fresh_1121/t9l62 prod_fresh_1121b/t9l62 prod_ia9_1123/t9l62 prod_ia11_1123/t9l62 prod_ia121_1123/t9l62 prod_rc3_1124/t9l62 prod_rc4v_1124/t9l62
prod_c4_1124/t9l62 sdo_1124/S5off sdo_1124/S5on sdo_1124/S5Foff sdo_1124/S5Fon sdo_1124/S5Goff sdo_1124/S5Gon sdo_1124/S4off sdo_1124/S4on t9hot2_1124/B t9hot2_1124/F
t9hot2_1124/K t9hot2_1124/P1 t9arms_1123/L t9arms_1123/S t9q_1123/LAY0 t9q_1123/LAYF8 t9q_1123/T0 t9q_1123/T0a t9q_1123/T0d t9q_1123/T0da t9q_1123/T0h t9q_1123/W
t9rays_1124/G1 t9rays_1124/G2b t9rays_1124/pre_G2b t9wk_1123/SB t9wk_1123/WA t9wk_1123/WB thinf_1124/C0t thinf_1124/T1 thinf_1124/T2 tests_1124/WD0 tests_1124/WD2 unif_1124/PDtB
unif_1124/PStA unif_1124/PStB unif_1124/U4A unif_1124/U4B unif_1124/U4C cand_1124/iA cand_1124/iB cand_1124/iC cand_1124/iD cand_1124/iE cand_1124/PL cand_1124/PO cand_1124/PQ
cand_1124/pre_t9l62 cand_1124/t9l62_F cand_1124/t9l62_R compton_1124/C0 compton_1124/C0b compton_1124/C1 compton_1124/C1ra compton_1124/C1rb compton_1124/C2 relsmoke_1124/t9l62_B
relsmoke_1124/t9l62_R ia9_t9diag_1123/A proj_1124/t9o4 t9l62_1113/prod t9l62_1113/s10 t9l62_1113/s20""".split()
def rd(f):
    names = None; rows = []
    for ln in open(f, errors='replace'):
        if ln.startswith('#'):
            if '[1]=' in ln:
                names = [m.split('=', 1)[1] for m in re.findall(r'\[\d+\]=\S+', ln)]
        elif ln.split() and names and len(ln.split()) == len(names):
            rows.append([float(x) for x in ln.split()])
    a = np.array(rows) if rows else np.zeros((0, len(names or [1])))
    nraw = len(a)
    if len(a):
        t = a[:, 0]; keep = t < np.minimum.accumulate(np.append(t[1:], np.inf)[::-1])[::-1]
        a = a[keep]
    return names, a, nraw
P = []; p = lambda *x: P.append(' '.join(str(y) for y in x))
p('run | user.hst rows (raw -> kept) | columns | t range [s]')
INFO = {}
for r in RUNS:
    f = J + r + '/t9l62.user.hst'
    if not os.path.exists(f):
        p('%s | MISSING user.hst' % r); continue
    n, a, nr = rd(f); INFO[r] = (n, a)
    c = J + r + '/chunks.txt'; rl = J + r + '/run.log'
    prov = ''
    if os.path.exists(c):
        L = open(c, errors='replace').read(); prov = ' ; '.join(sorted(set(re.findall(r'BIN=\S*/([^/\s]+)', L)))[:8]) + ' ; first start ' + str(re.findall(r'start=(-[ir])', L)[:1])
    elif os.path.exists(rl):
        prov = 'run.log: ' + ' / '.join(open(rl, errors='replace').read().splitlines()[:2])[:200]
    p('%s | %d -> %d | %d cols | %.4e .. %.4e | %s' % (r, nr, len(a), len(n or []), a[0, 0] if len(a) else 0, a[-1, 0] if len(a) else 0, prov))
open(OUT + '/runs.txt', 'w').write('\n'.join(P) + '\n')
# ------------------------------------------------ per-run series on a 2e4 s grid (runs spanning >= 1e5 s or > 1e3 rows)
LG = ['']
def col(n, a, k):
    return a[:, n.index(k)] if k in n else None
for r, (n, a) in INFO.items():
    if len(a) < 20 or a[-1, 0] - a[0, 0] < 3e4:
        continue
    t = a[:, 0]; Lt, Lb, Li = col(n, a, 'L_top'), col(n, a, 'L_bot'), col(n, a, 'L_in')
    KE, Kr, M, Mt, dt = col(n, a, 'KE_int'), col(n, a, 'KEr_int'), col(n, a, 'M_tot'), col(n, a, 'Mdot_top'), col(n, a, 'dt')
    LG.append('\n## %s  (L_top/L_in, L_bot/L_in, KE_int[erg], KEr_int, M_tot-M_tot(first)[g], Mdot_top[g/s, mean of rows], cumulative mean Mdot_top x dt[g], dt median [s]; row windows of 2e4 s: nearest-row averages)' % r)
    g0 = np.arange(np.floor(t[0] / 2e4) * 2e4, t[-1] + 2e4, 2e4); cum = 0
    M0 = M[0] if M is not None else np.nan
    for lo in g0:
        s = (t >= lo) & (t < lo + 2e4)
        if not s.any():
            continue
        v = lambda x: ('%.4g' % np.mean(x[s])) if x is not None else 'n/a'
        dtm = np.median(dt[s]) if dt is not None else np.nan
        if Mt is not None:
            cum += np.mean(Mt[s]) * (min(lo + 2e4, t[-1]) - max(lo, t[0]))
        mxdt = np.min(dt[s]) if dt is not None else np.nan
        LG.append('%9.3e %9.4f %9.4f %10.3e %10.3e %12.4e %11.3e %11.3e %9.2f (min %.2f)' % (lo, np.mean(Lt[s] / Li[s]) if Lt is not None else np.nan, np.mean(Lb[s] / Li[s]) if Lb is not None else np.nan, np.mean(KE[s]) if KE is not None else np.nan, np.mean(Kr[s]) if Kr is not None else np.nan, np.mean(M[s]) - M0 if M is not None else np.nan, np.mean(Mt[s]) if Mt is not None else np.nan, cum, dtm, mxdt))
open(OUT + '/series_2e4.txt', 'w').write('\n'.join(LG) + '\n')
# ------------------------------------------------ episodes: orion's dt definition
E = ['episodes: dt < 0.5 x %.2f s for > 10 consecutive cycles (orion zzzzc5, arbitrary). Rows where the run is in an episode; runs whose hst cadence is not 1 per cycle are marked (the cycle count cannot be used; see dt<0.5 share).' % QUIET_DT]
for r, (n, a) in INFO.items():
    if len(a) < 30 or 'dt' not in n:
        continue
    t = a[:, 0]; dt = col(n, a, 'dt'); ddt = np.median(np.diff(t)) / np.median(dt)
    low = dt < 0.5 * QUIET_DT
    LT, KE, Li = col(n, a, 'L_top'), col(n, a, 'KE_int'), col(n, a, 'L_in')
    ep = []; i = 0
    ext = lambda i, j: (t[min(j + 1, len(t) - 1)] - t[i]) if j + 1 < len(t) else (t[j] - t[i] + dt[j])
    while i < len(t):
        if low[i]:
            j = i
            while j + 1 < len(t) and low[j + 1]:
                j += 1
            ncyc = (ext(i, j) / np.mean(dt[i:j + 1])) if ddt > 1.5 else (j - i + 1)      # sampled hst: cycles estimated from elapsed time / mean dt
            if ncyc > 10:
                ep.append((i, j))
            i = j + 1
        else:
            i += 1
    E.append('\n## %s%s  rows %d, cadence %.2f rows/cycle-equivalent (median dt(t)/dt), share of rows with dt < %.1f: %.3f; n episodes %d' % (r, '' if ddt < 1.5 else '  [SAMPLED hst: one row per ~%.0f cycles; cycles estimated]' % ddt, len(t), 1 / ddt if ddt > 0 else np.nan, 0.5 * QUIET_DT, low.mean(), len(ep)))
    for i, j in ep:
        s = slice(i, j + 1)
        E.append('   t %.4e - %.4e s (rows %d), min dt %.2f, peak L_top/L_in %.4f, peak KE_int %.3e' % (t[i], t[j], j - i + 1, dt[s].min(), (LT[s] / Li[s]).max(), KE[s].max()))
    if not ep:
        E.append('   none')
    E.append('   whole run: min dt %.2f at t %.4e; max L_top/L_in %.4f at t %.4e; max KE_int %.3e at t %.4e; KE at first row %.3e' % (dt.min(), t[dt.argmin()], (LT / Li).max(), t[(LT / Li).argmax()], KE.max(), t[KE.argmax()], KE[0]))
open(OUT + '/episodes.txt', 'w').write('\n'.join(E) + '\n')
print('done', len(INFO))
