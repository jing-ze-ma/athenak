"""cr2ana.py (CR2 caltech_rol2u_1015; adapted from Raven rol2u_drift_1015/ana/drf.py): ROL2U vs BASE window comparison.
usage: python3 -I cr2ana.py --mode e1|hg [--prefix P] [--rph CM] [--L0 L] [--t0 T] [--tmax T] [--bc DIR_OF_bin_convert.py]
                           BASE_DIR NAME:DIR [NAME:DIR ...]
  NAME conventions: R = ROL2U (chain dir; only its part inside the window is used), H = BASE cfl/2 (optional yardstick),
                    P = production (FIN4) dir (optional; read-only use of its hst/bins/logs).
Window: [max of the arms' first hst time, min of their last hst time] (or --t0/--tmax); arms interpolated onto BASE hst times.
hst (user.hst, dedup over links, last occurrence wins): L_top/L0, L_bot/L0, KE_int, Picard (column 'Picard' if present).
bins (bin/*.hydro_w.*.bin at times common to all arms, inside the window): sin-weighted shell means as pulse_ana.py:
  e1: rho and v_r at 1.05 R_ph (lrho105 = log10 <rho>, vr105 km/s), vr110, rsurf (outermost <rho> > 1e-12, in R_ph), KEe (0.9-1.3 R_ph)
  hg: r12 (outermost <rho> > 1e-12, Rsun), lrho_ph / vr_ph at --rph (cm; the R_ph of the guard at the window start),
      lrho_12 / vr_12 at the r12 of that dump, vmax (max shell-mean v_r, km/s)
logs (run.log.* or run*.log of the dir, only cycles inside the window for s/cycle): s/cycle (median d elapsed / d cycle over ndiag lines),
  DT_RETRY lines, 'Picard NON-CONVERGED after' lines + sum of the NON-CONVERGED= report values, backward-Euler steps (sum of the report
  lines), floor clips (solved E <= e_floor cell-solves; written-back eint <= 0), lagged solves, FATAL.  For a long chain (R) the report
  lines cover whole links: pass --rlog <run.log of the window link only> semantics by giving a dir that holds only that link's log.
FLAG: an observable whose rms(R - BASE) over the window exceeds rms(H - BASE) (when H is given) is printed with FLAG.
"""
import argparse, glob, os, re, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--mode', choices=('e1', 'hg'), required=True)
ap.add_argument('--prefix', default=None); ap.add_argument('--rph', type=float, default=None)
ap.add_argument('--L0', type=float, default=None); ap.add_argument('--t0', type=float, default=None)
ap.add_argument('--tmax', type=float, default=None); ap.add_argument('--nwin', type=int, default=6)
ap.add_argument('--bc', default=os.path.dirname(os.path.abspath(__file__))); ap.add_argument('--nobin', action='store_true')
ap.add_argument('base'); ap.add_argument('arms', nargs='*')
a = ap.parse_args()
PFX = a.prefix or ('agcar3d' if a.mode == 'e1' else 'hegiant')
RSUN = 6.957e10
RPH = a.rph or (2.70140777e13 if a.mode == 'e1' else None)
sys.path.insert(0, a.bc)


def hst(d):
    f = os.path.join(d, PFX + '.user.hst')
    hdr = [l for l in open(f) if l.startswith('#  [')]
    names = re.findall(r'\[\d+\]=(\S+)', hdr[-1]) if hdr else []
    rows = [l.split() for l in open(f) if not l.startswith('#')]
    x = np.array([[float(v) for v in r] for r in rows if len(r) == len(names)])
    keep = []; tmax = np.inf
    for i in range(len(x) - 1, -1, -1):
        if x[i, 0] < tmax: keep.append(i); tmax = x[i, 0]
    x = x[keep[::-1]]
    return {n: x[:, k] for k, n in enumerate(names)}


arms = [('BASE', a.base)] + [tuple(s.split(':', 1)) for s in a.arms]
H = {n: hst(d) for n, d in arms}
t0 = max(H[n]['time'][0] for n, _ in arms) if a.t0 is None else a.t0
t1 = min(H[n]['time'][-1] for n, _ in arms)
if a.tmax: t1 = min(t1, a.tmax)
tb = H['BASE']['time']; tt = tb[(tb >= t0) & (tb <= t1)]
L0 = a.L0 or H['BASE']['L_top'][np.searchsorted(tb, t0)]
print(f'# mode {a.mode} window t0 {t0:.6e} -> t1 {t1:.6e} (dt {t1 - t0:.3e} s), {len(tt)} BASE hst rows; L0 {L0:.4e}'
      + ('' if a.L0 else ' (= BASE L_top at t0)'))
for n, d in arms: print(f'#   {n:5s} {d}')
Q = {}
for n, _ in arms:
    h = H[n]; ti = h['time']; q = {}
    for k, c, s in (('Ltop', 'L_top', L0), ('Lbot', 'L_bot', L0), ('KE_int', 'KE_int', 1.0), ('Picard', 'Picard', 1.0)):
        if c in h: q[k] = np.interp(tt, ti, h[c]/s)
    Q[n] = q


def logs(d):
    fs = sorted(set(glob.glob(os.path.join(d, 'run.log.*')) + glob.glob(os.path.join(d, 'run*.log'))))
    fs = [f for f in fs if not os.path.islink(f)]
    o = dict(retry=0, ncl=0, ncrep=0.0, be=0.0, fe=0.0, fg=0.0, lag=0.0, fatal=0, files=len(fs)); dts = []
    for f in fs:
        last = None
        for ln in open(f, errors='replace'):
            if ln.startswith('DT_RETRY cycle='): o['retry'] += 1
            elif 'Picard NON-CONVERGED after' in ln: o['ncl'] += 1
            if 'FATAL' in ln: o['fatal'] += 1
            m = re.search(r'elapsed=([0-9.eE+-]+) cycle=(\d+) time=([0-9.eE+-]+)', ln)
            if m:
                e, c, t = float(m.group(1)), int(m.group(2)), float(m.group(3))
                if last and t0 <= t <= t1 and c > last[1] and e > last[0]: dts.append((e - last[0])/(c - last[1]))
                last = (e, c)
            m = re.search(r'NON-CONVERGED=([0-9.eE+-]+)', ln)
            if m: o['ncrep'] += float(m.group(1))
            m = re.search(r'backward-Euler steps=([0-9.eE+-]+)', ln)
            if m: o['be'] += float(m.group(1))
            m = re.search(r'e_floor cell-solves=([0-9.eE+-]+) written-back gas eint <= 0 cell-solves=([0-9.eE+-]+)', ln)
            if m: o['fe'] += float(m.group(1)); o['fg'] += float(m.group(2))
            m = re.search(r'lagged solves=([0-9.eE+-]+)', ln)
            if m: o['lag'] += float(m.group(1))
    o['spc'] = float(np.median(dts)) if dts else np.nan
    return o


print('# logs (window cycles for s/cycle; counters = all logs in the dir: keep window arms in their own dirs)')
print('%-5s %8s %7s %6s %8s %8s %9s %10s %10s %10s %5s' % ('arm', 's/cyc', 'Picard', 'retry', 'NClines', 'NCrep', 'BEsteps',
                                                          'EfloorClip', 'eint<=0', 'lagged', 'FATAL'))
for n, d in arms:
    o = logs(d); pic = Q[n]['Picard'].mean() if 'Picard' in Q[n] else np.nan
    print('%-5s %8.4f %7.3f %6d %8d %8.0f %9.0f %10.0f %10.0f %10.0f %5d' % (n, o['spc'], pic, o['retry'], o['ncl'], o['ncrep'],
                                                                             o['be'], o['fe'], o['fg'], o['lag'], o['fatal']))

ed = np.linspace(t0, t1, a.nwin + 1)
keys = [k for k in ('Ltop', 'Lbot', 'KE_int', 'Picard') if all(k in Q[n] for n, _ in arms)]
print('# hst window means (Ltop/Lbot in L0)')
print('%-5s %-11s' % ('arm', 'win_end') + ''.join('%12s' % k for k in keys))
for n, _ in arms:
    for k in range(a.nwin):
        s = (tt >= ed[k]) & (tt <= ed[k + 1])
        if s.any(): print('%-5s %-11.5e' % (n, ed[k + 1]) + ''.join('%12.6g' % Q[n][q][s].mean() for q in keys))


def rmsd(n, q, rel=False):
    dd = Q[n][q] - Q['BASE'][q]
    if rel: dd = dd/np.mean(np.abs(Q['BASE'][q]))
    return np.sqrt(np.mean(dd**2)), (np.polyfit((tt - t0)/1e5, dd, 1)[0] if len(tt) > 3 else np.nan)


print('# hst: rms(arm - BASE) over the window [slope per 1e5 s]; KE_int relative; FLAG = R rms > H rms')
R = {}
for n, _ in arms[1:]:
    line = f'{n:5s}'
    for q in keys:
        r, sl = rmsd(n, q, rel=(q == 'KE_int')); R[(n, q)] = r
        line += f' {q} {r:.3g} [{sl:+.2g}]'
    print(line)
if 'H' in dict(arms) and 'R' in dict(arms):
    for q in keys:
        if q != 'Picard' and R[('R', q)] > R[('H', q)]: print(f'FLAG hst {q}: R-BASE rms {R[("R", q)]:.3g} > H-BASE rms {R[("H", q)]:.3g} (ratio {R[("R", q)]/R[("H", q)]:.2f})')

if a.nobin: sys.exit()
import bin_convert as bc


def shells(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, 'bin', '*.hydro_w.*.bin'))):
        h = bc.read_binary_as_athdf(f); t = float(h['Time'])
        if not (t0 - 1 <= t <= t1 + 1): continue
        th = np.asarray(h['x2v'], float); w = np.sin(th)[None, :, None]; w = w/np.sin(th).sum()/h['dens'].shape[0]
        rho = np.asarray(h['dens'], float); vr = np.asarray(h['velx'], float); r = np.asarray(h['x1v'], float)
        mr = (rho*w).sum(axis=(0, 1)); mv = (vr*w).sum(axis=(0, 1)); mrv2 = (rho*vr**2*w).sum(axis=(0, 1))
        j = np.where(mr > 1e-12)[0]; o = {'t': t}
        if a.mode == 'e1':
            x = r/RPH; wv = r**2*np.gradient(r)
            for xx, nn in ((1.05, '105'), (1.10, '110')):
                i = np.argmin(abs(x - xx)); o['vr' + nn] = mv[i]/1e5
                if nn == '105': o['lrho105'] = np.log10(mr[i])
            e = (x >= 0.9) & (x <= 1.3); o['KEe'] = (0.5*mrv2[e]*wv[e]).sum()
            o['rsurf'] = x[j.max()] if len(j) else np.nan
        else:
            i12 = j.max() if len(j) else len(r) - 1
            o['r12'] = r[i12]/RSUN; o['lrho_12'] = np.log10(mr[i12]); o['vr_12'] = mv[i12]/1e5
            if RPH:
                i = np.argmin(abs(r - RPH)); o['lrho_ph'] = np.log10(mr[i]); o['vr_ph'] = mv[i]/1e5
            o['vmax'] = mv.max()/1e5
        out[round(t, 3)] = o
    return out


S = {n: shells(d) for n, d in arms}
ts = sorted(set.intersection(*[set(S[n]) for n, _ in arms]))
sk = ('vr105', 'vr110', 'lrho105', 'rsurf', 'KEe') if a.mode == 'e1' else \
     tuple(k for k in ('r12', 'lrho_12', 'vr_12', 'lrho_ph', 'vr_ph', 'vmax') if not k.endswith('_ph') or RPH)
print('# shell metrics at common dump times (%d): t | arm values (%s)' % (len(ts), ' '.join(sk)))
for t in ts:
    print(f'{t:.6e} ' + ' '.join(f'| {n} ' + ' '.join(f'{S[n][t][k]:.5g}' for k in sk) for n, _ in arms))
print('# shells: rms(arm - BASE) over the common dumps (KEe relative); FLAG = R rms > H rms')
RS = {}
for n, _ in arms[1:]:
    r = []
    for k in sk:
        dd = np.array([S[n][t][k] - S['BASE'][t][k] for t in ts])
        if k == 'KEe': dd = dd/np.mean([S['BASE'][t][k] for t in ts])
        RS[(n, k)] = np.sqrt((dd**2).mean()) if len(dd) else np.nan; r.append(RS[(n, k)])
    print(f'{n:5s} ' + ' '.join(f'{k} {x:.3g}' for k, x in zip(sk, r)))
if 'H' in dict(arms) and 'R' in dict(arms):
    for k in sk:
        if RS[('R', k)] > RS[('H', k)]: print(f'FLAG shell {k}: R-BASE rms {RS[("R", k)]:.3g} > H-BASE rms {RS[("H", k)]:.3g} (ratio {RS[("R", k)]/RS[("H", k)]:.2f})')
