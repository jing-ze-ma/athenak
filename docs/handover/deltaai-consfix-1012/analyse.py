# analyse.py <runs/JOBID>: TASK-2026-10-12-deltaai-consfix-c2drag report numbers
import sys, os, re, numpy as np
S = sys.argv[1]
def hst(n):
    f = [x for x in (f'{S}/{n}/bsg3d.user.hst', f'{S}/{n}/agcar3d.user.hst') if os.path.exists(x)]
    if not f: return None
    hdr = [l for l in open(f[0]) if l.startswith('#')][-1]
    return np.loadtxt(f[0], ndmin=2), {m.group(2): int(m.group(1)) - 1 for m in re.finditer(r'\[(\d+)\]=(\S+)', hdr)}, hdr
R = {}
for n in ('bM', 'bC2', 'bC2g', 'aM', 'aC2'):
    p = f'{S}/{n}/out.txt'
    if not os.path.exists(p): continue
    o = open(p, errors='replace').read(); R[n] = o
    el = [(int(m.group(2)), float(m.group(1)), float(m.group(3))) for m in re.finditer(r'elapsed=(\S+) cycle=(\d+) time=\S+ dt=(\S+)', o)]
    d = {c1: t1 - t0 for (c0, t0, _), (c1, t1, _) in zip(el[:-1], el[1:]) if c1 == c0 + 1}
    last = el[-1][0] if el else -1
    st = [d[c] for c in d if c > 2 and c % 10 and c != last]
    na = re.findall(r'NOT ADMISSIBLE at cycle (\d+)', o)
    nn = len(re.findall(r'(?i)\bnan\b', o))
    print(f'=== {n}: FATAL {o.count("FATAL")}  NaN-lines {nn}  last cycle {last}')
    for l in o.splitlines():
        if re.search(r'implicit transport: solves|time_scheme=hesdirk2: stage steps|implicit_vimp positivity|gas coupling: newton', l): print('  ' + l.strip()[:300])
    print(f'  NOT ADMISSIBLE: {len(na)} at cycles {",".join(na)}')
    if st: print(f'  s/cycle steady (c>2, not x10, not last): mean {np.mean(st):.3f} median {np.median(st):.3f} (n {len(st)}); setup {el[0][1]:.1f} s')
    print('  dt: ' + ' '.join(f'{c}:{dt:.4g}' for c, _, dt in el))
    dbg = [m.start() for m in re.finditer('DBG_T2_ADMISS', o)]
    if dbg:
        print(f'  DBG_T2_ADMISS lines: {len(dbg)}; first 12:')
        for l in [l for l in o.splitlines() if 'DBG_T2_ADMISS' in l][:12]: print('    ' + l.strip()[:260])
    h = hst(n)
    if h: print('  L_top per row: ' + ' '.join(f'{v:.5e}' for v in h[0][:, h[1]['L_top']]))
def at_cycles(n):
    h = hst(n); t = h[0][:, h[1]['time']]
    el = {int(m.group(2)): float(m.group(1)) for m in re.finditer(r'cycle=(\d+) time=(\S+)', R[n]) for _ in [0]} if False else None
    return h
print('\n=== BSG: bC2 / bC2g vs bM at cycles 10, 20, 30 (hst row = cycle when dt is the same; else by row index)')
hb = hst('bM')
for n in ('bC2', 'bC2g'):
    h = hst(n)
    if not (hb and h): continue
    for c in (10, 20, 30):
        if c < len(hb[0]) and c < len(h[0]):
            s = ' '.join(f'{k} {abs(h[0][c, h[1][k]] - hb[0][c, hb[1][k]]) / abs(hb[0][c, hb[1][k]]):.3e}' for k in ('L_top', 'E_rad', 'Etot'))
            print(f'  {n} row {c} (t {h[0][c, h[1]["time"]]:.6e} vs {hb[0][c, hb[1]["time"]]:.6e}): {s}')
print('\n=== AG Car A: aC2 vs aM, per user.hst column, max relative difference over the rows')
a, b = hst('aM'), hst('aC2')
if a and b:
    m = min(len(a[0]), len(b[0]))
    for k, i in sorted(a[1].items(), key=lambda x: x[1]):
        x, y = a[0][:m, i], b[0][:m, i]; den = np.max(np.abs(x)) if np.max(np.abs(x)) > 0 else 1
        rel = np.max(np.abs(x - y) / np.maximum(np.abs(x), 1e-300)) if np.all(x != 0) else np.max(np.abs(x - y)) / den
        print(f'  {k:10s} {rel:.3e}{"" if np.all(x != 0) else " (normalised by max|aM|: column has zeros)"}')
