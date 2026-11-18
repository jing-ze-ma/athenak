#!/usr/bin/env python3
"""print '<hydro_w> <m1|-> <tag>' for every dump not yet reduced (or reduced without its m1 partner that has appeared since).
Static phase = STATIC bin dir (t <= 0.01, tag a_NNNNN), stream = RUN/stream/bin (non-symlink, tag b_NNNNN). m1 matched by header
time; files younger than AGE s are skipped (may still be written)."""
import os, re, glob, sys, time
STATIC, STREAM, DATA = sys.argv[1], sys.argv[2], sys.argv[3]
AGE = float(os.environ.get('AGE', '180'))
def ht(f):
    with open(f, 'rb') as fh:
        m = re.search(rb'time=\s*([0-9.eE+-]+)', fh.read(4096))
    return float(m.group(1)) if m else None
H, M = [], {}
for tag, d in (('a', STATIC), ('b', STREAM)):
    for f in sorted(glob.glob(f'{d}/ryper.hydro_w.*.bin')):
        if not os.path.islink(f): H.append((f'{tag}_' + f.split('.')[-2], f))
    for f in sorted(glob.glob(f'{d}/ryper.m1.*.bin')):
        if not os.path.islink(f): M[f] = ht(f)
now = time.time()
for tag, f in H:
    if now - os.path.getmtime(f) < AGE: continue
    t = ht(f)
    if tag.startswith('a') and (t is None or t > 0.01 + 1e-9): continue
    o = f'{DATA}/red_{tag}.npz'
    m = [x for x, tm in M.items() if tm is not None and abs(tm - t) < 1e-6 and now - os.path.getmtime(x) > AGE]
    if os.path.exists(o):
        import numpy as np
        try:
            z = np.load(o)
            if int(z['ver']) >= 2 and (bool(z['has_m1']) or not m): continue   # trc-1102: VER 2 adds the tracer reductions
        except Exception:
            pass
    print(f, m[0] if m else '-', tag)
