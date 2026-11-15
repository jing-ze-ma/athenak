"""idort_page Q2: bubbles metrics against the retarded (light-cone) and steady free-streaming references (bubref.py).
L1(E) = sum|E - E_ref| / sum E_ref ;  E_tot = sum E / sum E_ref ;  L1(F) = sum|F - F_ref| / sum|F_ref| (vector F, |F| <= cE)
on the reference sample (512^2 cell centres every <every> cells; outside the discs and inside, as the suite).
Runs: suiteF_1012 runs_fin (FIN) and runs_d (PRODd), periodic unit box; plus the Q2 2x2-box rerun (central unit box).
usage: python bubeval.py [every]   -> cache/bubeval_<every>.json and prints a table"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, '/raven/u/jinma/ATHENAK/wt-prodfinal/vis/python')
import bin_convert as bc  # noqa: E402

CACHE = '/raven/ptmp/jinma/idort_page/cache'
F = '/raven/ptmp/jinma/suiteF_1012'
RUNS = [('FIN_r100k5', F + '/runs_fin/bub/FIN_r100k5', '2', 1),
        ('FIN_r100k3', F + '/runs_fin/bub/FIN_r100k3', '2', 1),
        ('FIN_g8x32k3', F + '/runs_fin/bub/FIN_g8x32k3', '3', 1),
        ('PRODd_r100k5', F + '/runs_d/bub/PRODd_r100k5', '2', 1),
        ('PRODd_r100k3', F + '/runs_d/bub/PRODd_r100k3', '2', 1),
        ('PRODd_g8x32k3', F + '/runs_d/bub/PRODd_g8x32k3', '3', 1),
        ('FIN_r100k5_box2', '/raven/ptmp/jinma/idort_page/reruns/bub_FIN_r100k5_box2', '2', 2)]
TS = (0.005, 0.025, 0.1)


def load_at(d, t):
    best = None
    for f in sorted(glob.glob(os.path.join(d, 'bin', '*.m1.*.bin'))):
        r = bc.read_binary_as_athdf(f)
        if best is None or abs(r['Time'] - t) < abs(best['Time'] - t):
            best = r
    return best


def fields(r, box):
    E, F1, F2 = r['m1_e'][0], r['m1_f1'][0], r['m1_f2'][0]
    if box == 2:     # central unit box of the 1024^2 [-1,1]^2 run (same dx)
        s = slice(256, 768)
        E, F1, F2 = E[s, s], F1[s, s], F2[s, s]
    return E, F1, F2


def metrics(E, F1, F2, Er, Fxr, Fyr):
    return dict(L1=float(np.abs(E - Er).sum()/Er.sum()), Etot=float(E.sum()/Er.sum()),
                L1F=float((np.hypot(F1 - Fxr, F2 - Fyr)).sum()/np.hypot(Fxr, Fyr).sum()))


def main(every=4):
    R = np.load(os.path.join(CACHE, 'bubref_%d.npz' % every))
    out = {}
    for name, d, dim, box in RUNS:
        if not glob.glob(os.path.join(d, 'bin', '*.m1.*.bin')):
            continue
        o = {}
        for k, t in enumerate(TS):
            r = load_at(d, t)
            if abs(r['Time'] - t) > 1e-6:
                continue
            E, F1, F2 = fields(r, box)
            sl = (slice(every//2, None, every),)*2
            E, F1, F2 = E[sl], F1[sl], F2[sl]
            row = {'time': float(r['Time'])}
            refs = [('per_ret', 'per', k), ('per_steady', 'per', 3)]
            if box == 2:
                refs += [('open_ret', 'open', k)]
            for lab, bcn, kk in refs:
                row[lab] = metrics(E, F1, F2, R['E%s_%s' % (dim, bcn)][kk], R['F%sx_%s' % (dim, bcn)][kk],
                                   R['F%sy_%s' % (dim, bcn)][kk])
            o['%.3f' % t] = row
        out[name] = dict(ref_dim=dim, box=box, t=o)
    fn = os.path.join(CACHE, 'bubeval_%d.json' % every)
    json.dump(out, open(fn, 'w'), indent=1)
    for name, v in out.items():
        for t, row in v['t'].items():
            print('%-16s t=%s  ' % (name, t) + '  '.join(
                '%s L1 %.3f Etot %.3f L1F %.3f' % (lab, m['L1'], m['Etot'], m['L1F'])
                for lab, m in row.items() if lab != 'time'))
    print('wrote', fn)


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
