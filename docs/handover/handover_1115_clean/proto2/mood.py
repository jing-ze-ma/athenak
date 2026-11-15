"""mood.py: realizability sweep with the a posteriori per-face positivity limiter (no D floor).
Usage: python3 mood.py pos | thin <variant> <chi> | quick"""
import sys
import itertools
import numpy as np
from multiprocessing import Pool
import scheme as S
import tests as T

V = {'vl_lf': dict(cand='A3', recon='vl', mood=True),
     'vl_char': dict(cand='A3', recon='vl', mood=True, lo_mode='char'),
     'mc_char': dict(cand='A3', recon='mc', mood=True, lo_mode='char'),
     'vl_char_nb1': dict(cand='A3', recon='vl', mood=True, lo_mode='char', mood_nb=1),
     '1st_char': dict(cand='A3', recon=None, mood=True, lo_mode='char')}
if '--all' in sys.argv:      # with the last resort: whole step low order (BE, every face)
    sys.argv.remove('--all')
    V = {k + '_all': dict(v, mood_all=True) for k, v in V.items() if k != 'vl_lf'}
CASES = [(t, x, k, c) for t in ('beam', 'pulse', 'counter', 'slab')
         for x in (0.0, 0.01, 0.1, 1.0, 10.0, 100.0)
         for k in (('scat',) if x == 0 else ('scat', 'abs'))
         for c in (0.5, 4.0, 100.0, 1000.0)]


def one(arg):
    name, case = arg
    for k in S.STATS:
        S.STATS[k] = 0
    w = T.pos_case(*case, V[name])
    return name, case, w, dict(S.STATS)


def pos():
    with Pool(32) as p:
        out = p.map(one, list(itertools.product(V, CASES)), chunksize=1)
    for name in V:
        rows = [r for r in out if r[0] == name]
        bad = [r for r in rows if r[2][0] < -1e-12 or r[2][1] > 1e-12]
        st = {k: sum(r[3][k] for r in rows) for k in rows[0][3]}
        print('%-12s cases %d violations %d worst minE %.2e max(|F|-cE) %.2e | steps %d stages %d'
              ' marked faces %d (%.2f/step) re-solves %d (%.2f/stage) BE fallbacks %d (%.1f%% of'
              ' steps) all-faces steps %d BE still bad %d' % (name, len(rows), len(bad), min(r[2][0] for r in rows),
                                           max(r[2][1] for r in rows), st['steps'], st['stages'],
                                           st['marked'], st['marked']/st['steps'], st['resolves'],
                                           st['resolves']/st['stages'], st['be'],
                                           100.0*st['be']/st['steps'], st['allfaces'],
                                           st['befail']))
        for r in bad:
            print('    %-7s chi dx %-6g %-4s cfl %-6g minE %+.2e max(|F|-cE) %+.2e  marked %d'
                  ' BE %d all %d BEbad %d' % (*r[1], *r[2], r[3]['marked'], r[3]['be'],
                                              r[3]['allfaces'], r[3]['befail']))
        sys.stdout.flush()


if __name__ == '__main__':
    if sys.argv[1] == 'pos':
        pos()
    elif sys.argv[1] == 'thin':
        T.test_thin([V[sys.argv[2]]], chis=(float(sys.argv[3]),))
        print('limiter stats', S.STATS)
    elif sys.argv[1] == 'quick':
        sel = [V['vl_char'], V['mc_char']]
        for f in (T.test_thick, T.test_oddeven, T.test_cons, T.test_torder_thick,
                  T.test_torder, T.test_rwave):
            for k in S.STATS:
                S.STATS[k] = 0
            f(sel)
            print('  limiter stats after %s: %s' % (f.__name__, S.STATS))
            sys.stdout.flush()
