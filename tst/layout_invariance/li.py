#!/usr/bin/env python3
"""li.py: MeshBlock-layout invariance regression test.

For every case (cases.py) the same problem is run with several MeshBlock decompositions
and rank counts from identical initial states, and the GLOBAL state (double-precision .bin
output, every cell, every listed variable) is compared after cycle 1 and cycle N against
  * the base layout (bitwise reference), and
  * the base layout with a 1-ulp perturbation (default time/cfl_number*(1+2^-52)), which
    measures how far plain round-off spreads in N cycles.
Verdict per (layout, cycle, variable): BITWISE | ROUNDOFF (<= max(10*ulp_ref, 100 eps)) |
SUSPECT | FAIL (> max(1000*ulp_ref, 1e-10)).  A layout bug (the cs seam-halo clamp of
11-03) gives O(1e-4) after ONE cycle where round-off is 1e-16.

usage: li.py [--bins bins.cfg] [--work DIR] [--np-max 16] [--only case[,case]] [--list]
             [--compare-only] [--keep]
bins.cfg: lines "family = /path/to/athena" (families used in cases.py: bip, cs, sp, dhj,
he, acc).  Cases whose binary is missing are SKIPPED (reported)."""
import argparse
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import libin  # noqa: E402
import cases as CASES  # noqa: E402

EPS = np.finfo(np.float64).eps


# ---------------------------------------------------------------- athinput editing
def read_deck(fn):
    blocks, order, cur = {}, [], None
    for s in open(fn):
        t = s.split('#', 1)[0].strip()
        if not t:
            continue
        if t.startswith('<'):
            cur = t.strip('<>').strip()
            if cur not in blocks:
                blocks[cur] = []
                order.append(cur)
            continue
        if '=' in t and cur is not None:
            k, v = t.split('=', 1)
            blocks[cur].append([k.strip(), v.strip()])
    return blocks, order


def deck_set(deck, key, val):
    blocks, order = deck
    b, k = key.split('/', 1)
    if b not in blocks:
        blocks[b] = []
        order.append(b)
    for kv in blocks[b]:
        if kv[0] == k:
            kv[1] = str(val)
            return
    blocks[b].append([k, str(val)])


def deck_get(deck, key, default=None):
    b, k = key.split('/', 1)
    for kv in deck[0].get(b, []):
        if kv[0] == k:
            return kv[1]
    return default


def deck_write(deck, fn):
    blocks, order = deck
    with open(fn, 'w') as f:
        for b in order:
            if b in blocks:
                f.write('<%s>\n' % b)
                for k, v in blocks[b]:
                    f.write('%s = %s\n' % (k, v))
                f.write('\n')


def deck_drop_outputs(deck):
    blocks, order = deck
    for b in list(order):
        if re.match(r'output\d+$', b):
            del blocks[b]
            order.remove(b)


# ---------------------------------------------------------------- running
class Cores:
    """counting slots.  wallwb-1105: a run wider than the pool takes the whole pool (it runs
    oversubscribed) instead of waiting for ever: with --np-max 4 and a 6-rank layout, take(6)
    never returned and the suite hung silently for 11 h with no athena running"""
    def __init__(self, n):
        self.n, self.total, self.c = n, n, threading.Condition()

    def take(self, k):
        k = min(k, self.total)
        with self.c:
            while self.n < k:
                self.c.wait()
            self.n -= k

    def give(self, k):
        k = min(k, self.total)
        with self.c:
            self.n += k
            self.c.notify_all()


def build_deck(case, lay, perturb, rundir):
    deck = read_deck(os.path.join(HERE, 'decks', case['deck']))
    deck_drop_outputs(deck)
    for k, v in case.get('set', {}).items():
        deck_set(deck, k, v)
    mb = lay.get('mb')
    if mb is not None:
        for d, n in zip(('nx1', 'nx2', 'nx3'), mb):
            deck_set(deck, 'meshblock/' + d, n)
    for k, v in lay.get('set', {}).items():
        deck_set(deck, k, v)
    ncyc = case.get('ncyc', 10)
    deck_set(deck, 'time/nlim', ncyc)
    deck_set(deck, 'time/tlim', 1e30)
    if perturb and case.get('tolref'):
        # layinv-1105 (USER 11-05 option a): the reference spread of a deck whose solver
        # stops at a tolerance with a block-local preconditioner (mg / rbgs): the same
        # layout with the tolerances scaled (e.g. x 0.5), instead of the 1-ulp perturbation
        for pk, f in case['tolref'].items():
            deck_set(deck, pk, repr(float(deck_get(deck, pk)) * f))
    elif perturb:
        pk = case.get('perturb', 'time/cfl_number')
        x = float(deck_get(deck, pk))
        deck_set(deck, pk, repr(float(np.nextafter(x, np.inf))))
    for n, var in enumerate(case['vars']):
        b = 'output%d' % (n + 1)
        for k, v in (('file_type', 'bin'), ('variable', var), ('dcycle', 1),
                     ('bin_double', 'true'), ('file_id', 'li%d' % n)):
            deck_set(deck, b + '/' + k, v)
    deck_write(deck, os.path.join(rundir, 'run.athinput'))


def run_one(case, lay, perturb, binp, work, cores, opts):
    tag = lay['name'] + ('_ulp' if perturb else '')
    rd = os.path.join(work, case['name'], tag)
    if opts.compare_only:
        return rd, 'cached'
    if os.path.isdir(rd):
        shutil.rmtree(rd)
    os.makedirs(rd)
    for src, dst in case.get('links', {}).items():
        os.symlink(src, os.path.join(rd, dst))
    build_deck(case, lay, perturb, rd)
    np_ = lay.get('np', 1)
    cmd = ['nice', '-n', '10']
    if os.environ.get('LI_SRUN'):
        # layinv-1105: launch through srun into a held allocation with a memory cap, e.g.
        # LI_SRUN="srun --overlap --jobid=5719508 -w hpc-sm-01-05 -N1 -c1 --mpi=pmix"
        # (LI_MEM_GB per rank, default 6)
        mem = int(os.environ.get('LI_MEM_GB', '6'))*np_ + 4
        cmd += os.environ['LI_SRUN'].split() + ['-n', str(np_), '--mem=%dG' % mem]
    elif np_ > 1 or opts.always_mpirun:
        cmd += ['mpirun', '-np', str(np_), '--bind-to', 'none', '--oversubscribe']
    cmd += [binp, '-i', 'run.athinput']
    env = dict(os.environ, OMP_NUM_THREADS='1')
    env.update(case.get('env', {}))
    cores.take(np_)
    t0 = time.time()
    try:
        with open(os.path.join(rd, 'run.log'), 'w') as lf:
            r = subprocess.run(cmd, cwd=rd, stdout=lf, stderr=subprocess.STDOUT, env=env,
                               timeout=opts.run_timeout or case.get('timeout', 3600))
        st = 'rc=%d %.0fs' % (r.returncode, time.time() - t0)
    except subprocess.TimeoutExpired:
        st = 'TIMEOUT'
    finally:
        cores.give(np_)
    return rd, st


# ---------------------------------------------------------------- comparing
GROUPS = [('mom1', 'mom2', 'mom3'), ('bcc1', 'bcc2', 'bcc3'), ('m1_f1', 'm1_f2', 'm1_f3'),
          ('F1', 'F2', 'F3'), ('velx', 'vely', 'velz')]


def scales(B):
    """normalisation per variable: max |.| of the variable, or of its vector group (a
    component that is ~0 everywhere, e.g. v_phi of a blast at rest, is then measured
    against the vector's magnitude, not against its own round-off)"""
    sc = {v: float(np.nanmax(np.abs(B[v]))) for v in B}
    for g in GROUPS:
        m = max([sc[v] for v in g if v in sc] or [0.0])
        for v in g:
            if v in sc:
                sc[v] = m
    return sc


def reldiff(a, b, sc=None):
    d = np.abs(a - b)
    if sc is None:
        sc = np.nanmax(np.abs(a))
    if not np.isfinite(d).all():
        return np.inf, d
    return (d.max() / sc if sc > 0 else d.max()), d


def locate(d, lay, ng):
    """where the largest differences sit relative to the INTERNAL block ends of `lay`"""
    idx = np.unravel_index(np.argmax(d), d.shape)   # (panel, k, j, i)
    n1, n2, n3 = lay['n']
    N = d.shape[1:][::-1]   # (N1, N2, N3)
    nz = np.argwhere(d > 0)
    near = 0
    if len(nz):
        def nearend(c, n, Nn):
            r = c % n
            interior = (Nn // n) > 1
            return interior & ((r < ng) | (r >= n - ng))
        m = nearend(nz[:, 3], n1, N[0]) | nearend(nz[:, 2], n2, N[1]) | nearend(nz[:, 1], n3, N[2])
        near = m.mean()
    return 'max@p%d k%d j%d i%d; %d cells differ, %.0f%% within ng of an internal block end' % (
        idx[0], idx[1], idx[2], idx[3], len(nz), 100 * near)


def verdict(d, ref, tol=None):
    if d == 0:
        return 'BITWISE'
    if tol is not None:   # layinv-1105 tolerance gate: within tol x the tolerance spread
        return 'TOLERANCE' if d <= max(tol * ref, 100 * EPS) else 'FAIL'
    if d <= max(10 * ref, 100 * EPS):
        return 'ROUNDOFF'
    if d > max(1000 * ref, 1e-10):
        return 'FAIL'
    return 'SUSPECT'


def compare_case(case, work, out):
    base = case['layouts'][0]['name']
    bd = os.path.join(work, case['name'], base)
    ud = bd + '_ulp'
    fc = libin.files_by_cycle(bd)
    if not fc:
        out.append('  NO OUTPUT in %s (see run.log)' % bd)
        return 'ERROR'
    ncyc = max(fc)
    cycs = sorted(set([1, ncyc]) & set(fc))
    ng = int(case.get('ng', 2))
    worst = 'BITWISE'
    rank = ['BITWISE', 'ROUNDOFF', 'TOLERANCE', 'NOWORSE', 'SUSPECT', 'FAIL', 'ERROR']
    tolg = case.get('tol_factor', 1.0) if case.get('tolref') else None
    for c in cycs:
        B, _ = libin.load_cycle(bd, c)
        U, _ = libin.load_cycle(ud, c)
        refs = {}
        sc = scales(B)
        for v in B:
            refs[v] = reldiff(B[v], U[v], sc[v])[0] if v in U else np.nan
        out.append('  cycle %d  %s: %s' % (c, 'tol-ref' if tolg else 'ulp-ref', ' '.join('%s %.1e' % (v, refs[v]) for v in B)))
        for lay in case['layouts'][1:]:
            ld = os.path.join(work, case['name'], lay['name'])
            try:
                L, ll = libin.load_cycle(ld, c)
            except Exception as e:   # noqa: BLE001
                out.append('    %-14s ERROR %s' % (lay['name'], e))
                worst = 'ERROR'
                continue
            if not L:
                out.append('    %-14s ERROR no output at cycle %d' % (lay['name'], c))
                worst = 'ERROR'
                continue
            parts, wl, wd = [], 'BITWISE', None
            for v in B:
                if v not in L or L[v].shape != B[v].shape:
                    parts.append('%s SHAPE' % v)
                    wl = 'ERROR'
                    continue
                r, d = reldiff(B[v], L[v], sc[v])
                vd = verdict(r, refs[v], tolg)
                if case.get('ref_old') and vd != 'BITWISE':
                    # layinv-1105 (a1, USER 11-05): a block-local preconditioner (mg / rbgs):
                    # PASS when the layout pair is no worse (x2) than the OLD binary's pair
                    try:
                        ob, _ = libin.load_cycle(os.path.join(work, case['name'], 'old_' + base), c)
                        ol, _ = libin.load_cycle(os.path.join(work, case['name'],
                                                              'old_' + lay['name']), c)
                        ro = reldiff(ob[v], ol[v], sc[v])[0]
                        vd = 'NOWORSE' if r <= max(2.0 * ro, 100 * EPS) else 'FAIL'
                        r = (r, ro)
                    except Exception:   # noqa: BLE001
                        vd = 'ERROR'
                if isinstance(r, tuple):
                    parts.append('%s %.1e(old %.1e)%s' % (v, r[0], r[1], '' if vd == 'NOWORSE'
                                                         else '!'))
                else:
                    parts.append('%s %.1e%s' % (v, r, '' if vd in ('BITWISE', 'ROUNDOFF',
                                                                  'TOLERANCE') else '!'))
                if rank.index(vd) > rank.index(wl):
                    wl, wd = vd, (v, d)
            line = '    %-14s %-8s %s' % (lay['name'], wl, ' '.join(parts))
            if wd is not None and wl in ('SUSPECT', 'FAIL'):
                line += '\n                            %s: %s' % (wd[0], locate(wd[1], ll, ng))
            out.append(line)
            if rank.index(wl) > rank.index(worst):
                worst = wl
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bins', default=os.path.join(HERE, 'bins.cfg'))
    ap.add_argument('--work', default=os.path.join(os.getcwd(), 'li_runs'))
    ap.add_argument('--np-max', type=int, default=16)
    ap.add_argument('--only', default='')
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--compare-only', action='store_true')
    ap.add_argument('--always-mpirun', action='store_true')
    ap.add_argument('--run-timeout', type=int, default=0,
                    help='seconds per run, overrides the case timeouts (TIMEOUT -> ERROR in the report)')
    ap.add_argument('--allow-float32', action='store_true',
                    help='compare float32 .bin dumps anyway (default: refuse)')
    opts = ap.parse_args()
    libin.ALLOW_FLOAT32 = opts.allow_float32 or libin.ALLOW_FLOAT32
    bins = {}
    if os.path.exists(opts.bins):
        for s in open(opts.bins):
            s = s.split('#', 1)[0].strip()
            if '=' in s:
                k, v = s.split('=', 1)
                bins[k.strip()] = v.strip()
    sel = [c for c in CASES.CASES if not opts.only or c['name'] in opts.only.split(',')]
    if opts.list:
        for c in sel:
            print('%-24s %-4s %s | %s' % (c['name'], c['bin'], c['deck'],
                                         ' '.join(l['name'] for l in c['layouts'])))
        return 0
    os.makedirs(opts.work, exist_ok=True)
    cores = Cores(opts.np_max)
    jobs, skipped = [], []
    with ThreadPoolExecutor(max_workers=opts.np_max) as ex:
        for c in sel:
            bp = bins.get(c['bin'])
            if not bp or not os.path.exists(bp):
                skipped.append(c['name'])
                continue
            ob = bins.get(c['bin'] + '_old') if c.get('ref_old') else None
            for i, lay in enumerate(c['layouts']):
                jobs.append((c, lay, ex.submit(run_one, c, lay, False, bp, opts.work, cores, opts)))
                if ob and os.path.exists(ob):   # layinv-1105 (a1): the OLD binary's same layouts
                    lo = dict(lay, name='old_' + lay['name'])
                    jobs.append((c, lo, ex.submit(run_one, c, lo, False, ob, opts.work, cores, opts)))
                if i == 0:
                    jobs.append((c, lay, ex.submit(run_one, c, lay, True, bp, opts.work, cores, opts)))
        status = {}
        for c, lay, fu in jobs:
            rd, st = fu.result()
            status.setdefault(c['name'], []).append('%s:%s' % (os.path.basename(rd), st))
    summary = []
    rep = []
    for c in sel:
        if c['name'] in skipped:
            summary.append('%-24s SKIPPED (no binary for family %s)' % (c['name'], c['bin']))
            continue
        out = ['== %s  (%s)  %s' % (c['name'], c.get('what', ''), ' '.join(status.get(c['name'], [])))]
        try:
            w = compare_case(c, opts.work, out)
        except Exception as e:   # noqa: BLE001
            out.append('  COMPARE ERROR %r' % e)
            w = 'ERROR'
        exp = c.get('expect', 'PASS')
        ok = w in ('BITWISE', 'ROUNDOFF', 'TOLERANCE', 'NOWORSE')
        if c.get('require_bitwise'):   # layinv-1105: layout_invariant cases
            ok = (w == 'BITWISE')
        tagx = '' if (ok and exp == 'PASS') or (not ok and exp != 'PASS') else '  <-- UNEXPECTED'
        summary.append('%-24s %-8s expect %s%s' % (c['name'], w, exp, tagx))
        rep += out
    txt = '\n'.join(rep + ['', 'SUMMARY'] + summary)
    print(txt)
    with open(os.path.join(opts.work, 'REPORT.txt'), 'w') as f:
        f.write(txt + '\n')
    return 0 if all('UNEXPECTED' not in s and 'ERROR' not in s for s in summary) else 1


if __name__ == '__main__':
    sys.exit(main())
