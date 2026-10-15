#!/usr/bin/env python3
"""ridort-1014 reval.py <runs root> [family glob ...]: R-vs-REF difference of every run dir, set against that family's OWN
ladder differences (same deck family, neighbouring level).

Run dirs: <root>/<kind>/<case>/<name>/{ref,r,o} (env ARM = r | o picks the candidate arm); the level is the trailing _<int> of <name> (sp HFfZ_1841da_<l>, lad
FIN_hesdirk2_{T,X,C}_<l>, cart FIN3_n512 / _r4 ... -> no ladder unless the names differ only in a trailing number).
Data: last dump of every output stream (tab/*.tab columns, bin/*.bin variables via bin_convert, BINCONV=<athenak>/vis/python).
Metric per variable: L1rel = sum|a - b| / sum|ref| (same normalisation for the R diff and the ladder diff, so the ratio is
metric-consistent even for ~0-mean variables such as velocities in static tests).
  D_R(l)   = L1rel(R_l, REF_l)
  D_lad(l) = L1rel(REF_l, REF_l') with l' the next finer level (else next coarser); same grid -> direct (T ladders);
             1-D tabs on different grids -> the finer interpolated (np.interp in x1v) onto the coarser (X / C ladders);
             bin on different grids -> skipped (rlat X: use the IDORT evaluator).
  ratio    = D_R / D_lad.   Per dir prints the worst ratio over 'energy-type' variables (dens, eint, press, m1_e, E, T, r_e ...)
             and over the variables the ladder resolves (D_lad < 0.3; 'noise vars' = the others, e.g. transverse
             velocities / fluxes of 1-D tests, whose L1rel between ANY two runs is O(1)), with the variable names.
Output lines:  <dir> | lev | maxrel_any | E: var D_R D_lad ratio | res: var D_R D_lad ratio | noise vars n | Picard REF/R
"""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.environ.get('BINCONV', ''))
ROOT = sys.argv[1]
ARM = os.environ.get('ARM', 'r')   # r = option R, o = ROL2U
PAT = sys.argv[2:] or ['*/*']
SKIPC = {'gid', 'i', 'j', 'k', 'x1v', 'x2v', 'x3v'}
ETYPE = re.compile(r'^(dens|eint|press|pgas|m1_e|e|E|tgas|temp|T|r_e|er|rho|ener|eng)$', re.I)


def last_dumps(d):
    """stream -> last file"""
    out = {}
    for sub in ('tab', 'bin'):
        for f in sorted(glob.glob(os.path.join(d, sub, '*.' + sub))):
            b = os.path.basename(f).split('.')
            out['.'.join(b[:-2]) + '.' + sub] = f
    return out


_cache = {}


def load(f):
    if f in _cache:
        return _cache[f]
    if f.endswith('.tab'):
        names = None
        with open(f) as h:
            for ln in h:
                if ln.startswith('#') and ('x1v' in ln or ' i ' in ln):
                    names = ln[1:].split()
                if not ln.startswith('#'):
                    break
        a = np.loadtxt(f, comments='#', ndmin=2)
        for drop in ('k', 'x3v', 'j', 'x2v'):   # athenak tab headers may list k/x3v columns the data lacks (cf. od.py)
            if len(names) > a.shape[1] and drop in names:
                names.remove(drop)
        r = {'x1v': a[:, names.index('x1v')] if 'x1v' in names else None}
        r['vars'] = {n: a[:, c] for c, n in enumerate(names) if n not in SKIPC and c < a.shape[1]}
        r['grid'] = a.shape[0]
    else:
        import bin_convert
        x = bin_convert.read_binary(f)
        r = {'x1v': None, 'vars': {v: np.concatenate([np.ravel(q) for q in x['mb_data'][v]]) for v in x['var_names']}}
        r['grid'] = tuple(r['vars'][x['var_names'][0]].shape)
    _cache[f] = r
    return r


def l1(a, b, ref):
    s = np.sum(np.abs(ref))
    return np.sum(np.abs(a - b)) / s if s > 0 else (0.0 if np.sum(np.abs(a - b)) == 0 else np.inf)


def diffs(da, db, dref=None):
    """var -> L1rel(a, b) normalised by a (dref: dir whose data normalises); interpolates 1-D tabs if grids differ"""
    out = {}
    A, B = last_dumps(da), last_dumps(db)
    for s in A:
        if s not in B:
            continue
        x, y = load(A[s]), load(B[s])
        for v, a in x['vars'].items():
            if v not in y['vars']:
                continue
            b = y['vars'][v]
            if a.shape != b.shape:
                if x['x1v'] is None or y['x1v'] is None:
                    continue
                if len(np.unique(x['x1v'])) != len(x['x1v']) or len(np.unique(y['x1v'])) != len(y['x1v']):
                    continue   # multi-d tab: skip
                if len(a) < len(b):
                    b = np.interp(x['x1v'], y['x1v'], b)
                else:
                    a = np.interp(y['x1v'], x['x1v'], a)
                    out[s + ':' + v] = l1(a, b, b)
                    continue
            out[s + ':' + v] = l1(a, b, a)
    return out


def maxrel(da, db):
    m = 0.0
    A, B = last_dumps(da), last_dumps(db)
    for s in A:
        if s in B:
            x, y = load(A[s]), load(B[s])
            for v, a in x['vars'].items():
                if v in y['vars'] and a.shape == y['vars'][v].shape:
                    b = y['vars'][v]
                    sc = np.maximum(np.abs(a), np.abs(b))
                    m = max(m, float(np.max(np.where(sc > 0, np.abs(a - b) / np.where(sc > 0, sc, 1), 0))))
    return m


def pic(d):
    for f in ('all.log', 'run.log'):
        p = os.path.join(d, f)
        if os.path.exists(p):
            t = open(p, errors='replace').read()
            m = re.findall(r'Picard iterations mean=(\S+)', t)
            return float(m[-1]) if m else float('nan')
    return float('nan')


def rc(d):
    try:
        return open(os.path.join(d, 'rc.txt')).read().strip()
    except OSError:
        return '-'


fams = {}
for pat in PAT:
    for d in sorted(glob.glob(os.path.join(ROOT, pat))):
        if not (os.path.isdir(os.path.join(d, 'ref')) and os.path.isdir(os.path.join(d, ARM))):
            continue
        name = os.path.basename(d)
        m = re.match(r'^(.*?)(\d+)$', name)
        key = (os.path.dirname(d), m.group(1) if m else name)
        fams.setdefault(key, []).append((int(m.group(2)) if m else 0, d))
for (case, pre), lst in sorted(fams.items()):
    lst.sort()
    for i, (lev, d) in enumerate(lst):
        rr, rf = os.path.join(d, ARM), os.path.join(d, 'ref')
        if rc(rr) != '0' or rc(rf) != '0':
            print('%s | rc=%s/%s FAIL' % (os.path.relpath(d, ROOT), rc(rf), rc(rr)))
            continue
        dr = diffs(rf, rr)
        dl = {}
        nb = lst[i + 1][1] if i + 1 < len(lst) else (lst[i - 1][1] if i > 0 else None)
        if nb is not None and rc(os.path.join(nb, 'ref')) == '0':
            dl = diffs(rf, os.path.join(nb, 'ref'))
        best_e, best_a, nnoise = None, None, 0
        for v, x in dr.items():
            y = dl.get(v)
            q = (x / y) if (y is not None and y > 0) else (np.inf if x > 0 else 0.0)
            if y is None:
                q = float('nan')
            t = (q if q == q else -1, v, x, y)
            vn = v.split(':')[-1]
            if y is not None and y < 0.3 and (best_a is None or t[0] > best_a[0]):
                best_a = t   # 'res': variables the ladder itself resolves (D_lad < 0.3; excludes ~0-mean noise fields)
            if y is not None and y >= 0.3:
                nnoise += 1
            if ETYPE.match(vn) and (best_e is None or t[0] > best_e[0]):
                best_e = t
        f = lambda t: ('-' if t is None else '%s %.2e %s %s' % (t[1], t[2], '-' if t[3] is None else '%.2e' % t[3],
                                                                   '-' if t[0] < 0 else '%.2e' % t[0]))
        print('%s | l=%d | maxrel %.2e | E: %s | res: %s | noise vars %d | Picard %.3f/%.3f' % (
            os.path.relpath(d, ROOT), lev, maxrel(rf, rr), f(best_e), f(best_a), nnoise, pic(rf), pic(rr)))
        sys.stdout.flush()
    _cache.clear()
