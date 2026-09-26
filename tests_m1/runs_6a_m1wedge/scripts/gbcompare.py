#!/usr/bin/env python3
"""hipgate_0926: bitwise comparison of the HIP gate runs.
  python3 compare.py [prob ...]
(1) every output file of each arm's fresh run vs the reference arm old46 (cmp of whole files);
(2) restart: the rr run's last rst vs the fresh run's last rst (bytes after <par_end>), and
    the rr hst rows vs the fresh hst rows at the same cycles/times."""
import os
import sys
import glob
import numpy as np

P = '/viper/ptmp2/jinma/sprhd_0926/gpu/gb'
REF = 'base'


def files(d):
    out = []
    for root, _, fs in os.walk(d):
        for f in fs:
            if f.endswith('.log') or f.startswith('column_used'):
                continue
            out.append(os.path.relpath(os.path.join(root, f), d))
    return sorted(out)


def payload(fn):
    b = open(fn, 'rb').read()
    i = b.find(b'<par_end>')
    return b[i:] if i >= 0 else b


def hst_rows(fn):
    return [ln for ln in open(fn) if not ln.startswith('#')]


def hst_diff(a, b):
    ra, rb = hst_rows(a), hst_rows(b)
    n = min(len(ra), len(rb))
    for k in range(n):
        if ra[k] != rb[k]:
            x = np.array(ra[k].split(), float)
            y = np.array(rb[k].split(), float)
            rel = np.abs(x - y) / np.maximum(np.abs(x), 1e-300)
            return 'first diff row %d (of %d/%d), cols %s, max rel %.3e' % (
                k, len(ra), len(rb), list(np.nonzero(rel)[0]), rel.max())
    return 'rows equal (%d vs %d rows)' % (len(ra), len(rb))


def bytes_diff(a, b):
    x, y = open(a, 'rb').read(), open(b, 'rb').read()
    if len(x) != len(y):
        return 'size %d vs %d' % (len(x), len(y))
    xa = np.frombuffer(x, np.uint8)
    ya = np.frombuffer(y, np.uint8)
    idx = np.nonzero(xa != ya)[0]
    return '%d bytes differ, first at %d of %d' % (len(idx), idx[0], len(x))


probs = sys.argv[1:] or sorted(os.listdir(P))
for pb in probs:
    arms = sorted(os.listdir(os.path.join(P, pb)))
    print('==== %s  arms %s' % (pb, arms))
    rd = os.path.join(P, pb, REF, 'fresh')
    ref_files = files(rd) if os.path.isdir(rd) else []
    for arm in arms:
        if arm == REF:
            continue
        ad = os.path.join(P, pb, arm, 'fresh')
        af = files(ad)
        same, diff = [], []
        for f in ref_files:
            g = os.path.join(ad, f)
            if not os.path.exists(g):
                diff.append((f, 'missing'))
            elif open(os.path.join(rd, f), 'rb').read() == open(g, 'rb').read():
                same.append(f)
            else:
                why = (hst_diff(os.path.join(rd, f), g) if f.endswith('.hst')
                       else bytes_diff(os.path.join(rd, f), g))
                diff.append((f, why))
        extra = sorted(set(af) - set(ref_files))
        tag = 'BITWISE' if not diff and not extra else 'NOT BITWISE'
        print('  fresh %s vs %s: %s  (%d same, %d diff, extra %s)'
              % (arm, REF, tag, len(same), len(diff), extra))
        for f, why in diff:
            print('     DIFF %s: %s' % (f, why))
    for arm in arms:
        fd = os.path.join(P, pb, arm, 'fresh')
        rr = os.path.join(P, pb, arm, 'rr')
        fr = sorted(glob.glob(fd + '/rst/*.rst'))
        rs = sorted(glob.glob(rr + '/rst/*.rst'))
        if not fr or not rs:
            print('  restart %s: missing rst (fresh %d, rr %d)' % (arm, len(fr), len(rs)))
            continue
        ok = payload(fr[-1]) == payload(rs[-1])
        h1 = glob.glob(fd + '/*.hst')
        h2 = glob.glob(rr + '/*.hst')
        hs = ''
        if h1 and h2:
            r1 = hst_rows(h1[0])
            r2 = hst_rows(h2[0])
            common = [r for r in r2 if r in r1]
            hs = ', hst rr rows %d, found verbatim in fresh %d' % (len(r2), len(common))
        print('  restart %s: last rst %s vs %s payload %s%s' % (
            arm, os.path.basename(rs[-1]), os.path.basename(fr[-1]),
            'BITWISE' if ok else 'DIFF', hs))
