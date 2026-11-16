#!/usr/bin/env python3
"""compare.py [--runs DIR] [group ...]: GPU runs of runs.txt against the Orion CPU reference.

Per run (S1-S3): rc, FATAL/NaN, CLEAN WARNING, DT_RETRY, CLEAN cycle lines and their solves,
and the final dump of every bin and tab output (ref/dumps/<tag>/) as the field maxrel
max|gpu - cpu| / max|cpu| over all variables.  Gates:
  - rc 0, no FATAL/NaN, CLEAN WARNING and DT_RETRY no more than the CPU's (0 everywhere);
  - CLEAN cycle lines within 10 % of the CPU's (the solve sequence is the same);
  - field maxrel <= max(1e-6, 10 x the CPU's own 1- vs 4-thread maxrel, ref/ref.json).
Prints one table and "all gates passed" or the list of GATE FAILURES."""
import glob, json, os, re, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, H + '/ana')
import bin_convert as bc

args = sys.argv[1:]
runs = H + '/runs'
if args[:1] == ['--runs']:
    runs, args = os.path.abspath(args[1]), args[2:]
groups = args or ['S1', 'S2', 'S3']
ref = json.load(open(H + '/ref/ref.json'))


def logstat(f):
    t = open(f, errors='replace').read()
    rc = re.findall(r'^rc=(\d+)', t, re.M)
    cl = re.findall(r'^CLEAN cycle .* solves (\d+)', t, re.M)
    return dict(rc=int(rc[-1]) if rc else -1,
                bad=len(re.findall(r'### FATAL|(^|[^A-Za-z_])-?(nan|NaN)([^A-Za-z_]|$)', t)),
                warn=t.count('CLEAN WARNING'), retry=len(re.findall(r'^DT_RETRY', t, re.M)),
                lines=len(cl), solves=sum(int(x) for x in cl))


def last_bins(d):
    out = {}
    for f in glob.glob(d + '/bin/*.bin') + glob.glob(d + '/tab/*.tab'):
        stem, num = os.path.basename(f).rsplit('.', 2)[0], os.path.basename(f).rsplit('.', 2)[1]
        if stem not in out or num > out[stem][0]:
            out[stem] = (num, f)
    return {k: v[1] for k, v in out.items()}


def maxrel(a, b):
    grp = lambda k: re.sub(r'(m1_f|mom|m1_vh|bcc|vel)[123xyz]?$', r'\1', k)
    if a.endswith('.tab'):   # (formatted tables: the data columns after gid, i, x1v)
        p, q = np.loadtxt(a, comments='#', ndmin=2), np.loadtxt(b, comments='#', ndmin=2)
        if p.shape != q.shape:
            return float('inf')
        nd = q.shape[1] - 3
        names = open(b).readlines()[1].split()[1:][-nd:]
        sc = {}
        for n, c in zip(names, range(3, q.shape[1])):
            sc[grp(n)] = max(sc.get(grp(n), 0.0), float(np.abs(q[:, c]).max()))
        r = 0.0
        for n, c in zip(names, range(3, q.shape[1])):
            s = sc[grp(n)]
            r = max(r, float(np.abs(p[:, c] - q[:, c]).max()/s) if s > 0 else 0.0)
        return r
    x, y = bc.read_binary(a), bc.read_binary(b)
    # (vector components share one scale, the largest component: a component that is 0 up
    # to round-off is not compared with itself)
    sc = {}
    for k in y['var_names']:
        sc[grp(k)] = max(sc.get(grp(k), 0.0), float(np.abs(np.array(y['mb_data'][k])).max()))
    r = 0.0
    for k in y['var_names']:
        p, q = np.array(x['mb_data'][k]), np.array(y['mb_data'][k])
        s = sc[grp(k)]
        r = max(r, float(np.abs(p - q).max()/s) if s > 0 else float(np.abs(p).max()))
    return r


fails = []
print('| run | rc | warn (cpu) | DT_RETRY (cpu) | CLEAN lines (cpu) | solves (cpu) | field maxrel | gate |')
print('|---|---|---|---|---|---|---|---|')
for ln in open(H + '/runs.txt'):
    w = ln.split()
    if not w or w[0].startswith('#') or w[0] not in groups:
        continue
    tag = w[1]
    d = '%s/%s' % (runs, tag)
    if not os.path.exists(d + '/run.log'):
        print('| %s | missing |' % tag); fails.append(tag + ': missing'); continue
    s, c = logstat(d + '/run.log'), ref[tag]
    mr, miss = 0.0, []
    rb = last_bins('%s/ref/dumps/%s' % (H, tag))
    gb = last_bins(d)
    for k, f in rb.items():
        if k not in gb or os.path.basename(gb[k]) != os.path.basename(f):
            miss.append(k); continue
        mr = max(mr, maxrel(gb[k], f))
    gate = max(1e-6, 10.0*c['spread'])
    bad = []
    if s['rc'] != 0 or s['bad']: bad.append('rc/FATAL/NaN')
    if s['warn'] > c['warn'] or s['retry'] > c['retry']: bad.append('warnings/retries')
    if abs(s['lines'] - c['lines']) > 0.1*c['lines'] + 1: bad.append('CLEAN lines')
    if miss: bad.append('dumps missing: ' + ','.join(miss))
    if mr > gate: bad.append('maxrel %.1e > %.1e' % (mr, gate))
    print('| %s | %d | %d (%d) | %d (%d) | %d (%d) | %d (%d) | %.2e | %.1e %s |' % (
        tag, s['rc'], s['warn'], c['warn'], s['retry'], c['retry'], s['lines'], c['lines'],
        s['solves'], c['solves'], mr, gate, 'ok' if not bad else 'FAIL'))
    fails += ['%s: %s' % (tag, b) for b in bad]
print()
print('all gates passed' if not fails else 'GATE FAILURES:\n' + '\n'.join(fails))
