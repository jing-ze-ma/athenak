#!/usr/bin/env python3
"""resf.py D NH (FRESHPREP): fresh-start smoke summary. a = straight, b = NH cycles + restart -> end."""
import sys, glob, os, re
D, NH = sys.argv[1], int(sys.argv[2])
def lines(f): return open(f, errors='ignore').read().splitlines()
L = lines(D + '/a/run.log')
el = [(float(m.group(1)), int(m.group(2)), float(m.group(3))) for l in L for m in [re.search(r'elapsed=(\S+) cycle=(\d+) time=\S+ dt=(\S+)', l)] if m]
if len(el) > 2:
    e = [x[0] for x in el]; d = [e[i+1]-e[i] for i in range(len(e)-1)]; st = d[1:-1] if len(d) > 2 else d
    print('a cycles', el[0][1], '-', el[-1][1], 'dt first %.6g last %.6g' % (el[0][2], el[-1][2]), 's/cycle steady %.3f' % (sum(st)/len(st)), 'first elapsed %.1f' % e[0])
for k in ['FATAL', 'CLEAN WARNING', 'DT_RETRY', 'WARNING', 'march', 'clean_init steady:', 'he_star_m1: GM', 'he_seed_band = n2 BAND', 'seed = ', 'direction maps', 'he_ic_equil']:
    hits = [l for l in L if k in l]; print('a', repr(k), len(hits), (hits[-1][:300] if hits else ''))
for x in sorted(glob.glob(D + '/b/run*.log')):
    B = lines(x); print('b', os.path.basename(x), 'CW', sum('CLEAN WARNING' in l for l in B), 'DT_RETRY', sum('DT_RETRY' in l for l in B), 'FATAL', sum('FATAL' in l for l in B))
mem = {}
if os.path.exists(D + '/mem.log'):
    for l in open(D + '/mem.log', errors='ignore'):
        p = l.replace(',', ' ').split()
        if len(p) == 3 and p[2].isdigit(): k = p[0]+':'+p[1]; mem[k] = max(mem.get(k, 0), int(p[2]))
if mem: print('GPU mem peak MiB max %d min %d (n=%d) = %.1f %% of 143771 MiB' % (max(mem.values()), min(mem.values()), len(mem), 100*max(mem.values())/143771.))
def rows(f): return [l for l in open(f) if not l.startswith('#')]
for h in sorted(glob.glob(D + '/a/*.hst')):
    b = os.path.join(D, 'b', os.path.basename(h)); ra = rows(h); rb = rows(b) if os.path.exists(b) else []
    ta = {l.split()[0]: l for l in ra}; nb = [l for l in rb if l.split()[0] in ta]
    tail_b = rb[-(len(ra)-NH-1):] if len(ra) > NH+1 else []
    same = len(tail_b) > 0 and ra[-len(tail_b):] == tail_b
    print('rt', os.path.basename(h), 'rows a', len(ra), 'b', len(rb), 'last %d rows' % len(tail_b), 'IDENTICAL' if same else 'DIFF')
def data(f):
    r = open(f, 'rb').read(); i = r.find(b'<par_end>'); return r[i:] if i >= 0 else r
last = {}
for f in sorted(glob.glob(D + '/a/bin/*.bin')): last[f.split('.')[-3]] = f
for k, f in sorted(last.items()):
    g = os.path.join(D, 'b', 'bin', os.path.basename(f)); print('rt', os.path.basename(f), 'same' if os.path.exists(g) and data(f) == data(g) else 'DIFF/missing')
