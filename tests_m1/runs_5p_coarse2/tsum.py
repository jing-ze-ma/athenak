"""usage: python3 tsum.py <runs dir> [c0=10] [c1=30]
ms/cycle over cycles c0..c1 (elapsed), BiCGStab iterations per solve over steps >= c0
(implicit_picard_log), NON-CONVERGED; per run, then the mean per arm."""
import collections
import os
import re
import sys

R = sys.argv[1]
c0 = int(sys.argv[2]) if len(sys.argv) > 2 else 10
c1 = int(sys.argv[3]) if len(sys.argv) > 3 else 30
arms = collections.OrderedDict()
for a in sorted(os.listdir(R)):
    p = os.path.join(R, a, 'run.log')
    if not os.path.exists(p):
        continue
    t = {}
    nin = []
    nc = 'nan'
    for ln in open(p):
        m = re.search(r'elapsed=([0-9.e+-]+) cycle=(\d+)', ln)
        if m:
            t[int(m.group(2))] = float(m.group(1))
        m = re.search(r'plog step=(\d+) pass=\d+ .* nin=(\d+)', ln)
        if m and int(m.group(1)) >= c0:
            nin.append(int(m.group(2)))
        m = re.search(r'NON-CONVERGED=([0-9.e+-]+)', ln)
        if m:
            nc = m.group(1)
    ms = (t[c1] - t[c0])/(c1 - c0)*1e3 if (c0 in t and c1 in t) else float('nan')
    solves = [n for n in nin if n > 0]
    ips = sum(solves)/max(len(solves), 1)
    ipst = sum(nin)/max(c1 - c0 + 1, 1)
    name = a.rsplit('_', 1)[0]
    arms.setdefault(name, []).append((ms, ips, max(solves or [0]), ipst))
    print(f'{a:16s} {ms:7.2f} ms/cycle  it/solve {ips:5.1f} max {max(solves or [0]):3d}'
          f'  it/step {sum(nin)/max(1, len(set(range(c0, 41)))):5.1f}  NONCONV {nc}')
print('arm              ms/cycle (reps)          it/solve  max')
for k, v in arms.items():
    print(f'{k:16s} ' + ', '.join(f'{x[0]:.2f}' for x in v).ljust(22)
          + f'  {sum(x[1] for x in v)/len(v):5.1f}  {max(x[2] for x in v):3d}')
