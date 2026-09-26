# usage: lev.py <g5 dir>: ms/cycle, wall per simulated second (cycles 20..end of each arm),
# Picard and BiCGStab inner iterations per solve, NON-CONVERGED
import sys, glob, re, os, numpy as np
R = sys.argv[1]; res = {}
for d in sorted(glob.glob(R + '/*_r[123]')):
    a, rep = os.path.basename(d).rsplit('_', 1)
    L = open(d + '/run.log').read()
    rows = [tuple(float(x) for x in m) for m in re.findall(r'elapsed=(\S+) cycle=(\d+) time=(\S+)', L)]
    if len(rows) < 3: print(a, rep, 'FAILED'); continue
    (e0, c0, t0), (e1, c1, t1) = rows[1], rows[-1]
    pic = re.search(r'Picard iterations mean=(\S+)', L); inn = re.search(r'inner iterations mean=(\S+)', L)
    nc = re.search(r'NON-CONVERGED=(\S+)', L)
    res.setdefault(a, []).append((1e3*(e1-e0)/(c1-c0), (e1-e0)/(t1-t0), float(pic.group(1)),
                                  float(inn.group(1)), float(nc.group(1))))
ref = np.mean([r[1] for r in res['ref']])
print('%-9s %16s %22s %8s %8s %4s' % ('arm', 'ms/cycle r1/r2', 'wall/sim-s r1/r2', 'vs ref', 'Picard', 'inner'))
for a, v in res.items():
    w = np.mean([r[1] for r in v])
    print('%-9s %7.2f / %6.2f %10.4f / %8.4f %+7.1f%% %8.3f %8.3f NC=%g' % (a, v[0][0], v[-1][0], v[0][1], v[-1][1],
          100*(w/ref-1), v[0][2], v[0][3], sum(r[4] for r in v)))
