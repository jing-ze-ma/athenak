"""ms/cycle (cycles 10-90: a 9-17 s stall between cycles 90 and 100 hits random runs, all arms) per arm, repeats a and b, + solver statistics.
usage: python3 tsum.py runs/<tag> [runs/<tag> ...]"""
import glob
import os
import re
import sys

for D in sys.argv[1:]:
    res = {}
    order = []
    for d in sorted(glob.glob(D + '/*_[ab]')):
        arm, rep = os.path.basename(d).rsplit('_', 1)
        if arm not in res:
            res[arm] = {}
            order.append(arm)
        try:
            L = open(d + '/run.log').read()
        except OSError:
            continue
        el = {int(c): float(e) for e, c in re.findall(r'elapsed=(\S+) cycle=(\d+)', L)}
        g = lambda p: (float(re.search(p, L).group(1)) if re.search(p, L) else float('nan'))
        if 10 not in el or 90 not in el:
            res[arm][rep] = None
            continue
        res[arm][rep] = (1e3 * (el[90] - el[10]) / 80,
                         g(r'Picard iterations mean=(\S+)'), g(r'NON-CONVERGED=(\S+)'),
                         g(r'inner iterations mean=(\S+)'), g(r'breakdowns=(\S+)'),
                         g(r'fallbacks=(\S+)'))
    print(D)
    print('  %-9s %-17s %-7s %-7s %-8s %-6s %s' % ('arm', 'ms/cycle a, b', 'Picard',
                                                  'nonconv', 'inner', 'break', 'fallb'))
    for a in order:
        v = [res[a].get(r) for r in 'ab']
        cy = ', '.join('%6.1f' % x[0] if x else '  fail' for x in v)
        ok = [x for x in v if x]
        if ok:
            x = ok[0]
            print('  %-9s %-17s %-7.3f %-7g %-8.2f %-6g %g' % (a, cy, x[1], x[2], x[3],
                                                            x[4], x[5]))
        else:
            print('  %-9s %-17s' % (a, cy))
