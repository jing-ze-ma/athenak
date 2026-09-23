"""SC ms/call and ms/cycle (cycles 10-100) per arm, repeats a and b.
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
        sc = re.search(r'SC seconds=\S+ \((\S+) per call\)', L)
        pic = re.search(r'Picard iterations mean=(\S+)', L)
        nc = re.search(r'NON-CONVERGED=(\S+)', L)
        if 10 not in el or 100 not in el:
            res[arm][rep] = None
            continue
        res[arm][rep] = (1e3 * float(sc.group(1)) if sc else float('nan'),
                         1e3 * (el[100] - el[10]) / 90,
                         float(pic.group(1)) if pic else float('nan'),
                         float(nc.group(1)) if nc else float('nan'))
    print(D)
    print('  %-9s %-17s %-17s %-7s %s' % ('arm', 'SC ms/call a, b', 'ms/cycle a, b',
                                          'Picard', 'nonconv'))
    for a in order:
        v = [res[a].get(r) for r in 'ab']
        sc = ', '.join('%6.2f' % x[0] if x else '  fail' for x in v)
        cy = ', '.join('%6.1f' % x[1] if x else '  fail' for x in v)
        ok = [x for x in v if x]
        print('  %-9s %-17s %-17s %-7s %s' % (a, sc, cy,
                                              '%.3f' % ok[0][2] if ok else '-',
                                              '%g' % ok[0][3] if ok else '-'))
