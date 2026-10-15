#!/usr/bin/env python3
"""ridort-1014 summ.py <reval.txt> ...: one line per test family (<kind>/<case>/<arm prefix>): n dirs, n bitwise (maxrel 0),
worst E-type and worst resolved-variable ratio D_R/D_lad (with the level), D_R trend over levels (finest/coarsest), Picard REF -> R.
Verdict: PASS if every ratio < 1 (R's change below the family's own ladder step), INERT if all dirs are bitwise, CHECK otherwise."""
import re
import sys

fam = {}
for fn in sys.argv[1:]:
    for ln in open(fn):
        w = [x.strip() for x in ln.split('|')]
        if len(w) < 6:
            if 'FAIL' in ln:
                print('FAIL', ln.strip())
            continue
        d = w[0]
        m = re.match(r'^(.*?)(\d+)$', d)
        key = m.group(1) if m else d
        lev = int(w[1].split('=')[1])
        mr = float(w[2].split()[1])

        def parse(s):
            t = s.split()
            if len(t) < 5 or t[4] == '-':
                return None
            return (float(t[4]), t[1], float(t[2]), lev)
        e = parse(w[3])
        r = parse(w[4])
        pr = re.findall(r'Picard (\S+)/(\S+)', ln)
        f = fam.setdefault(key, {'n': 0, 'bit': 0, 'e': None, 'r': None, 'p0': [], 'p1': [], 'dr': []})
        f['n'] += 1
        f['bit'] += (mr == 0.0)
        for k, v in (('e', e), ('r', r)):
            if v is not None and (f[k] is None or v[0] > f[k][0]):
                f[k] = v
        if e is not None or w[3].split()[1:2]:
            try:
                f['dr'].append((lev, float(w[3].split()[2])))
            except (IndexError, ValueError):
                pass
        if pr:
            f['p0'].append(float(pr[0][0]))
            f['p1'].append(float(pr[0][1]))
for k in sorted(fam):
    f = fam[k]
    worst = max([x[0] for x in (f['e'], f['r']) if x is not None], default=None)
    if f['bit'] == f['n']:
        v = 'INERT(bitwise)'
    elif worst is None:
        v = 'NO-LADDER'
    else:
        v = 'PASS' if worst < 1.0 else 'CHECK'
    dr = sorted(f['dr'])
    tr = ('%.1e->%.1e' % (dr[0][1], dr[-1][1])) if len(dr) > 1 else ('%.1e' % dr[0][1] if dr else '-')
    fe = '-' if f['e'] is None else '%s %.1e (l%d)' % (f['e'][1].split(':')[-1], f['e'][0], f['e'][3])
    fr = '-' if f['r'] is None else '%s %.1e (l%d)' % (f['r'][1].split(':')[-1], f['r'][0], f['r'][3])
    p0 = sum(f['p0']) / len(f['p0']) if f['p0'] else float('nan')
    p1 = sum(f['p1']) / len(f['p1']) if f['p1'] else float('nan')
    print('%-40s n=%-2d bit=%-2d E-ratio %-26s res-ratio %-26s D_R(E) %-17s Picard %.2f->%.2f  %s' % (
        k, f['n'], f['bit'], fe, fr, tr, p0, p1, v))
