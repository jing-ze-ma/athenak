#!/usr/bin/env python3
"""m1-positivity decisive arms: per arm the time reached, NaN/abort signs and the rad_m1
end-of-run counters (run.log of runs/<arm>).  usage: analyze.py [arm ...]"""
import glob, os, re, sys
W = '/viper/ptmp2/jinma/m1pos_0930/runs'
arms = sys.argv[1:] or sorted(os.path.basename(d) for d in glob.glob(W + '/*') if os.path.isdir(d))
keys = ['Picard iterations', 'NON-CONVERGED', 'm1-positivity', 'vimp positivity', 'floor clips',
        'time_scheme=hesdirk2:', 'guard=', 'bicgstab: breakdowns']
for a in arms:
    f = os.path.join(W, a, 'run.log')
    if not os.path.exists(f):
        print(f'== {a}: no run.log'); continue
    tl, nan, term, summ = None, 0, [], {}
    nadm = 0
    with open(f, errors='replace') as fh:
        for ln in fh:
            m = re.match(r'cycle=(\d+) time=([0-9.eE+-]+) dt=([0-9.eE+-]+)', ln)
            if m: tl = (int(m.group(1)), float(m.group(2)), float(m.group(3)))
            m = re.search(r'elapsed=\S+ cycle=(\d+) time=([0-9.eE+-]+) dt=([0-9.eE+-]+)', ln)
            if m: tl = (int(m.group(1)), float(m.group(2)), float(m.group(3)))
            if re.search(r'\bnan\b|\binf\b', ln, re.I) and 'max=inf' not in ln: nan += 1
            if 'NOT ADMISSIBLE' in ln: nadm += 1
            if 'Terminating' in ln or 'FATAL' in ln or '###' in ln: term.append(ln.strip())
            for k in keys:
                if k in ln and ln.startswith('<rad_m1>'): summ[k] = ln.strip()[:260]
    print(f'== {a}: last (cycle, time, dt) = {tl}; lines with nan/inf = {nan}; '
          f'NOT ADMISSIBLE = {nadm}')
    for t in term[-3:]: print('   ', t[:200])
    for k in keys:
        if k in summ: print('   ', summ[k])
