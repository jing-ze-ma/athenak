import sys, re, numpy as np
Fin = 2.475202e15
P = {'pm': r'Picard iterations mean=(\S+)', 'px': r'Picard iterations mean=\S+ max=(\S+)',
     'nc': r'NON-CONVERGED=(\S+)', 'pf': r'positivity fallbacks=(\S+)',
     'me': r'min E from the solve=(\S+)'}
for d in sys.argv[1:]:
    log = open(d + '/run.log').read()
    el = [(float(a), int(b)) for a, b in re.findall(r'elapsed=(\S+) cycle=(\d+)', log)]
    (e0, c0), (e1, c1) = el[1], el[-1]
    v = {k: (re.search(p, log).group(1) if re.search(p, log) else '-') for k, p in P.items()}
    u = np.loadtxt(d + '/m1slab.user.hst')
    ms = 1e3 * (e1 - e0) / (c1 - c0)
    print('%-5s cyc %d t=%.1f ms/cyc(c%d-%d) %.1f Picard %s/%s nonconv %s posfb %s minE %s '
          'F1top/Fin %.4f [%.4f/%.4f] dt %.5g [min %.5g]'
          % (d, c1, u[-1, 0], c0, c1, ms, v['pm'], v['px'], v['nc'], v['pf'], v['me'],
             u[-1, 2] / Fin, u[:, 2].min() / Fin, u[:, 2].max() / Fin, u[-1, 1], u[:, 1].min()))
    for l in log.split('\n'):
        if 'fallback' in l.lower() or 'BiCGStab' in l or 'Newton' in l:
            print('      ' + l[:170])
