"""vimpbc: cmp.py runF runT  -> per dump: per variable max rel diff and the radial index (from top) where it sits;
top-cell velx sign; hst L_top/L_top0, Mdot_top, Picard; run.log Picard summary, s/cycle."""
import sys, glob, re, os, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import bin_convert as bc
a, b = sys.argv[1], sys.argv[2]
def summ(d):
    log = open(d + '/run.log').read()
    m = re.search(r'Picard iterations mean=(\S+) max=(\S+) NON-CONVERGED=(\S+)', log)
    el = re.findall(r'elapsed=(\S+) cycle=(\d+)', log)
    ic = re.findall(r'cpu time used\s*=\s*(\S+)', log)
    zc = re.search(r'zone-cycles/cpu_second\s*=\s*(\S+)', log)
    rc = open(d + '/rc.txt').read().strip() if os.path.exists(d + '/rc.txt') else '?'
    s = '%s %s Picard %s/%s NONCONV %s' % (os.path.basename(d), rc, m.group(1) if m else '-', m.group(2) if m else '-', m.group(3) if m else '-')
    if len(el) > 2:
        e0, c0 = float(el[1][0]), int(el[1][1]); e1, c1 = float(el[-1][0]), int(el[-1][1])
        s += '  s/cycle %.3f (c%d-%d)' % ((e1 - e0)/max(c1 - c0, 1), c0, c1)
    if zc: s += ' zc/s ' + zc.group(1)
    print(s)
summ(a); summ(b)
def hst(d):
    f = d + '/agcar3d.user.hst'; x = np.loadtxt(f); hdr = [l for l in open(f) if l.startswith('#  [')][0]
    return {n: x[:, i] for i, n in enumerate(re.findall(r'\[\d+\]=(\S+)', hdr))}
Ha, Hb = hst(a), hst(b); n = min(len(Ha['time']), len(Hb['time']))
print('hst rows', n)
for k in ('L_top', 'Mdot_top', 'Min_top', 'Etot', 'M_tot', 'Picard'):
    if k not in Ha: continue
    xa, xb = Ha[k][:n], Hb[k][:n]
    idx = sorted(set([1, 10, 30, 60, 100, 200, n - 1]) & set(range(n)))
    sc = np.abs(xa).max() or 1.0
    print('  %-9s F: %s' % (k, ' '.join('%.5g' % xa[i] for i in idx)))
    print('  %-9s T: %s   max|dT-F|/max|F| %.3e' % (k, ' '.join('%.5g' % xb[i] for i in idx), np.abs(xb - xa).max()/sc))
for kind in ('hydro_w', 'm1'):
    fa = sorted(glob.glob(a + '/bin/agcar3d.%s.*.bin' % kind))
    for f in fa:
        g = f.replace(a, b, 1)
        if not os.path.exists(g): continue
        ha, hb = bc.read_binary_as_athdf(f), bc.read_binary_as_athdf(g)
        c = int(ha['NumCycles']); out = []
        for v in [k for k in ha if k not in ('x1f','x1v','x2f','x2v','x3f','x3v','Time','NumCycles','MaxLevel')]:
            x, y = np.asarray(ha[v], float), np.asarray(hb[v], float)
            dd = np.abs(y - x); sc = np.abs(x).max() or 1.0
            if dd.max() == 0: out.append('%s 0' % v); continue
            rel = dd/(np.abs(x) + 1e-30*sc + 1e-12*sc)  # relative to local, floored at 1e-12 of max
            k3, j3, i3 = np.unravel_index(dd.argmax(), dd.shape)
            nx = dd.shape[2]
            # outermost radial index where diff > 1e-3 of its max
            rows = np.where(dd.max(axis=(0, 1)) > 1e-3*dd.max())[0]
            out.append('%s %.2e(i%d/top-%d; span i%d-%d)' % (v, dd.max()/sc, i3, nx-1-i3, rows.min(), rows.max()))
        extra = ''
        if kind == 'hydro_w':
            vx = np.asarray(ha['velx'], float)[:, :, -1]; vxb = np.asarray(hb['velx'], float)[:, :, -1]
            extra = ' | top velx F [%.3g,%.3g] T [%.3g,%.3g] km/s' % (vx.min()/1e5, vx.max()/1e5, vxb.min()/1e5, vxb.max()/1e5)
        print('c%-4d %s: %s%s' % (c, kind, '  '.join(out), extra))
