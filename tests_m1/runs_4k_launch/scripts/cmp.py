# usage: python3 cmp.py <dir> ref arm...   (hst of runs/<ref> vs runs/<arm>, or cpu/...)
# prints: bitwise identical?, max |d|/colmax over all rows (hydro, user), last-row relative
# diff of totE (hydro col 6), KE1 (7), F1top (user col 2), and the common last time.
import sys, os, filecmp, numpy as np
D = sys.argv[1]
def ld(r, f):
    return np.loadtxt(os.path.join(D, r, 'm1slab.%s.hst' % f))
ref = sys.argv[2]
for a in sys.argv[3:]:
    out = [a]
    same = all(filecmp.cmp(os.path.join(D, ref, 'm1slab.%s.hst' % f),
                           os.path.join(D, a, 'm1slab.%s.hst' % f), shallow=False)
               for f in ('hydro', 'user'))
    out.append('BITWISE' if same else 'differs')
    for f in ('hydro', 'user'):
        x, y = ld(ref, f), ld(a, f)
        n = min(len(x), len(y))
        x, y = x[:n], y[:n]
        cm = np.maximum(np.abs(x).max(axis=0), 1e-300)
        out.append('%s all %.1e' % (f, (np.abs(x - y)/cm).max()))
    h0, h1 = ld(ref, 'hydro'), ld(a, 'hydro')
    u0, u1 = ld(ref, 'user'), ld(a, 'user')
    n = min(len(h0), len(h1)); m = min(len(u0), len(u1))
    rel = lambda p, q: abs(p - q)/max(abs(p), 1e-300)
    out.append('t=%.2f totE %.1e KE1 %.1e KE2 %.1e F1top %.1e' % (
        h0[n-1, 0], rel(h0[n-1, 6], h1[n-1, 6]), rel(h0[n-1, 7], h1[n-1, 7]),
        rel(h0[n-1, 8], h1[n-1, 8]), rel(u0[m-1, 2], u1[m-1, 2])))
    print('  '.join(out))
