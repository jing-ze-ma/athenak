"""suite-1011: IDORT 4.2.2 multi-source bubbles, exact free-streaming reference + metrics.
20 opaque discs (IDORT create.py; ana/bub_sources.json) held at T_i in vacuum, periodic unit box, c = 10.
Exact E(x, t) = <I>, I = T_i^4 of the first disc hit along the backward ray if its 3-D path length
s_inplane / sin(theta) <= c t (3-D directions), or s_inplane <= c t (2-D, in-plane rays, IDORT's set).
Inside a disc E = T^4.  usage: python bub.py exact ; python bub.py <rundir> ...
"""
import sys, os, glob, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
C = 10.0
SRC = np.array(json.load(open(os.path.join(HERE, 'bub_sources.json'))))   # x, y, r, T


def first_hit(px, py, phi):
    """in-plane distance to and B of the first disc hit going backward from (px, py) along -u(phi)"""
    ux, uy = np.cos(phi), np.sin(phi)
    best = np.full(phi.shape, np.inf)
    B = np.zeros(phi.shape)
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            for (x, y, r, T) in SRC:
                cx, cy = x + sx - px, y + sy - py
                if cx * cx + cy * cy > (1.6 + r) ** 2:
                    continue
                b = -(cx * ux + cy * uy)          # point P - s u ; |P - s u - C|^2 = s^2 + 2 s b' + |c|^2, b' = c.u
                bb = cx * ux + cy * uy
                cc = cx * cx + cy * cy - r * r
                disc = bb * bb - cc
                ok = disc > 0
                s = np.where(ok, -bb - np.sqrt(np.where(ok, disc, 0)), np.inf)
                # backward ray: P - s u; with s>=0 the hit is at s = -(c.u) - sqrt  (c.u < 0 for discs behind)
                s = np.where(s >= 0, s, np.inf)
                hit = s < best
                best = np.where(hit, s, best)
                B = np.where(hit, T ** 4, B)
    return best, B


def exact_point(px, py, ts, nphi=2048, nmu=16):
    for (x, y, r, T) in SRC:
        dx = (px - x + 0.5) % 1 - 0.5
        dy = (py - y + 0.5) % 1 - 0.5
        if dx * dx + dy * dy < r * r:
            return [T ** 4] * len(ts), [T ** 4] * len(ts)
    phi = (np.arange(nphi) + 0.5) * 2 * np.pi / nphi
    s, B = first_hit(px, py, phi)
    mu, wmu = np.polynomial.legendre.leggauss(nmu)
    st = np.sqrt(1 - mu ** 2)
    e2, e3 = [], []
    for t in ts:
        e2.append(np.mean(np.where(s <= C * t, B, 0.0)))
        e3.append(np.sum(wmu[:, None] / 2 * np.where(s[None, :] <= C * t * st[:, None], B[None, :], 0.0)) / nphi)
    return e2, e3


TS = (0.005, 0.025, 0.1)


def make_exact(n=512, every=8):
    xs = -0.5 + (np.arange(n) + 0.5) / n
    sub = xs[every // 2::every]
    E2 = np.zeros((len(TS), len(sub), len(sub)))
    E3 = np.zeros_like(E2)
    for j, y in enumerate(sub):
        for i, x in enumerate(sub):
            a, b = exact_point(x, y, TS)
            E2[:, j, i] = a
            E3[:, j, i] = b
    np.savez(os.path.join(HERE, 'bub_exact_%d_%d.npz' % (n, every)), x=sub, ts=np.array(TS), E2=E2, E3=E3)
    print('done', E2.max(), E3.max())


def load_at(d, t, var='m1'):
    sys.path.insert(0, '/raven/u/jinma/ATHENAK/wt-littests/vis/python')
    import bin_convert as bc
    fs = sorted(glob.glob(os.path.join(d, 'bin', '*.%s.*.bin' % var)))
    best = None
    for f in fs:
        r = bc.read_binary_as_athdf(f)
        if best is None or abs(r['Time'] - t) < abs(best['Time'] - t):
            best = r
    return best


def metrics(d, every=None, ref='E3'):
    n = int([l.split('=')[1] for l in open(os.path.join(d, 'cmd.txt')).read().split() if l.startswith('mesh/nx1=')][0]) \
        if 'mesh/nx1=' in open(os.path.join(d, 'cmd.txt')).read() else 512
    every = every or (8 if n == 512 else 2)
    ex = np.load(os.path.join(HERE, 'bub_exact_%d_%d.npz' % (n, every)))
    out = {}
    fields = {}
    for k, t in enumerate(TS):
        r = load_at(d, t)
        E = r['m1_e'][0][every // 2::every, every // 2::every]
        Ee = ex[ref][k]
        out['t%.3f' % t] = 'L1=%.3f' % (np.sum(np.abs(E - Ee)) / np.sum(Ee))
        out['tt%.3f' % t] = '%.4f' % r['Time']
        out['Etot%.3f' % t] = '%.3f' % (np.sum(E) / np.sum(Ee))
        out['lit%.3f' % t] = '%.2f/%.2f' % (np.mean(E > 1.0), np.mean(Ee > 1.0))
        fields[t] = (r['m1_e'][0], Ee)
    return out, fields


if __name__ == '__main__':
    if sys.argv[1] == 'exact':
        make_exact()
    else:
        ref = os.environ.get('BUBREF', 'E3')
        for d in sys.argv[1:]:
            o, _ = metrics(d, ref=ref)
            print(os.path.basename(d.rstrip('/')), ' '.join('%s %s' % (k, v) for k, v in o.items()))
