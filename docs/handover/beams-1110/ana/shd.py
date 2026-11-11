"""littests-1011 static shadow (Ma+2025 4.1.4 / J21 4.3 / JSD12 6.5).
Exact (beam part): E_b sum over beams of exp(-tau) along the straight ray back to x1min,
sigma = rho^2 T^-3.5 = rho^5.5 analytic (T = 1/rho), periodic x2, 4 substeps per cell.
The thermal part (S = a T^4 <= 1 in the background, vs E_b = 200 per beam) is NOT in
the exact field; it matters only in the umbra.
usage: shd.py <rundir> ...  ->  metrics; cache exact in ana/shd_exact_<nx1>.npz"""
import os
import sys
import numpy as np
from lcommon import load, runinfo, W
from cb import readin, beams


def rho(p, x, y):
    ax, ay = float(p['problem/shadow_ax']), float(p['problem/shadow_ay'])
    dl, r0, rt = (float(p['problem/shadow_delta']), float(p['problem/shadow_rho0']),
                  float(p['problem/shadow_ratio']))
    d2 = (x/ax)**2 + (y/ay)**2
    ex = np.clip(dl*(d2 - 1.0), -700, 700)
    return r0*(1.0 + (rt - 1.0)/(1.0 + np.exp(ex)))


def sigma(p, x, y):
    d = rho(p, x, y)
    return d**5.5   # rho kappa, kappa = rho T^-3.5, T = 1/rho


def exact(p, x1, x2, each=False):
    nx1 = len(x1)
    fn = '%s/ana/shd_exact_%d.npz' % (W, nx1)
    if os.path.exists(fn):
        z = np.load(fn)
        return [z['E0'], z['E1']] if each else z['E0'] + z['E1']
    x1a = float(p['mesh/x1min'])
    y2a, L = float(p['mesh/x2min']), float(p['mesh/x2max']) - float(p['mesh/x2min'])
    dx1 = x1[1] - x1[0]
    res = []
    for b in beams(p):
        tn = np.tan(b['ang'])
        E = np.zeros((len(x2), nx1))
        for i, xc in enumerate(x1):
            n = 8*(i + 1)
            s = x1a + (np.arange(n) + 0.5)*(xc - x1a)/n
            ds = (xc - x1a)/n/np.cos(b['ang'])
            yy = x2[:, None] - (xc - s[None, :])*tn
            yy = (yy - y2a) % L + y2a
            tau = sigma(p, s[None, :], yy).sum(axis=1)*ds
            E[:, i] = b['e']*np.exp(-tau)
        res.append(E)
    np.savez(fn, E0=res[0], E1=res[1])
    return res if each else sum(res)


def geom(ang_deg=15.0, a=0.1, b=0.06):
    """umbra tip and the x2 of the lit edges at x for the ellipse a, b (midpoint edge)"""
    t = np.tan(np.radians(ang_deg))
    c = np.sqrt(a*a*t*t + b*b)   # tangent line y = t x + c
    return t, c


def edges(x2, v, lo=True):
    """x2 positions where v crosses half of its max on the x2 > 0 side (sorted)"""
    s = x2 > 0
    y, w = x2[s], v[s]
    h = 0.5*w.max()
    out = []
    for k in range(len(y) - 1):
        if (w[k] - h)*(w[k+1] - h) < 0:
            out.append(y[k] + (h - w[k])*(y[k+1] - y[k])/(w[k+1] - w[k]))
    return out


def metrics(d, cuts=(0.25, 0.4)):
    p = readin(d + '/in.athinput')
    v = load(d)
    x1, x2 = v['x1'], v['x2']
    Eex = exact(p, x1, x2)
    out = {'run': d.split('runs/')[-1]}
    out.update(runinfo(d))
    lit = np.abs(x1 - (-0.5)) < 1  # all
    for nm in ('E', 'J'):
        a = v.get(nm)
        if a is None:
            continue
        out[nm + '_L1'] = np.abs(a - Eex).sum()/Eex.sum()
        for xc in cuts:
            i = int(np.argmin(np.abs(x1 - xc)))
            out[nm + '_edges_%.2f' % xc] = edges(x2, a[:, i])
            out[nm + '_umbra_%.2f' % xc] = a[np.argmin(np.abs(x2)), i]
            out[nm + '_wing_%.2f' % xc] = a[:, i].max()
    for xc in cuts:
        i = int(np.argmin(np.abs(x1 - xc)))
        out['ex_edges_%.2f' % xc] = edges(x2, Eex[:, i])
        out['ex_umbra_%.2f' % xc] = Eex[np.argmin(np.abs(x2)), i]
        out['ex_wing_%.2f' % xc] = Eex[:, i].max()
    return out


if __name__ == '__main__':
    for d in sys.argv[1:]:
        m = metrics(d)
        f = lambda a: ' '.join('%.4f' % q for q in a)  # noqa: E731
        s = '%-22s L1 E %.3f J %.3f' % (m['run'], m['E_L1'], m['J_L1'])
        for xc in (0.25, 0.4):
            k = '%.2f' % xc
            s += ' | x=%s edges E [%s] J [%s] ex [%s] umbra E %.1f J %.1f ex %.1f' % (
                k, f(m['E_edges_' + k]), f(m['J_edges_' + k]), f(m['ex_edges_' + k]),
                m['E_umbra_' + k], m['J_umbra_' + k], m['ex_umbra_' + k])
        s += ' | picard %.1f nc %d %.0fs' % (m.get('picard', -1), m.get('nonconv', -1),
                                           m.get('sec', -1))
        print(s)
