"""suite-1011: radiation shocks against the IDORT/J21 transport reference (Rad_shock_sol in the IDORT Zenodo
package 15032138; Ferguson+17-type S_N solution, 4096 cells, nondimensional: rho0 = T0 = 1, p = rho T/gamma,
E = (T_r/T0)^4).  Our runs are in the cgs scaling of t7_radshock.py (rho0 = 5.69, T0 = 2.18e6 K, L = 1.000337 cm,
tunit = 8.2499208e7, arad_code = 1.6332399e-46).  Profiles are aligned on the density jump (max |d rho/dx|).
usage: python shockcmp.py <M: 1.2|2|3> <rundir> ...
"""
import sys, os, glob
import numpy as np
REF = '/raven/ptmp/jinma/suite_1011/ext/idort_mini/x/IDORT_mini/shock_1d/Rad_shock_sol/M%s'
RHO0, T0, L, TUN, AR, GAM = 5.69, 2.18e6, 1.000337, 8.2499208e7, 1.6332399e-46, 5.0 / 3.0


def ref(m):
    lab = {'1.2': '1p2', '2': '2', '3': '3'}[m]
    r = {k: np.loadtxt(os.path.join(REF % lab, '%s_4096.dat' % k)) for k in ('x', 'rho', 'p', 'E', 'u')}
    r['tg'] = GAM * r['p'] / r['rho']
    r['tr'] = np.maximum(r['E'], 0) ** 0.25
    return r


def rdtab(f):
    hdr = [c for c in open(f).readlines()[1].lstrip('#').split() if c not in ('k', 'x3v', 'j', 'x2v')]
    a = np.loadtxt(f)
    return {h: a[:, i] for i, h in enumerate(hdr)}


def run(d, last=True):
    fs = sorted(glob.glob(os.path.join(d, 'tab', '*.m1.*.tab')))
    f = fs[-1] if last else fs[0]
    m = rdtab(f); h = rdtab(f.replace('.m1.', '.hydro_w.'))
    t = float(open(f).readline().split('time=')[1].split()[0])
    j = h['j'] if 'j' in h else None
    # average over the x2 cells
    x = np.unique(h['x1v'])
    def avg(v):
        return np.array([v[h['x1v'] == xx].mean() for xx in x])
    rho = avg(h['dens']) / RHO0
    tg = (GAM - 1) * avg(h['eint']) / (avg(h['dens']) * TUN) / T0
    tr = (np.maximum(avg(m['m1_e']), 0) / AR) ** 0.25 / TUN / T0
    return t, x / L, rho, tg, tr


def jump(x, rho):
    g = np.abs(np.diff(rho))
    i = np.argmax(g)
    return 0.5 * (x[i] + x[i + 1])


def precursor(x, tg, xs):
    """length ahead of the shock where T_gas exceeds T0 + 0.01 (T(shock-) - T0)"""
    up = x < xs
    t1 = tg[up][-1]
    lev = 1 + 0.01 * (t1 - 1)
    xi = x[up][np.argmax(tg[up] > lev)]
    return xs - xi


def compare(m, d):
    R = ref(m)
    xsr = jump(R['x'], R['rho'])
    out = {}
    for tag, last in (('t0', False), ('end', True)):
        t, x, rho, tg, tr = run(d, last)
        xs = jump(x, rho)
        xr = R['x'] - xsr
        xx = x - xs
        sel = (xx > xr.min()) & (xx < xr.max()) & (np.abs(xx) > 3 * (x[1] - x[0]))   # skip the jump cells
        e = []
        for v, k in ((rho, 'rho'), (tg, 'tg'), (tr, 'tr')):
            rv = np.interp(xx[sel], xr, R[k])
            e.append(np.sum(np.abs(v[sel] - rv)) / np.sum(np.abs(rv)))
        out[tag] = (t, e, precursor(x, tg, xs), tg.max())
    pr = precursor(R['x'], R['tg'], xsr)
    return out, pr, R['tg'].max()


if __name__ == '__main__':
    m = sys.argv[1]
    for d in sys.argv[2:]:
        o, pr, spk = compare(m, d)
        s = os.path.basename(os.path.dirname(d.rstrip('/'))) + '/' + os.path.basename(d.rstrip('/'))
        for tag in ('t0', 'end'):
            t, e, p, tmax = o[tag]
            print('%-20s %-3s t=%.3e L1 rho %.4f Tg %.4f Tr %.4f | precursor %.4f (ref %.4f, %+.1f%%) | Tg max %.4f (ref %.4f)'
                  % (s, tag, t, e[0], e[1], e[2], p, pr, 100 * (p / pr - 1), tmax, spk))
