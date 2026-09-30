"""ana.py: velocities from trans.py spectra the way the observers measure them.

  python ana.py <run>_<tag> [R=70000]      -> printed table (also out/<run>_<tag>.txt)

Fe: CCF(v) = sum_w m_w * excess absorption of window w at v (weighted binary mask; m_w =
windless/rotationless template line depth, out/<run>_template.npz), after convolution with a
Gaussian LSF of resolving power R (ESPRESSO 4-UT MR ~ 70000).  Centre = Gaussian + constant
fit over +-30 km/s; also the first moment over +-25 km/s.  Na D: the same per line (D2, D1).
RV sign: positive = redshift.  Contribution functions: 16/50/84 % points in log p (bar).
"""
import sys
import numpy as np
from scipy.optimize import curve_fit

W = '/viper/ptmp2/jinma/w121prod_0929/synth2_trans/out/'


def lsf(y, u, R):
    if R <= 0:
        return y
    sig = 2.99792458e10/R/2.3548
    du = u[1] - u[0]
    k = np.arange(-int(4*sig/du), int(4*sig/du) + 1)*du
    g = np.exp(-0.5*(k/sig)**2)
    g /= g.sum()
    return np.apply_along_axis(lambda a: np.convolve(a, g, 'same'), -1, y)


def gfit(u, y, lim=30e5):
    m = np.abs(u) <= lim
    uu, yy = u[m]/1e5, y[m]
    i = np.argmax(yy)
    p0 = [yy[i] - yy.min(), uu[i], 3.0, yy.min()]
    try:
        p, _ = curve_fit(lambda x, A, x0, s, c0: A*np.exp(-0.5*((x - x0)/s)**2) + c0,
                         uu, yy, p0=p0, maxfev=20000)
        return p[1], abs(p[2]), p[0]
    except RuntimeError:
        return np.nan, np.nan, np.nan


def moment(u, y, lim=25e5):
    m = np.abs(u) <= lim
    return (u[m]*y[m]).sum()/y[m].sum()/1e5


def cfq(cf, lpb):
    c = np.cumsum(cf)
    c /= c[-1]
    return [np.interp(q, c, lpb + 0.05) for q in (0.16, 0.5, 0.84)]


def analyse(name, R=70000.0, ret=False):
    z = np.load(W + name + '.npz')
    run = name.split('_')[0]
    mw = np.load(W + run + '_template.npz')['mw']
    al = z['alpha']
    UW, UNA = z['UW'], z['UNA']
    Efe = z['Afe'] - z['Afe_c'][..., None]          # (ep, limb, w, u)
    Ena = z['Ana'] - z['Ana_c'][..., None]          # (ep, limb, q, u)
    ccf = lsf(np.einsum('w,elwu->elu', mw, Efe), UW, R)
    na = lsf(Ena, UNA, R)
    res = {}
    sets = {'all': al == al, 'first half': al < 0, 'second half': al > 0}
    for sn, sm in sets.items():
        if not sm.any():
            continue
        for ln, lsl in (('morning', [0]), ('evening', [1]), ('both', [0, 1])):
            y = ccf[sm][:, lsl].sum((0, 1))/sm.sum()
            g = gfit(UW, y)
            yn = na[sm][:, lsl].sum((0, 1))/sm.sum()   # (q, u), per epoch
            gn = [gfit(UNA, yn[q]) for q in range(2)]
            res[(sn, ln)] = dict(fe=g[0], fe_w=g[1], fe_m=moment(UW, y), fe_amp=g[2]/mw.sum(),
                                 na=[x[0] for x in gn], na_w=[x[1] for x in gn],
                                 na_m=[moment(UNA, yn[q]) for q in range(2)],
                                 na_amp=[x[2] for x in gn])
    lpb = z['lpb']
    cf = {}
    for ln, lsl in (('morning', [0]), ('evening', [1]), ('both', [0, 1])):
        for sp in ('fe', 'na'):
            c = z['CF' + sp][:, lsl].sum((0, 1))
            top = z['CF%s_top' % sp][:, lsl].sum()
            win = (-5, -3) if sp == 'fe' else (-6, -4)
            fin = c[(lpb >= win[0]) & (lpb < win[1])].sum()/c.sum()
            cf[(sp, ln)] = cfq(c, lpb) + [top/c.sum(), fin]
    lines = ['# %s  opt %s  LSF R=%g' % (name, z['opt'], R),
             '# set          limb      Fe CCF RV [km/s] gauss / moment  (sigma, amp ppm)'
             '   Na D2 / D1 RV gauss   (moment)   Na amp ppm']
    for k, v in res.items():
        lines.append('%-12s %-8s Fe %+6.2f / %+6.2f (%4.1f, %6.1f)   Na %+6.2f %+6.2f  '
                     '(%+6.2f %+6.2f)  %6.0f %6.0f' % (
                         k[0], k[1], v['fe'], v['fe_m'], v['fe_w'], v['fe_amp']*1e6,
                         v['na'][0], v['na'][1], v['na_m'][0], v['na_m'][1],
                         v['na_amp'][0]*1e6, v['na_amp'][1]*1e6))
    lines.append('# contribution (marginal line absorption) log10 p [bar]: 16/50/84 %,'
                 ' fraction above x1max, fraction in the old window')
    for k, v in cf.items():
        lines.append('CF %-2s %-8s %6.2f %6.2f %6.2f   top %.3f   oldwin %.3f' % (k + tuple(v)))
    txt = '\n'.join(lines)
    if ret:
        return res, cf, txt
    open(W + name + '.txt', 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    R = 70000.0
    for a in sys.argv[2:]:
        if a.startswith('R='):
            R = float(a[2:])
    analyse(sys.argv[1], R)
