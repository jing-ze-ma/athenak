"""littests-1011 stage 2 analysis: eq | nlte | sdif | ddif   usage: s2.py <kind> <rundir>..."""
import glob
import math
import os
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/blend_1022/wt/tests_m1')
import common as cm  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lcommon import runinfo  # noqa: E402
from cb import readin  # noqa: E402


def files(d, var, ext):
    return sorted(glob.glob('%s/**/*.%s.*.%s' % (d, var, ext), recursive=True)) or \
        sorted(glob.glob('%s/**/*.out1.*.%s' % (d, ext), recursive=True))


def prof(dump, name):
    a = dump.data[name]
    return a.reshape(-1, a.shape[-1]).mean(axis=0)


def nlte(d):
    p = readin(d + '/in.athinput')
    eps = float(p['rad_m1/kappa_p'])
    f = sorted(glob.glob(d + '/**/*.tab', recursive=True))[-1]
    du = cm.load_dump(f)
    x = du.x1
    E = prof(du, 'm1_e') if 'm1_e' in du.data else prof(du, 'E')
    tau = 1.0e-3*(np.exp(10.0 - x) - 1.0)
    ex = 1.0 - np.exp(-np.sqrt(3*eps)*tau)/(1.0 + np.sqrt(eps))
    s = (tau > 1e-3) & (tau < 1e5)
    rel = np.abs(E - ex)/ex
    return dict(eps=eps, x=x, tau=tau, E=E, ex=ex, L1=float(np.abs(E-ex)[s].sum()/ex[s].sum()),
                max=float(rel[s].max()), top=float(E[-1]/ex[-1]), t=du.time)


def gauss(x, t, D, v=0.0, L=None, bg=math.exp(-10.0)):
    s = 1.0 + 160.0*D*t
    xx = x - v*t
    if L is not None:
        xx = (xx + 0.5*L) % L - 0.5*L
    return np.exp(-40.0*xx*xx/s)/np.sqrt(s) + bg


def dif(d, kind):
    p = readin(d + '/in.athinput')
    c = float(p['rad_m1/c_light'])
    D = c/(3.0*float(p['rad_m1/kappa_s']))
    v = float(p.get('problem/gas_vx', 0.0))
    L = 20.0 if kind == 'ddif' else None
    out = []
    for f in sorted(glob.glob(d + '/**/*.tab', recursive=True)):
        du = cm.load_dump(f)
        if du.time <= 0:
            continue
        E = prof(du, 'm1_e')
        ex = gauss(du.x1, du.time, D, v, L)
        out.append(dict(t=du.time, x=du.x1, E=E, ex=ex,
                        L1=float(np.abs(E - ex).sum()/np.abs(ex - math.exp(-10)).sum()),
                        peak=float(E.max()/ex.max()), xc=float((du.x1*(E - E.min())).sum()/(E - E.min()).sum())))
    return out


def eq(d):
    p = readin(d + '/in.athinput')
    c = float(p['rad_m1/c_light'])
    ka = float(p['rad_m1/kappa_p'])
    out = []
    fm = files(d, 'm1', 'bin')
    fh = files(d, 'hydro_w', 'bin')
    for a, b in zip(fm, fh):
        m = cm.load_dump(a)
        h = cm.load_dump(b)
        E = float(np.mean(m.data['m1_e']))
        dn = float(np.mean(h.data['dens']))
        T = float(np.mean(h.data['eint']))*(2.0/3.0)/dn if 'eint' in h.data else \
            float(np.mean(h.data['press']))/dn
        vx = float(np.mean(h.data['velx']))
        out.append(dict(t=m.time, E=E, T=T, v=vx))
    return dict(c=c, ka=ka, hist=out)


def eq_ref(T0, E0, ka, c, ts, rho=1.0, gm1=2.0/3.0):
    """exact static relaxation: de/dt = -c rho ka (T^4 - E), e = rho T/gm1, e + E const"""
    from scipy.integrate import solve_ivp
    tot = rho*T0/gm1 + E0

    def f(t, y):
        T = y[0]
        E = tot - rho*T/gm1
        return [-c*rho*ka*(T**4 - E)*gm1/rho]
    s = solve_ivp(f, (0, max(ts)), [T0], t_eval=ts, rtol=1e-11, atol=1e-12, method='Radau')
    T = s.y[0]
    return T, tot - rho*T/gm1


if __name__ == '__main__':
    k = sys.argv[1]
    for d in sys.argv[2:]:
        r = runinfo(d)
        tag = d.split('runs/')[-1]
        if k == 'nlte':
            m = nlte(d)
            print('%-24s eps %.0e L1 %.4f maxrel %.4f E_top/ex %.4f t %.3g picard %.1f' % (
                tag, m['eps'], m['L1'], m['max'], m['top'], m['t'], r.get('picard', -1)))
        elif k in ('sdif', 'ddif'):
            for m in dif(d, k):
                print('%-24s t %7.2f L1 %.4f peak/ex %.4f xc %.4f' % (tag, m['t'], m['L1'], m['peak'], m['xc']))
        elif k == 'eq':
            m = eq(d)
            h = m['hist']
            print('%-24s final t %.3g E %.5f T %.5f v %.5f' % (tag, h[-1]['t'], h[-1]['E'], h[-1]['T'], h[-1]['v']))


def eq_table(d):
    p = readin(d + '/in.athinput')
    m = eq(d)
    T0 = float(p['problem/gas_eint'])*(2.0/3.0)
    E0 = float(p['problem/e_rad'])
    h = m['hist']
    ts = np.array([q['t'] for q in h])
    Ec = np.array([q['E'] for q in h])
    Tc = np.array([q['T'] for q in h])
    Tr, Er = eq_ref(T0, E0, m['ka'], m['c'], ts)
    s = ts > 0
    return dict(t=ts, E=Ec, T=Tc, Er=Er, Tr=Tr, v=h[-1]['v'],
                errE=float(np.max(np.abs(Ec[s]/Er[s] - 1))) if s.any() else np.nan,
                errT=float(np.max(np.abs(Tc[s]/Tr[s] - 1))) if s.any() else np.nan)
