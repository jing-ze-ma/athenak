"""patm / patmT evaluation against the 1-D spherical reference (ref/patm_ref.npz):
per tau band, mean / max |E/J_ref - 1| (cell-volume weighted), F r^2 - 1, the radius
where E = J_ref(tau = 2/3) (photosphere shift), and orders vs N.
usage: patm_eval.py <case> <arm> [levels...]"""
import os
import sys

import numpy as np

import spl

R = np.load(spl.W + '/ref/patm_ref.npz')
rr, JJ = R['br'], R['bJ']


def tau(r):
    return 128.0*(r**-7 - 5.0**-7)


def jref(r):
    return np.exp(np.interp(r, rr, np.log(JJ)))


BANDS = [('tau<0.1', 0.0, 0.1), ('0.1-1', 0.1, 1.0), ('1-10', 1.0, 10.0),
         ('10-20', 10.0, 20.0)]


def ev(case, arm, n):
    d = os.path.join(spl.RUNS, case, '%s_%d' % (arm, n))
    fs = spl.tabs(d)
    t = spl.tab(fs[-1])
    e = t['m1_e']
    f = t['m1_f1']
    poly = case.endswith('_s')
    # the tab x1v is the LOGICAL (unstretched) midpoint: use the sp volume centroid of
    # the stretched cells
    rf = spl.rfaces(1.0, 5.0, len(e), poly)
    r = spl.rcent(rf)
    vol = np.diff(rf**3)
    rel = e/jref(r) - 1.0
    ta = tau(r)
    out = {'n': n, 'st': spl.logstats(d), 'cyc': t['cycle']}
    for nm, a, b in BANDS:
        m = (ta >= a) & (ta < b)
        out[nm] = (np.sum(vol[m]*np.abs(rel[m]))/np.sum(vol[m]), np.max(np.abs(rel[m])))
    out['Fr2'] = np.max(np.abs(f*r*r - 1.0))
    # photosphere: radius where E = J_ref(tau = 2/3) (interpolated, E decreasing)
    rph = rr[np.argmin(np.abs(tau(rr) - 2.0/3.0))]
    jph = jref(rph)
    i = np.where(e < jph)[0][0]
    x = r[i-1] + (r[i] - r[i-1])*(np.log(e[i-1]) - np.log(jph))/(np.log(e[i-1]) - np.log(e[i]))
    out['dph'] = (x - rph)/rph
    out['qtop'] = f[-1]/e[-1]
    out['e_ph'] = np.interp(rph, r, rel)
    out['rel'] = (r, rel)
    # thin envelope (r > 1.2 R_ph): volume-weighted SIGNED mean of E/J_ref - 1
    m = r > 1.2*2.1185
    out['env'] = np.sum(vol[m]*rel[m])/np.sum(vol[m])
    return out


def main():
    case, arm = sys.argv[1:3]
    levs = [int(x) for x in sys.argv[3:]]
    res = [ev(case, arm, n) for n in levs]
    print('%s %s' % (case, arm))
    for o in res:
        st = o['st']
        print('  N %5d cyc %4d Picard %.2f/%g nonconv %g | ' % (
            o['n'], o['cyc'], st.get('pic', -1), st.get('picmax', -1),
            st.get('nonconv', -1)) +
            ' '.join('%s %.2e(%.2e)' % (nm, o[nm][0], o[nm][1]) for nm, _, _ in BANDS) +
            ' | E/J-1 at R_ph %+.2e  dR_ph/R_ph %+.2e  max|Fr2-1| %.1e  q_top %.4f  env %+.2e' % (
                o['e_ph'], o['dph'], o['Fr2'], o['qtop'], o['env']))
    for nm, _, _ in BANDS:
        e = [o[nm][0] for o in res]
        print('  orders %-8s %s' % (nm, spl.fmto(spl.order(e))))


if __name__ == '__main__':
    main()
