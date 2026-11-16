"""moveq-1011: time series of the uniform moving-equilibrium runs (tab outputs).
usage: eqts.py <rundir>...   prints t, E, F1, v, T, total momentum and energy."""
import glob
import sys
import numpy as np

C = 10.0
GM1 = 2.0/3.0


def rd(f):
    with open(f) as fh:
        l1 = fh.readline()
        t = float(l1.split('time=')[1].split()[0])
        hdr = [h for h in fh.readline().lstrip('#').split() if h not in ('k', 'x3v', 'j', 'x2v')]
    a = np.loadtxt(f, ndmin=2)
    return t, {h: a[:, n].mean() for n, h in enumerate(hdr)}


def series(d):
    out = []
    for fm in sorted(glob.glob(d + '/tab/*.m1.*.tab')):
        fh = fm.replace('.m1.', '.hydro_w.')
        t, m = rd(fm)
        _, h = rd(fh)
        rho, v = h['dens'], h['velx']
        eint = h['eint']
        T = eint*GM1/rho
        E, F = m['m1_e'], m['m1_f1']
        out.append((t, E, F, v, T, rho*v + F/C**2, 0.5*rho*v*v + eint + E))
    return np.array(out)


if __name__ == '__main__':
    for d in sys.argv[1:]:
        s = series(d)
        print('== %s  (n=%d)' % (d, len(s)))
        print('     t        E         F1        v         T       P_tot     E_tot')
        idx = sorted(set([0, 1, 2, 4, 10, 20, 40, len(s)//2, len(s)-1]))
        for i in idx:
            if i < len(s):
                print('%7.3f ' % s[i, 0] + ' '.join('%9.6f' % x for x in s[i, 1:]))
