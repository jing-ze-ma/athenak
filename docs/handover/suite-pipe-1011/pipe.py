"""suite-1011: JSD12 Sect 5.2 crooked pipe (tophat), Cartesian, probe temperature histories.
Units: length cm, time 1e-8 s (c = 300), T in 0.5 keV, a_r T^4 = 1 at T = 1, gamma - 1 = 1/58.32 (c_v = 1e15 erg/g/keV).
T_gas = (gamma-1) eint / rho  (hydro_w eint is the internal energy density).
Probes (Gentile 2001, (z, r) -> (x, y)): (0.25, 0), (2.75, 0), (3.5, 1.25), (4.25, 0), (6.75, 0).
usage: python pipe.py <rundir> ...   -> <rundir>/probes.txt and a summary line
"""
import sys, os, glob
import numpy as np
sys.path.insert(0, os.environ.get('ATHENAK_VIS', '/raven/u/jinma/ATHENAK/wt-littests/vis/python'))
import bin_convert as bc
GM1 = 0.017147
PROBES = [(0.25, 0.0), (2.75, 0.0), (3.5, 1.25), (4.25, 0.0), (6.75, 0.0)]


def probes(d):
    fs = sorted(glob.glob(os.path.join(d, 'bin', '*.hydro_w.*.bin')) or glob.glob(os.path.join(d, '*.hydro_w.*.bin')))
    rows = []
    for f in fs:
        r = bc.read_binary_as_athdf(f)
        T = GM1 * r['eint'][0] / r['dens'][0]
        x, y = r['x1v'], r['x2v']
        vals = []
        for (px, py) in PROBES:
            i = np.argmin(np.abs(x - px)); j = np.argmin(np.abs(y - py))
            jj = j - 1 if y[j] > py else j + 1          # average the two rows around y = py
            vals.append(0.5 * (T[j, i] + T[jj, i]))
        rows.append([r['Time']] + vals + [float(np.abs(r['velx'][0]).max())])
    a = np.array(rows)
    np.savetxt(os.path.join(d, 'probes.txt'), a, header='t T(0.25,0) T(2.75,0) T(3.5,1.25) T(4.25,0) T(6.75,0) max|vx|')
    return a


if __name__ == '__main__':
    for d in sys.argv[1:]:
        a = probes(d)
        if len(a) == 0:
            print(d, 'no dumps'); continue
        # time at which each probe first reaches T = 0.5 (0.25 keV)
        th = []
        for k in range(1, 6):
            idx = np.where(a[:, k] >= 0.5)[0]
            th.append('%.2f' % a[idx[0], 0] if len(idx) else '>%.1f' % a[-1, 0])
        print(os.path.basename(d.rstrip('/')), 't_end=%.2f' % a[-1, 0], 'T_end=' + ' '.join('%.3f' % v for v in a[-1, 1:6]),
              't(T>=0.5)=' + ' '.join(th), 'max|v|=%.1e' % a[:, 6].max())
