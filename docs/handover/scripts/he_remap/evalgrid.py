# flake8: noqa
"""evaluate a StretchRPoly parameter vector against the pp2 70 ks profiles (and the old 470 grid)"""
import importlib.util
import sys
import os
import numpy as np
os.environ.setdefault('FAC', '0.85')
sys.argv = sys.argv[:1] + sys.argv[1:]
spec = importlib.util.spec_from_file_location(
    'f3', '/viper/ptmp2/jinma/he_mltpp_1002/remap128/grid/fit3.py')
f3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f3)
R = f3.R
q = np.load('/viper/ptmp2/jinma/he_mltpp_1002/remap128/grid/prof70.npz')
Ht = q['Ht'].min(axis=0)
H5 = q['Hr5'].min(axis=0)
rvq = q['rv']


def ev(p, N, tag):
    rf = f3.faces(p, N)
    dr = np.diff(rf)
    rc = 0.5 * (rf[1:] + rf[:-1])
    h = np.interp(rc, rvq, Ht)
    h5 = np.interp(rc, rvq, H5)
    x = rc / R
    print('%s: N %d, dr_min %.5f R at %.4f R, max adjacent ratio %.3f' % (tag, N, dr.min(
    ) / R, x[np.argmin(dr)], np.max(np.maximum(dr[1:] / dr[:-1], dr[:-1] / dr[1:]))))
    for a, b in ((0.40, 0.45), (0.45, 0.55), (0.55, 0.60), (0.60, 0.70),
                 (0.70, 0.85), (0.85, 1.0), (1.0, 1.05), (1.05, 3.0)):
        m = (x >= a) & (x < b)
        print('   %.2f-%.2f R: %4d cells, dr/R %.5f-%.5f, min cells/Ht %.2f (at %.3f), min cells/Hray5 %.2f'
              % (a, b, m.sum(), dr[m].min() / R, dr[m].max() / R, (h / dr)[m].min(), x[m][np.argmin((h / dr)[m])],
                 (h5 / dr)[m].min()))
    return rf


if __name__ == '__main__':
    old = np.load('/viper/ptmp2/jinma/he_mltpp_1002/grid/p_470.npy')
    ev(old, 470, 'OLD p_470')
    for a in sys.argv[1:]:
        p = np.load(a)
        N = int(''.join(ch for ch in a.split('_')[1].split('.')[0] if ch.isdigit()))
        rf = ev(p, N, a)
        np.save(a.replace('p_', 'faces_'), rf)
