"""scho-1009 shadow metrics on E (M1) and J_fs (vet_sc formal solution) vs exact
usage: shdfs.py <exact.npz> <rundir> [...]   (env SHD_YC centre y, SHD_LIT lit offset)
prints per x: depth = v(yc)/v(yc+lit), 10-90 edge width on y>yc side of v, both for E, Jfs
and the exact J; also abs J_fs/J_exact at the lit point and min J_fs (positivity)"""
import glob
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.environ['REPO'], 'vis/python'))
import bin_convert as bc  # noqa: E402
YC = float(os.environ.get('SHD_YC', '1.5'))
LIT = float(os.environ.get('SHD_LIT', '0.3'))
YMAX = float(os.environ.get('SHD_YMAX', '9.0'))


def edge(y, v, lo, hi):
    s = (y > YC) & (y < YC + YMAX)
    yy, vv = y[s], v[s]
    vn = (vv - vv.min())/max(vv.max() - vv.min(), 1e-300)

    def cross(level):
        k = np.argmax(vn >= level)
        if k == 0:
            return yy[0]
        return yy[k-1] + (level - vn[k-1])*(yy[k] - yy[k-1])/(vn[k] - vn[k-1])
    return cross(hi) - cross(lo)


ex = np.load(sys.argv[1])
J, x, y = ex['J'], ex['x'], ex['y']
for d in sys.argv[2:]:
    f = sorted(glob.glob(d + '/bin/*.m1_vsc.*.bin'))[-1]
    v = bc.read_binary_as_athdf(f)

    Jf = v['m1_vj'][0]
    m = bc.read_binary_as_athdf(sorted(glob.glob(d + '/bin/*.m1.*.bin'))[-1])
    E = m['m1_e'][0]
    print(d, os.path.basename(f), 'minJfs %.3e' % Jf.min())
    for xc in (0.65, 0.85):
        i = int(np.argmin(np.abs(x - xc)))
        jl = int(np.argmin(np.abs(y - (YC + LIT))))
        jc = int(np.argmin(np.abs(y - YC)))
        out = '  x=%.2f' % x[i]
        for nm, a in (('E', E), ('Jfs', Jf), ('ex', J)):
            out += '  %s depth %.3f edge %.3f' % (nm, a[jc, i]/a[jl, i],
                                                   edge(y, a[:, i], 0.1, 0.9))
        out += '  Jfs/ex lit %.3f E/ex lit %.3f' % (Jf[jl, i]/J[jl, i], E[jl, i]/J[jl, i])
        print(out)
