"""Helpers for the S1 (implicit M1 on a spherical-polar wedge) gates.

load(): assemble an AthenaK .bin dump into whole-mesh arrays (k, j, i).
rgrid(): the radial faces and cell centroids exactly as CoordSphericalPolar builds
them (uniform, or the polynomial stretch use_grid_stretch_r_poly).
"""
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/wt_m1sp2/vis/python')
sys.path.insert(0, '/viper/ptmp2/jinma/wt_m1sp2/tests_m1')
import bin_convert  # noqa: E402
import common  # noqa: E402


def load_tab(fn):
    """a DOUBLE-precision radial slice (file_type = tab): arrays over i"""
    d = common.load_tab(fn)
    return {k: v[0, 0, :] for k, v in d.data.items()} | {'time': d.time}


def load(fn):
    d = bin_convert.read_binary(fn)
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    b1, b2, b3 = d['nx1_out_mb'], d['nx2_out_mb'], d['nx3_out_mb']
    out = {}
    for v in d['var_names']:
        a = np.zeros((n3, n2, n1))
        for m, blk in enumerate(d['mb_data'][v]):
            lx1, lx2, lx3 = d['mb_logical'][m][:3]
            a[lx3*b3:(lx3+1)*b3, lx2*b2:(lx2+1)*b2, lx1*b1:(lx1+1)*b1] = blk
        out[v] = a
    out['time'] = d['time']
    out['cycle'] = d['cycle']
    return out


def rgrid(r0, r1, n, poly=None):
    """faces rf (n+1) and centroids x1v (n) of the radial grid"""
    xi = np.arange(n + 1)/n
    u = xi.copy()
    if poly is not None:
        xik = xi.copy()
        for c in poly:
            u += c*xik*(1.0 - xi)
            xik = xik*xi
    rf = r0 + (r1 - r0)*u
    rl, rr = rf[:-1], rf[1:]
    q = rl/rr
    x1v = 0.25*(q*q + 1.0)/((q*q + q + 1.0)/3.0)*(rr + rl)
    return rf, x1v


HE4_POLY = [0.134559, 4.942505, -8.950106, 4.403042]
