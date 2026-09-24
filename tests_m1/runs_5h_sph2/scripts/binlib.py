"""assemble AthenaK .bin dumps into whole-mesh arrays (k, j, i)"""
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_sph2/vis/python')
import bin_convert  # noqa: E402


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
