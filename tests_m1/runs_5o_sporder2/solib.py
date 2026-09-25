"""runs_5o_sporder2: helpers for the spherical-polar wedge convergence studies.

mkinput(): an athinput for the sph_shell / sph_atm pgen on the wedge from dicts.
run(): one CPU run (mpirun -np np) in its own directory.
tab(): a double-precision radial slice.  binall(): a float32 bin dump as (k, j, i).
restrict1(): volume-weighted 2:1 restriction along r (nested grids, uniform or poly).
order(): orders from successive errors.
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'vis', 'python'))
import common  # noqa: E402

W = '/viper/ptmp2/jinma/sporder2_0925'
POLY = {'f_stretch_r_c1': '+0.134559', 'f_stretch_r_c2': '+4.942505',
        'f_stretch_r_c3': '-8.950106', 'f_stretch_r_c4': '+4.403042'}


def mkinput(path, mesh, mb, time, rad, prob, outs, extra=None):
    """write an athinput; every argument is an ordered dict of key: value"""
    blocks = [('job', {'basename': 'so'}), ('mesh', mesh), ('meshblock', mb),
              ('time', time),
              ('hydro', {'eos': 'ideal', 'reconstruct': 'plm', 'rsolver': 'hllc',
                         'gamma': '1.666666666666667'}),
              ('rad_m1', rad), ('problem', prob)]
    if extra:
        blocks += extra
    for n, o in enumerate(outs):
        blocks.append(('output%d' % (n + 1), o))
    with open(path, 'w') as f:
        for b, d in blocks:
            f.write('<%s>\n' % b)
            for k, v in d.items():
                f.write('%-24s = %s\n' % (k, v))
            f.write('\n')


def run(binary, inp, rdir, np_=1, args=()):
    os.makedirs(rdir, exist_ok=True)
    cmd = ['mpirun', '-np', str(np_), '--oversubscribe', '--bind-to', 'none', binary,
           '-i', inp, '-d', rdir] + list(args)
    with open(os.path.join(rdir, 'log'), 'w') as lg:
        return subprocess.call(cmd, stdout=lg, stderr=subprocess.STDOUT)


def tab(fn):
    d = common.load_tab(fn)
    out = {k: v[0, 0, :] for k, v in d.data.items()}
    out['time'] = d.time
    out['x1'] = d.x1
    return out


def binall(fn):
    import bin_convert
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
    return out


def rfaces(r0, r1, n, poly=False):
    """radial faces as CoordSphericalPolar builds them (uniform or poly stretch)"""
    xi = np.arange(n + 1)/n
    u = xi.copy()
    if poly:
        xik = xi.copy()
        for m in range(4):
            u += float(POLY['f_stretch_r_c%d' % (m + 1)])*xik*(1.0 - xi)
            xik = xik*xi
    return r0 + (r1 - r0)*u


def restrict_r(a, rf):
    """2:1 volume-weighted restriction along the last axis (r); rf fine faces"""
    w = np.diff(rf**3)
    s = a*w
    return (s[..., 0::2] + s[..., 1::2])/(w[0::2] + w[1::2])


def order(e):
    e = np.asarray(e, dtype=float)
    return np.log2(e[:-1]/e[1:])
