#!/usr/bin/env python3
"""phi-averaged latitude profiles of v_r and v_theta (km/s) in the FeCZ shells of the
he_star_m1 3-D wedge (hydro_w bin dumps), plus the rank-0 theta index of the Picard
NON-CONVERGED worst cells in the run logs (2 ranks, 4 MeshBlocks: rank 0 = m 0, 1).
usage: lat_prof.py INPUT RUNDIR[,RUNDIR2] DUMPNUM [DUMPNUM ...] [--rlo 0.7 --rhi 0.93]"""
import glob
import re
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/vis/python')
sys.path.insert(0, '/viper/ptmp2/jinma/wt_hepresn/tests_hepresn')
from bin_convert import read_binary  # noqa: E402
import analyze_gate_mlt as A  # noqa: E402
import analyze_3d as B  # noqa: E402

JS = [0, 4, 8, 16, 24, 31, 32, 40, 48, 56, 59, 63]


def arg(name, default):
    return type(default)(sys.argv[sys.argv.index(name) + 1]) \
        if name in sys.argv else default


def prof(fn, sel):
    d = read_binary(fn)
    geo = np.asarray(d['mb_geometry'])
    v1 = np.asarray(d['mb_data']['velx'])
    v2 = np.asarray(d['mb_data']['vely'])
    nmb, n3, n2, n1 = v1.shape
    t0 = sorted(set(np.round(geo[:, 2], 8)))
    ntb = len(t0) * n2
    pr, pt, cnt = np.zeros(ntb), np.zeros(ntb), np.zeros(ntb)
    for m in range(nmb):
        j0 = t0.index(np.round(geo[m, 2], 8)) * n2
        pr[j0:j0 + n2] += v1[m][:, :, sel].mean(axis=(0, 2))
        pt[j0:j0 + n2] += v2[m][:, :, sel].mean(axis=(0, 2))
        cnt[j0:j0 + n2] += 1
    return d['time'], pr / cnt / 1e5, pt / cnt / 1e5


def nc_theta(run):
    jg = []
    for fn in glob.glob(run + '/run*.log'):
        for m, j in re.findall(r'Picard NON-CONVERGED.*?\(m,k,j,i\)=\((\d+),\d+,(\d+),',
                               open(fn, errors='replace').read()):
            jg.append(int(m) * 32 + int(j) - 3)
    return np.bincount(np.array(jg, int), minlength=64) if jg else None


def main():
    inp = sys.argv[1]
    runs = sys.argv[2].split(',')
    nums = [int(a) for a in sys.argv[3:] if a.isdigit()]
    par = B.read_mesh(inp)
    re_ = A.edges(par, int(par['nx1']))
    x = A.centroid(re_[:-1], re_[1:]) / A.RSTAR
    sel = (x > arg('--rlo', 0.7)) & (x < arg('--rhi', 0.93))
    print('theta index:', JS)
    for n in nums:
        for run in runs:
            fn = '%s/bin/hepresn.hydro_w.%05d.bin' % (run, n)
            try:
                t, pr, pt = prof(fn, sel)
            except (OSError, KeyError):
                continue
            name = run.split('/')[-1]
            print('%-10s dump %2d t %6.0f v_r  %s'
                  % (name, n, t, np.round(pr[JS], 1)))
            print('%-10s                  v_th %s' % ('', np.round(pt[JS], 1)))
    for run in runs:
        h = nc_theta(run)
        if h is not None:
            nz = ' '.join('%d:%d' % (q, c) for q, c in enumerate(h) if c)
            print('%s NC worst-cell theta index (rank 0): %s' % (run, nz))


if __name__ == '__main__':
    main()
