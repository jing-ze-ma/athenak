#!/usr/bin/env python3
"""Part 3: independent shell-mean resolved flux from a 3-D dump vs the closure's
F_conv_res at the same time.  One dump at a time (login-node memory)."""
import sys
import numpy as np

sys.path.insert(0, '/viper/u2/jinma/ATHENAK/athenak/vis/python')
import bin_convert as bc                                       # noqa: E402

D = '/viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r11/wedge9d/'
TURN = 4705.0
R = 2.3717e11


def read_blocks(path):
    out, t, rows = [], None, []
    for line in open(path):
        if line.startswith('# t '):
            if t is not None:
                out.append((t, np.array(rows)))
            t = float(line.split()[2])
            rows = []
        elif line.strip() and not line.startswith('#'):
            rows.append([float(x) for x in line.split()])
    out.append((t, np.array(rows)))
    return out


MLT = read_blocks(D + 'mltfaces_wedge9d.txt.ts')


def shell_sums(fn):
    fd = bc.read_binary(fn)
    nx1 = fd['Nx1']
    acc = np.zeros((6, nx1))
    for q in range(fd['n_mbs']):
        d = fd['mb_data']['dens'][q].astype(np.float64)
        vr = fd['mb_data']['velx'][q].astype(np.float64)
        v2 = (vr**2 + fd['mb_data']['vely'][q].astype(np.float64)**2
              + fd['mb_data']['velz'][q].astype(np.float64)**2)
        e = fd['mb_data']['eint'][q].astype(np.float64)
        eps = e / d
        acc[0] += d.sum(axis=(0, 1))
        acc[1] += e.sum(axis=(0, 1))
        acc[2] += (d * vr).sum(axis=(0, 1))
        acc[3] += (vr * e).sum(axis=(0, 1))          # = <rho vr eps>
        acc[4] += eps.sum(axis=(0, 1))
        acc[5] += (0.5 * d * vr * v2).sum(axis=(0, 1))
        n = d.shape[0] * d.shape[1]
    return fd['time'], acc / (n * fd['n_mbs']), nx1


for fn in sys.argv[1:]:
    t, s, nx1 = shell_sums(fn)
    j = int(np.argmin([abs(b[0] - t) for b in MLT]))
    tm, a = MLT[j]
    r = a[:, 1]
    ii = a[:, 0].astype(int)              # face index in code (is+1..ie), is = 3
    # chi = p/e per shell from the closure's face pressure and the dump's <e>
    ecell = s[1]
    # face-index i in the table corresponds to the face between cells i-1 and i;
    # dump arrays are 0-based cell index ic = i - 3 (nghost = 3)
    ic = ii - 3
    chi_f = a[:, 3] / (0.5 * (ecell[ic] + ecell[ic - 1]))
    # cell-level resolved flux, two brackets for chi
    print('\n### %s   t = %.1f s = %.4f turnovers  (mlt block t = %.1f)'
          % (fn.split('/')[-1], t, t / TURN, tm))
    print('%7s %7s %10s %10s %10s %10s %10s %9s' %
          ('r/R', 'chi', 'Fres_1/3', 'Fres_2/3', 'Fkin', 'Fres_clos',
           'F_req', 'ratio(2/3)'))
    for x in [0.55, 0.60, 0.65, 0.72, 0.80, 0.85, 0.88, 0.91, 0.93, 0.96, 0.99]:
        k = int(np.argmin(np.abs(r / R - x)))
        icc = ic[k]
        out = []
        for chi in (1.0 / 3.0, 2.0 / 3.0):
            fl = (1 + chi) * (s[3][icc - 1] - s[4][icc - 1] * s[2][icc - 1])
            fr = (1 + chi) * (s[3][icc] - s[4][icc] * s[2][icc])
            out.append(0.5 * (fl + fr))
        fkin = 0.5 * (s[5][icc - 1] + s[5][icc])
        fcl = a[k, 14]
        freq = a[k, 10]
        print('%7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %9.2f'
              % (r[k] / R, chi_f[k], out[0] + fkin, out[1] + fkin, fkin, fcl, freq,
                 (out[1] + fkin) / fcl if fcl != 0 else np.nan))
