#!/usr/bin/env python3
"""tests_r17x: follow-ups to dt3d.py -- gas/radiation pressure fraction at the limiting
shell, the density-floor census, and whether the dt-setting structure is grid scale or
block-edge aligned.  READ-ONLY."""
import sys

import numpy as np

import dt3d as D

DF = D.DFLOOR


def extra(idx):
    d, g = D.load(idx)
    rc, dr, thc, dth, phc, dph = D.grid(d)
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    rho, eint = g['dens'], g['eint']
    print('\n########## dump %d  t=%.1f = %.4f turnover' % (idx, d['time'],
                                                            d['time']/D.TURN))
    # ---- density-floor census
    print('  global min rho = %.4e  (dfloor %.1e, ratio %.3f)'
          % (rho.min(), DF, rho.min()/DF))
    k, j, i = np.unravel_index(np.argmin(rho), rho.shape)
    print('    at i=%d r/R=%.4f  th=%.2f  ph=%.2f' % (i, rc[i]/D.RS,
                                                      np.degrees(thc[j]),
                                                      np.degrees(phc[k])))
    for f in (1.05, 1.2, 1.5, 2.0, 5.0):
        m = rho <= f*DF
        print('    cells rho <= %4.2f*dfloor : %8d' % (f, m.sum()), end='')
        if m.any():
            ii = np.where(m.any(axis=(0, 1)))[0]
            print('   shells i=%d..%d  r/R %.4f..%.4f' % (ii[0], ii[-1],
                                                          rc[ii[0]]/D.RS,
                                                          rc[ii[-1]]/D.RS))
        else:
            print()
    # ---- the deep base: gas vs radiation pressure and grid-scale test
    for i in (11, 13, 14):
        r2 = rho[:, :, i]
        T, p, g1, cs = D.eos_state(r2, eint[:, :, i])
        x = np.log10(r2)
        w, _ = D.weight(x)
        prad = w*D.ARAD*T**4/3.0
        beta = 1.0 - prad/p
        print('  shell %3d r/R %.4f : p_gas/p  min %.4f med %.4f max %.4f ;'
              ' Gamma1 min %.4f med %.4f ; T range %.4e..%.4e'
              % (i, rc[i]/D.RS, beta.min(), np.median(beta), beta.max(),
                 g1.min(), np.median(g1), T.min(), T.max()))
        # grid-scale content: rms of the 2-cell (Nyquist) component in theta and phi
        q = r2/r2.mean()
        for nm, ax in (('theta', 1), ('phi', 0)):
            nyq = 0.25*(np.roll(q, 1, ax) - 2*q + np.roll(q, -1, ax))
            print('      %-6s rms(q)=%.4f  rms(Nyquist part)=%.4f  ratio %.3f'
                  % (nm, q.std(), nyq.std(), nyq.std()/q.std()))
        # block-edge alignment of the 1% lowest-density cells
        lo = q <= np.percentile(q, 1.0)
        kk, jj = np.where(lo)
        for nm, a, nb in (('theta j%24', jj, 24), ('phi   k%24', kk, 24)):
            h = np.bincount(a % nb, minlength=nb)/len(a)*nb
            print('      %s occupancy of the 1%% lowest-rho cells: edge(0,23)='
                  '%.2f,%.2f  interior mean %.2f  max %.2f'
                  % (nm, h[0], h[-1], h[1:-1].mean(), h.max()))


if __name__ == '__main__':
    for a in ([int(x) for x in sys.argv[1:]] or [4, 5, 6]):
        extra(a)
