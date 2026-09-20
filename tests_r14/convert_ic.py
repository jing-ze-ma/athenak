#!/usr/bin/env python3
"""tests_r14: convert an ic_profile built with the OLD EOS radiation taper to the NEW,
density-only, very-low-window taper.

The IC file stores (r, rho, eint) with
    eint = e_gas(rho,T) + w_old(rho,T) a T^4,
w_old = max(w_rho, w_T) the gated taper (rho 4.97e-10..1.93e-9, T 4.52e4..6.33e4 K).
The rho and T of the star are NOT changed: T is recovered from (rho, eint) with the OLD
weight (relax_ic.temp_of, which matches the run to dT 7e-5), and the file is rewritten
with
    eint_new = e_gas(rho,T) + w_new(rho) a T^4,
w_new the SAME cubic smoothstep of rad_taper::Weight in log10 rho, but with the new
window and NO temperature gate (eos_rad_t_hi = 0).

usage: convert_ic.py [in.txt] [out.txt] [--lo RHO_LO] [--hi RHO_HI]
      --lo / --hi override the new window (defaults 1e-11 / 1e-10).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'tests_r12'))
import relax_ic as R   # noqa: E402   (egas, taper_w, temp_of, A_RAD, RSTAR)

RHO_LO_NEW, RHO_HI_NEW = 1e-11, 1e-10


def w_new(rho):
    """rad_taper::Weight(log10 rho, log10 RHO_LO_NEW, log10 RHO_HI_NEW): 0 at and below
    rho_lo, 1 at and above rho_hi, the C1 cubic s^2(3-2s) in between."""
    s = np.clip((np.log10(rho) - np.log10(RHO_LO_NEW))
                / (np.log10(RHO_HI_NEW) - np.log10(RHO_LO_NEW)), 0.0, 1.0)
    return s*s*(3.0 - 2.0*s)


def main():
    global RHO_LO_NEW, RHO_HI_NEW
    av = sys.argv[1:]
    for key, nm in (('--lo', 'RHO_LO_NEW'), ('--hi', 'RHO_HI_NEW')):
        if key in av:
            k = av.index(key)
            globals()[nm] = float(av[k+1])
            del av[k:k+2]
    print('new window: rho_lo %.3e  rho_hi %.3e' % (RHO_LO_NEW, RHO_HI_NEW))
    fin = av[0] if len(av) > 0 else R.IC0
    fout = av[1] if len(av) > 1 else os.path.join(HERE, 'ic_he4_neww.txt')
    d = np.loadtxt(fin, comments='#')
    r, rho, eint = d[:, 0], d[:, 1], d[:, 2]

    T = R.temp_of(rho, eint)
    wo = R.taper_w(rho, T)
    wn = w_new(rho)
    eg = R.egas(rho, T)
    en = eg + wn*R.A_RAD*T**4
    # round-trip of the UNCHANGED rows: eint_old rebuilt from the recovered T
    eo = eg + wo*R.A_RAD*T**4

    both1 = (wo >= 1.0) & (wn >= 1.0)
    chg = np.abs(en/eint - 1.0)
    print('rows %d ; w_old==w_new==1 on %d of them' % (len(r), both1.sum()))
    print('  max |eint_new/eint_old - 1| where both weights are 1 : %.3e' % chg[both1].max())
    print('  (round-trip only, e_gas+w_old aT^4 vs the file)       : %.3e'
          % np.abs(eo[both1]/eint[both1] - 1.0).max())
    m = chg > 1e-12
    if m.any():
        print('  eint CHANGES on %d rows, r/R = %.4f .. %.4f (rho %.3e .. %.3e)'
              % (m.sum(), r[m].min()/R.RSTAR, r[m].max()/R.RSTAR,
                 rho[m].min(), rho[m].max()))
        print('  relative change there: min %.3e  max %.3e' % (chg[m].min(), chg[m].max()))
        for s in (0.95, 0.98, 0.99, 1.00, 1.01):
            k = int(np.argmin(np.abs(r/R.RSTAR - s)))
            print('    r/R %.3f  rho %.3e  T %.4e  w_old %.4f  w_new %.4f  '
                  'eint_new/eint_old %.6f' % (r[k]/R.RSTAR, rho[k], T[k], wo[k], wn[k],
                                              en[k]/eint[k]))
        k = len(r) - 1
        print('    top      rho %.3e  T %.4e  w_old %.4f  w_new %.4f  ratio %.6f'
              % (rho[k], T[k], wo[k], wn[k], en[k]/eint[k]))
    # where w_new < 1 at all
    mw = wn < 1.0
    if mw.any():
        print('  w_new < 1 above r/R = %.4f (%d rows)' % (r[mw].min()/R.RSTAR, mw.sum()))

    with open(fout, 'w') as fh:
        fh.write('# tests_r14/convert_ic.py: %s converted to the NEW EOS radiation\n'
                 '# taper, DENSITY ONLY, window rho %.3e (w=0) .. %.3e (w=1) g/cm^3,\n'
                 '# temperature gate OFF (eos_rad_t_hi = 0).  rho and T are unchanged;\n'
                 '# T was recovered from (rho, eint) with the OLD gated taper.\n'
                 '# r[cm]  rho[g/cm^3]  eint[erg/cm^3]\n'
                 % (os.path.basename(fin), RHO_LO_NEW, RHO_HI_NEW))
        for a, b, c in zip(r, rho, en):
            fh.write('%.10e %.10e %.10e\n' % (a, b, c))
    print('wrote', fout)


if __name__ == '__main__':
    main()
