#!/usr/bin/env python3
"""AG Car IC column extended INWARD (agcmodel-1012, worker AC): a deeper inner wall below the Fe convection zone.

The rows of --base (e.g. ic_agcar_A_ge.txt) are kept VERBATIM.  From its first (deepest) row (r0, T0, rho0) the same
equations the column was built with (make_ic_agcar_mltpp.py rhs2: dPtot/dr = -rho GM/r^2, dlnT/dlnPtot = nabla, nabla
= MLT++ gradient, = nabla_rad where radiative) are integrated INWARD on the log grid continued below r0 (dln r = --dlnr)
down to --rin.  New rows: r rho eint F_r E=aT^4 T F_MLT/F (as the file).  Same stellar model (L, M, Teff, EOS, kappa_R,
MLT alpha, MLT++ f1/f2): only the stable radiative buffer below the CZ is made thicker.  Prints H_p counts.

usage: make_ic_agcar_deep.py --base FILE --rin R[cm] --name OUT [--lum --gm --teff --rstar]
"""
import argparse
import os
import sys
import types

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_ic_mlt_star as M                     # noqa: E402
from make_ic_agcar_mltpp import StarPP           # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--base', required=True)
    p.add_argument('--rin', type=float, required=True)
    p.add_argument('--name', required=True)
    p.add_argument('--out', default='.')
    p.add_argument('--lum', type=float, default=3.41170859e39)      # AG Car A (file header)
    p.add_argument('--gm', type=float, default=7.29918420e27)
    p.add_argument('--teff', type=float, default=9000.0)
    p.add_argument('--rstar', type=float, default=2.70140777e13)
    p.add_argument('--eos', default='table:/raven/ptmp/jinma/agcar_1008/shakeB_ana/eos_dump.txt')
    p.add_argument('--ross', default='/raven/ptmp/jinma/agcar_1008/files/rosseland_tops_gs98_x0.36_z0.02.txt')
    p.add_argument('--alpha', type=float, default=1.5)
    p.add_argument('--f1', type=float, default=1.0e-4)
    p.add_argument('--f2', type=float, default=1.0e-3)
    p.add_argument('--dlnr', type=float, default=2.0e-4)
    a = p.parse_args()
    st = StarPP(types.SimpleNamespace(gm=a.gm, lum=a.lum, rstar=a.rstar, teff=a.teff, eos=a.eos, ross=a.ross,
                                      ld_ext=None, alpha=a.alpha, mltpp=1, f1=a.f1, f2=a.f2))
    GM, LUM, eos = st.GM, st.LUM, st.eos
    hdr, rows = [], []
    with open(a.base) as fh:
        for ln in fh:
            (hdr if ln.startswith('#') else rows).append(ln)
    r0, rho0 = [float(v) for v in rows[0].split()[:2]]
    T0 = float(rows[0].split()[5])

    def rhs2(rr, y):
        T_, rho_ = np.exp(y[0]), np.exp(y[1])
        nab = st.mlt_gradient(T_, rho_, rr)[0]
        P_, pg_, chir, chit, cp, dlt, nad_ = st.thermo(T_, rho_)
        dlnP = -rho_*GM/rr**2/P_
        dlnT = nab*dlnP
        return [dlnT, (dlnP - chit*dlnT)/chir]
    nk = int(np.floor(np.log(r0/a.rin)/a.dlnr)) + 1
    rn = r0*np.exp(-a.dlnr*np.arange(nk + 1))
    s = solve_ivp(rhs2, [r0, rn[-1]], [np.log(T0), np.log(rho0)], t_eval=rn, method='LSODA', rtol=1e-10, atol=1e-13)
    assert s.success, s.message
    r = s.t[1:][::-1]
    T = np.exp(s.y[0][1:])[::-1]
    rho = np.exp(s.y[1][1:])[::-1]
    fm = np.zeros(r.size); ffr = np.ones(r.size); Hp = np.zeros(r.size); nrad = np.zeros(r.size); nad = np.zeros(r.size)
    for i in range(r.size):
        o = st.mlt_gradient(T[i], rho[i], r[i])
        ffr[i], fm[i], Hp[i], nad[i], nrad[i] = o[1], o[2], o[6], o[7], o[8]
    eint = rho*10.0**eos.le(np.log10(T), np.log10(rho))
    Fr = ffr*LUM/(4*np.pi*r**2)
    E = M.A_RAD_*T**4
    nH = np.cumsum((np.gradient(r)/Hp)[::-1])[::-1]
    print('base row r0 %.10e T0 %.6e rho0 %.6e; new rows %d down to r %.10e (%.2f Rsun): T %.4e rho %.4e'
          % (r0, T0, rho0, r.size, r[0], r[0]/M.RSUN, T[0], rho[0]))
    print('new rows: max nabla_rad/nabla_ad %.3f (convective rows %d), max F_MLT/F %.3e; H_p(r_in) %.3e (%.2f Rsun);'
          ' H_p from r_in to r0: %.3f' % ((nrad/nad).max(), int((nrad > nad).sum()), fm.max(), Hp[0], Hp[0]/M.RSUN, nH[0]))
    for x in np.linspace(r[0], r0, 8):
        i = int(np.argmin(abs(r - x)))
        print('  r %.4e (%.1f Rsun) T %.4e rho %.4e nrad/nad %.3f Hp %.2f Rsun' % (r[i], r[i]/M.RSUN, T[i], rho[i],
                                                                                 nrad[i]/nad[i], Hp[i]/M.RSUN))
    fn = os.path.join(a.out, a.name)
    with open(fn, 'w') as fh:
        fh.write(hdr[0])
        for h in hdr[1:-1]:
            fh.write(h)
        fh.write('# agcmodel-1012 (make_ic_agcar_deep.py): %d rows r < %.10e added INWARD to %.10e (%.2f Rsun) with the '
                 'column equations (dPtot/dr = -rho GM/r^2, MLT++ nabla; radiative there, max F_MLT/F %.1e); rows from r0 '
                 'up copied verbatim from %s\n' % (r.size, r0, r[0], r[0]/M.RSUN, fm.max(), a.base))
        fh.write(hdr[-1])
        for row in zip(r, rho, eint, Fr, E, T, fm):
            fh.write('%.10e %.10e %.10e %.10e %.10e %.10e %.6e\n' % row)
        for ln in rows:
            fh.write(ln)
    print('wrote', fn, 'F_r(r_in) = %.10e  (implicit_flux_x1min for x1min = r_in: L/(4 pi r^2) = %.10e)'
          % (Fr[0], LUM/(4*np.pi*r[0]**2)))


if __name__ == '__main__':
    main()
