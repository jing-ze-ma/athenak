#!/usr/bin/env python3
"""Steady grey spherical atmosphere on the wedge (sph_atm, pure scattering), Eddington,
against the ANALYTIC continuum solution (flux in at r_in, Marshak out at r_out):

  E(r) = 3 rho0 k F_in r_in^(n+2)/(c (n+1)) (r^-(n+1) - r_out^-(n+1)) + F_out/(c q)

usage: steady.py run|eval <tag> <case> [n ...];  case = edd_u | edd_s (stretched r)
r = 1..5, rho kappa_s = 100 (r/r_in)^-2 (tau 80), c = 100, F_in = 1, q = 0.5.
Error: L1_V(E/E_exact - 1) at the cell centroids; also the top-cell error.
Extra rad keys: env SO_RAD="k=v,k=v".
"""
import os
import sys

import numpy as np

import solib as so

R0, R1, C, K, Q, FIN, NRHO = 1.0, 5.0, 100.0, 100.0, 0.5, 1.0, 2.0


def mk(tag, case, n):
    poly = case.endswith('_s')
    d = os.path.join(so.W, 'cpu', 'steady', tag,
                     '%s_%d%s' % (case, n, os.environ.get('SO_SUF', '')))
    os.makedirs(d, exist_ok=True)
    mesh = {'nghost': 2, 'nx1': n, 'x1min': R0, 'x1max': R1,
            'ix1_bc': 'outflow', 'ox1_bc': 'outflow',
            'nx2': 4, 'x2min': 1.4707963267948966, 'x2max': 1.6707963267948966,
            'ix2_bc': 'periodic', 'ox2_bc': 'periodic',
            'nx3': 4, 'x3min': 0.0, 'x3max': 0.2, 'ix3_bc': 'periodic',
            'ox3_bc': 'periodic', 'use_spherical_polar': 'true'}
    if poly:
        mesh['use_grid_stretch_r_poly'] = 'true'
        mesh.update(so.POLY)
    mb = {'nx1': n, 'nx2': 4, 'nx3': 4}
    time = {'evolution': 'dynamic', 'integrator': 'rk2', 'cfl_number': 0.3,
            'nlim': int(os.environ.get('SO_NLIM', 40)), 'tlim': '1.0e30',
            'ndiag': 100000}
    rad = {'c_light': C, 'chat_over_c': 1.0, 'closure': 'eddington',
           'e_floor': '1.0e-30', 'opacity': 'const', 'kappa_p': 0.0, 'kappa_e': 0.0,
           'kappa_f': 0.0, 'kappa_s': K, 'arad': 1.0, 'coupling': 'true',
           'gas_feedback': 'false', 'transport': 'implicit', 'time_scheme': 'be',
           'implicit_cfl': '1.0e6', 'implicit_solver': 'bicgstab',
           'implicit_maxit': 100, 'implicit_precond': 'line',
           'implicit_tol': '1.0e-10', 'implicit_lin_tol': '1.0e-10',
           'implicit_bc_x1min': 'flux', 'implicit_flux_x1min': FIN,
           'implicit_bc_x1max': 'marshak', 'implicit_ebath_x1max': 0.0,
           'marshak_q': Q}
    for kv in filter(None, os.environ.get('SO_RAD', '').split(',')):
        k, v = kv.split('=')
        rad[k] = v
    prob = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_atm', 'atm_rho0': 1.0,
            'atm_rho_n': NRHO, 'atm_temp': '1.0e-3', 'atm_init': 'eddington',
            'e_out': '1.0e-6', 'user_srcs': 'true'}
    outs = [{'file_type': 'tab', 'variable': 'm1', 'data_format': '%24.16e',
             'slice_x2': 1.5707963267948966, 'slice_x3': 0.1,
             'dcycle': int(os.environ.get('SO_NLIM', 40))}]
    inp = os.path.join(d, 'in.athinput')
    so.mkinput(inp, mesh, mb, time, rad, prob, outs)
    return d, inp


def exact(r):
    n = NRHO
    return (3.0*K*FIN*R0**(n + 2)/(C*(n + 1))*(r**-(n + 1) - R1**-(n + 1))
            + FIN*R0**2/(R1**2*C*Q))


def main():
    act, tag, case = sys.argv[1:4]
    ns = [int(x) for x in sys.argv[4:]] or [32, 64, 128, 256, 512]
    if act == 'run':
        b = os.path.join(so.W, 'bin', 'athena_%s_none_cpu' % tag)
        for n in ns:
            d, inp = mk(tag, case, n)
            print(case, n, 'rc', so.run(b, inp, d), flush=True)
        return
    errs, tops = [], []
    for n in ns:
        d, _ = mk(tag, case, n)
        fs = sorted(f for f in os.listdir(os.path.join(d, 'tab')) if f.endswith('.tab'))
        t = so.tab(os.path.join(d, 'tab', fs[-1]))
        rf = so.rfaces(R0, R1, n, case.endswith('_s'))
        rl, rr = rf[:-1], rf[1:]
        rc = 0.75*(rr**4 - rl**4)/(rr**3 - rl**3)
        # the exact CELL AVERAGE (volume), by 8-point Gauss in r^2 dr
        g, wg = np.polynomial.legendre.leggauss(8)
        ea = np.zeros(n)
        for gg, ww in zip(g, wg):
            r = 0.5*(rl + rr) + 0.5*(rr - rl)*gg
            ea += ww*0.5*(rr - rl)*exact(r)*r*r
        ea /= (rr**3 - rl**3)/3.0
        w = rr**3 - rl**3
        e = t['m1_e']
        errs.append(np.sum(w*np.abs(e/ea - 1.0))/np.sum(w))
        tops.append(abs(e[-1]/ea[-1] - 1.0))
        del rc
    print('%s %s: L1 %s' % (tag, case, ' '.join('%.3e' % x for x in errs)))
    print('   orders %s' % ' '.join('%.2f' % p for p in so.order(errs)))
    print('   top-cell %s' % ' '.join('%.3e' % x for x in tops))


if __name__ == '__main__':
    main()
