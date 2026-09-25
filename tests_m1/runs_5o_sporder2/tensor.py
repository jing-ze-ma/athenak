#!/usr/bin/env python3
"""vet_col tensor accuracy on the wedge: f_K of the FIRST build (one step, the analytic
initial state) against a fine reference build, per radial resolution.

usage: tensor.py run|eval <tag> <case> [variant ...]
case: mid (Gaussian shell, rho kappa_s = 1, S = E, reflect r) | ex (sph_atm, chi =
10/r^2, flux in, Marshak out: the runs_5d T-S6 (ii) state)
variant: nc<ncore>p<nsub>[Q2] (Q2: vet_col_order2 = true); the reference is n = 2048,
nc64p4 (env SO_REFV) of the same tag.
"""
import os
import sys

import numpy as np

import solib as so

NS = [32, 64, 128, 256, 512]


def mk(tag, case, n, var):
    nc = int(var[2:var.index('p')])
    rest = var[var.index('p') + 1:]
    npp = int(rest.split('Q')[0])
    d = os.path.join(so.W, 'cpu', 'tensor', tag, '%s_%s_%d' % (case, var, n))
    os.makedirs(d, exist_ok=True)
    r1 = 3.0 if case == 'mid' else 5.0
    mesh = {'nghost': 2, 'nx1': n, 'x1min': 1.0, 'x1max': r1,
            'ix1_bc': 'outflow', 'ox1_bc': 'outflow',
            'nx2': 4, 'x2min': 1.4707963267948966, 'x2max': 1.6707963267948966,
            'ix2_bc': 'periodic', 'ox2_bc': 'periodic',
            'nx3': 4, 'x3min': 0.0, 'x3max': 0.2, 'ix3_bc': 'periodic',
            'ox3_bc': 'periodic', 'use_spherical_polar': 'true'}
    mb = {'nx1': n, 'nx2': 4, 'nx3': 4}
    time = {'evolution': 'dynamic', 'integrator': 'rk2', 'cfl_number': 0.3,
            'nlim': 1, 'tlim': '1.0e30', 'ndiag': 100000}
    rad = {'c_light': 1.0 if case == 'mid' else 100.0, 'chat_over_c': 1.0,
           'closure': 'vet_col', 'e_floor': '1.0e-30', 'opacity': 'const',
           'kappa_p': 0.0, 'kappa_e': 0.0,
           'kappa_f': 0.0 if case == 'mid' else 10.0,
           'kappa_s': 1.0 if case == 'mid' else 0.0, 'arad': 1.0,
           'coupling': 'true', 'gas_feedback': 'false', 'transport': 'implicit',
           'time_scheme': 'be', 'implicit_cfl': '1.0e-3',
           'implicit_solver': 'bicgstab', 'implicit_precond': 'line',
           'implicit_tol': '1.0e-10', 'implicit_lin_tol': '1.0e-10',
           'vet_col_ncore': nc, 'vet_col_nsub': npp,
           'vet_col_dump': os.path.join(d, 'dump'), 'vet_col_dump_every': 1}
    if 'Q' in rest:
        rad['vet_col_order2'] = 'true'
    if case == 'mid':
        rad.update({'implicit_bc_x1min': 'reflect', 'implicit_bc_x1max': 'reflect'})
        prob = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_shell', 'gas_rho': 1.0,
                'gas_temp': '1.0e-16', 'e_out': '1.0e-2', 'e_amp': 1.0, 'r0': 2.0,
                'width': 0.2}
    else:
        rad.update({'implicit_bc_x1min': 'flux', 'implicit_flux_x1min': 1.0,
                    'implicit_bc_x1max': 'marshak', 'implicit_ebath_x1max': 0.0,
                    'marshak_q': 0.5})
        prob = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_atm', 'atm_rho0': 1.0,
                'atm_rho_n': 2.0, 'atm_temp': '1.0e-3', 'atm_init': 'eddington',
                'e_out': '1.0e-6', 'user_srcs': 'true'}
    inp = os.path.join(d, 'in.athinput')
    so.mkinput(inp, mesh, mb, time, rad, prob, [])
    return d, inp


def load(d):
    a = np.loadtxt(os.path.join(d, 'dump.00000.txt'))
    return a[:, 0], a[:, 6]


def main():
    act, tag, case = sys.argv[1:4]
    vs = sys.argv[4:] or ['nc8p1']
    if act == 'run':
        b = os.path.join(so.W, 'bin', 'athena_%s_none_cpu' % tag)
        jobs = ([(2048, os.environ.get('SO_REFV', 'nc64p4'))]
                + [(n, v) for v in vs for n in NS])
        for n, v in jobs:
            d, inp = mk(tag, case, n, v)
            if not os.path.exists(os.path.join(d, 'dump.00000.txt')):
                print(n, v, 'rc', so.run(b, inp, d), flush=True)
        return
    rr, fr = load(mk(tag, case, 2048, os.environ.get('SO_REFV', 'nc64p4'))[0])
    for v in vs:
        e1, ei = [], []
        for n in NS:
            r, f = load(mk(tag, case, n, v)[0])
            e = f - np.interp(r, rr, fr)
            e1.append(np.mean(np.abs(e)))
            ei.append(np.max(np.abs(e)))
        print('%s %s %-10s L1 %s  ord %s' % (tag, case, v,
              ' '.join('%.2e' % x for x in e1),
              ' '.join('%.2f' % x for x in so.order(e1))))
        print('%26s Linf %s  ord %s' % ('', ' '.join('%.2e' % x for x in ei),
              ' '.join('%.2f' % x for x in so.order(ei))))


if __name__ == '__main__':
    main()
