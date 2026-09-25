#!/usr/bin/env python3
"""Radial convergence on the wedge (spherically symmetric Gaussian shell, sph_shell).

usage: radial.py run|eval <binary-tag> <case> [n ...]
cases (r = 1..3, theta = pi/2 +- 0.4 reflect, phi 0..0.8 periodic, 4 x 4 angular):
  thick_u / thick_s   rho kappa_s = 1000, t = 200, implicit_cfl 50 (n = 32)
  thin_u / thin_s     rho kappa_s = 1e-3, t = 0.6, implicit_cfl 0.5 (n = 32)
  *_vc                the same with closure = vet_col (default eddington)
dt scales with dx (implicit_cfl fixed): hesdirk2 is second order in time, so the
measured order is min(space, time).  Self-convergence: e_n = L1_V(E_n - R E_2n) / L1_V(E).
Extra rad keys: env SO_RAD="k=v,k=v".
"""
import os
import sys

import numpy as np

import solib as so

CASES = {
    'thick': dict(kap='1000.0', T=200.0, cfl=50.0, e_out='1.0e-2'),
    'thin': dict(kap='1.0e-3', T=0.6, cfl=0.5, e_out='1.0e-2'),
    'mid': dict(kap='1.0', T=1.0, cfl=0.5, e_out='1.0e-2'),
    # moving gas (the enthalpy flux): v_r = 0.01 sin(pi (r-1)/2), T_gas = 1e-3
    # the same with a Marshak (vacuum) outer face, the column top vet_col assumes
    'midm': dict(kap='1.0', T=1.0, cfl=0.5, e_out='1.0e-2',
                 rad={'implicit_bc_x1max': 'marshak', 'marshak_q': 0.5,
                      'implicit_ebath_x1max': 0.0}),
    'vthick': dict(kap='1000.0', T=50.0, cfl=50.0, e_out='1.0e-2',
                   prob={'gas_vr': 0.01, 'gas_temp': '1.0e-3'}),
}


def params(case):
    base, rest = case.split('_', 1)
    poly = rest.startswith('s')
    vc = rest.endswith('vc')
    return CASES[base], poly, vc


def mk(tag, case, n):
    c, poly, vc = params(case)
    d = os.path.join(so.W, 'cpu', 'radial', tag,
                     '%s_%d%s' % (case, n, os.environ.get('SO_SUF', '')))
    os.makedirs(d, exist_ok=True)
    mesh = {'nghost': 2, 'nx1': n, 'x1min': 1.0, 'x1max': 3.0,
            'ix1_bc': 'outflow', 'ox1_bc': 'outflow',
            'nx2': 4, 'x2min': 1.1707963267948966, 'x2max': 1.9707963267948966,
            'ix2_bc': 'reflect', 'ox2_bc': 'reflect',
            'nx3': 4, 'x3min': 0.0, 'x3max': 0.8, 'ix3_bc': 'periodic',
            'ox3_bc': 'periodic', 'use_spherical_polar': 'true'}
    if poly:
        mesh['use_grid_stretch_r_poly'] = 'true'
        mesh.update(so.POLY)
    if case.endswith('_c'):
        # the Cartesian twin (x1 = r, periodic transverse)
        mesh.update({'x2min': 0.0, 'x2max': 0.8, 'ix2_bc': 'periodic',
                     'ox2_bc': 'periodic', 'use_spherical_polar': 'false'})
    if 'prob' in c:
        # moving gas: reflecting hydro walls (v_r = 0 there, odd about the wall)
        mesh.update({'ix1_bc': 'reflect', 'ox1_bc': 'reflect'})
    mb = {'nx1': n, 'nx2': 4, 'nx3': 4}
    time = {'evolution': 'dynamic', 'integrator': 'rk2', 'cfl_number': 0.3,
            'nlim': 100000, 'tlim': c['T'], 'ndiag': 100000}
    rad = {'c_light': 1.0, 'chat_over_c': 1.0,
           'closure': 'vet_col' if vc else 'eddington', 'e_floor': '1.0e-30',
           'opacity': 'const', 'kappa_p': 0.0, 'kappa_f': 0.0, 'kappa_s': c['kap'],
           'arad': 1.0, 'coupling': 'true', 'gas_feedback': 'false',
           'transport': 'implicit', 'time_scheme': 'hesdirk2',
           'implicit_cfl': c['cfl'],
           'implicit_solver': 'bicgstab', 'implicit_maxit': 200,
           'implicit_precond': 'line', 'implicit_tol': '1.0e-12',
           'implicit_lin_tol': '1.0e-13',
           'implicit_bc_x1min': 'reflect', 'implicit_bc_x1max': 'reflect'}
    if os.environ.get('SO_DT'):
        # a FIXED dt for every n (time error frozen: the differences are spatial)
        rf = so.rfaces(1.0, 3.0, n, poly)
        rad['implicit_cfl'] = '%.12g' % (float(os.environ['SO_DT'])
                                         / min(np.diff(rf).min(), 0.18))
        d = d + '_dt' + os.environ['SO_DT']
        os.makedirs(d, exist_ok=True)
    rad.update(c.get('rad', {}))
    for kv in filter(None, os.environ.get('SO_RAD', '').split(',')):
        k, v = kv.split('=')
        rad[k] = v
    prob = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_shell', 'gas_rho': 1.0,
            'gas_temp': '1.0e-16', 'e_out': c['e_out'], 'e_amp': 1.0, 'r0': 2.0,
            'width': 0.2}
    prob.update(c.get('prob', {}))
    for kv in filter(None, os.environ.get('SO_PROB', '').split(',')):
        k, v = kv.split('=')
        prob[k] = v
    outs = [{'file_type': 'tab', 'variable': 'm1', 'data_format': '%24.16e',
             'slice_x2': 0.05 if case.endswith('_c') else 1.5207963267948967,
             'slice_x3': 0.05, 'dt': c['T']}]
    if 'prob' in c:
        outs.append(dict(outs[0], variable='hydro_w'))
    inp = os.path.join(d, 'in.athinput')
    so.mkinput(inp, mesh, mb, time, rad, prob, outs)
    return d, inp


def final(d):
    fs = sorted(f for f in os.listdir(os.path.join(d, 'tab')) if f.endswith('.tab'))
    return so.tab(os.path.join(d, 'tab', fs[-1]))


def main():
    act, tag, case = sys.argv[1:4]
    ns = [int(x) for x in sys.argv[4:]] or [32, 64, 128, 256, 512]
    _, poly, _ = params(case)
    if act == 'run':
        b = os.path.join(so.W, 'bin', 'athena_%s_none_cpu' % tag)
        for n in ns:
            d, inp = mk(tag, case, n)
            rc = so.run(b, inp, d)
            print(case, n, 'rc', rc, flush=True)
        return
    E = {}
    for n in ns:
        d, _ = mk(tag, case, n)
        t = final(d)
        E[n] = t
    errs = []
    for a, b in zip(ns[:-1], ns[1:]):
        rf = so.rfaces(1.0, 3.0, b, poly)
        rfa = so.rfaces(1.0, 3.0, a, poly)
        ec = E[a]['E'] if 'E' in E[a] else E[a][list(E[a])[0]]
        key = [k for k in E[a] if k not in ('time', 'x1')][0]
        ec = E[a][key]
        ef = so.restrict_r(E[b][key], rf)
        w = np.diff(rfa**3)
        errs.append(np.sum(w*np.abs(ec - ef))/np.sum(w*np.abs(ec)))
    print('%s %s key=%s t=%.6g' % (tag, case, key, E[ns[-1]]['time']))
    print('  L1 self-diff', ' '.join('%.3e' % e for e in errs))
    print('  orders      ', ' '.join('%.2f' % p for p in so.order(errs)))


if __name__ == '__main__':
    main()
