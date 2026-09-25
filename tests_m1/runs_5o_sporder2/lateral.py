#!/usr/bin/env python3
"""Lateral (theta / phi) convergence on the wedge: a smooth lateral mode on a uniform
background, r = 1..3 reflecting, theta = pi/2 +- 0.4 reflecting, phi 0..0.8 periodic.

usage: lateral.py run|eval <tag> <case> [s ...]   (s = 1 2 4 8: n1 = n_lat = 16 s)
cases: th_<reg>[_vc] (lat_l2 = 1, theta refined, 4 phi cells) | ph_<reg>[_vc]
       (lat_l3 = 1, phi refined, 4 theta cells); reg = thick (rho kappa_s 1000,
       t = 400) | thin (1e-3, t = 0.6) | mid (1, t = 1)
Self-convergence of the float32 bin dump: e_s = L1_V(E_s - R E_2s)/L1_V(E_s - <E_s>).
"""
import os
import sys

import numpy as np

import solib as so

REG = {'thick': dict(kap='1000.0', T=400.0, dtc=16.0),
       # the T-S4 atmosphere (sph_atm, pure scattering rho kappa_s = 100 (r/r_in)^-2,
       # r = 1..5, c = 100, flux in, Marshak out, Eddington start) + the lateral mode,
       # 400 steps of dt = 2e-4 (t = 0.08)
       'atm': dict(kap='100.0', T=0.08, dtc=2.0e-4*16.0, atm=True),
       'thin': dict(kap='1.0e-3', T=0.6, dtc=0.025),
       'mid': dict(kap='1.0', T=1.0, dtc=0.025)}
TH0, TH1, PH1 = 1.1707963267948966, 1.9707963267948966, 0.8


def cfg(case):
    p = case.split('_')
    return p[0], REG[p[1]], (len(p) > 2 and p[2] == 'vc')


def mk(tag, case, s):
    ax, rg, vc = cfg(case)
    n1 = 16*s
    n2 = 16*s if ax == 'th' else 4
    n3 = 16*s if ax == 'ph' else 4
    npr = 4 if s >= 8 else 1
    d = os.path.join(so.W, 'cpu', 'lateral', tag,
                     '%s_%d%s' % (case, s, os.environ.get('SO_SUF', '')))
    os.makedirs(d, exist_ok=True)
    mesh = {'nghost': 2, 'nx1': n1, 'x1min': 1.0, 'x1max': 3.0,
            'ix1_bc': 'outflow', 'ox1_bc': 'outflow',
            'nx2': n2, 'x2min': TH0, 'x2max': TH1, 'ix2_bc': 'reflect',
            'ox2_bc': 'reflect', 'nx3': n3, 'x3min': 0.0, 'x3max': PH1,
            'ix3_bc': 'periodic', 'ox3_bc': 'periodic', 'use_spherical_polar': 'true'}
    mb = {'nx1': n1, 'nx2': n2//npr if ax == 'th' else n2,
          'nx3': n3//npr if ax == 'ph' else n3}
    time = {'evolution': 'dynamic', 'integrator': 'rk2', 'cfl_number': 0.3,
            'nlim': 100000, 'tlim': rg['T'], 'ndiag': 100000}
    # dt = implicit_cfl * (smallest physical width)/c; the smallest width is the
    # lateral one at r = 1: 0.8/(16 s) (theta) or 0.8 sin(0.4+..)/(16 s) (phi)
    dxmin = 0.8/(16*s) if ax == 'th' else 0.8*np.sin(TH0)/(16*s)
    cfl = rg['dtc']/s/dxmin
    if os.environ.get('SO_DT'):
        cfl = float(os.environ['SO_DT'])/dxmin   # a FIXED dt for every s
    atm = rg.get('atm', False)
    rad = {'c_light': 1.0, 'chat_over_c': 1.0,
           'closure': 'vet_col' if vc else 'eddington', 'e_floor': '1.0e-30',
           'opacity': 'const', 'kappa_p': 0.0, 'kappa_f': 0.0, 'kappa_s': rg['kap'],
           'arad': 1.0, 'coupling': 'true', 'gas_feedback': 'false',
           'transport': 'implicit', 'time_scheme': 'hesdirk2',
           'implicit_cfl': '%.12g' % cfl, 'implicit_solver': 'bicgstab',
           'implicit_maxit': 200, 'implicit_precond': 'line',
           'implicit_tol': '1.0e-11', 'implicit_lin_tol': '1.0e-12',
           'implicit_bc_x1min': 'reflect', 'implicit_bc_x1max': 'reflect'}
    if atm:
        mesh['x1max'] = 5.0
        cfl = 100.0*rg['dtc']/s/dxmin
        rad.update({'c_light': 100.0, 'implicit_cfl': '%.12g' % cfl,
                    'implicit_bc_x1min': 'flux', 'implicit_flux_x1min': 1.0,
                    'implicit_bc_x1max': 'marshak', 'implicit_ebath_x1max': 0.0,
                    'marshak_q': 0.5})
    for kv in filter(None, os.environ.get('SO_RAD', '').split(',')):
        k, v = kv.split('=')
        rad[k] = v
    prob = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_shell', 'gas_rho': 1.0,
            'gas_temp': '1.0e-16', 'e_out': 1.0, 'lat_amp': 0.5,
            'lat_l2': 1.0 if ax == 'th' else 0.0, 'lat_l3': 1.0 if ax == 'ph' else 0.0}
    if atm:
        prob = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_atm', 'atm_rho0': 1.0,
                'atm_rho_n': 2.0, 'atm_temp': '1.0e-3', 'atm_init': 'eddington',
                'e_out': '1.0e-6', 'user_srcs': 'true', 'lat_amp': 0.5,
                'lat_l2': prob['lat_l2'], 'lat_l3': prob['lat_l3']}
    outs = [{'file_type': 'bin', 'variable': 'm1', 'dt': rg['T']}]
    inp = os.path.join(d, 'in.athinput')
    so.mkinput(inp, mesh, mb, time, rad, prob, outs)
    return d, inp, npr


def restrict(a, ax, s, r1=3.0):
    """2:1 volume restriction in r and in the refined lateral axis (uniform r)"""
    rf = np.linspace(1.0, r1, 16*s + 1)
    a = so.restrict_r(a, rf)
    if ax == 'th':
        tf = np.linspace(TH0, TH1, 16*s + 1)
        w = -np.diff(np.cos(tf))[None, :, None]
        b = a*w
        return (b[:, 0::2, :] + b[:, 1::2, :])/(w[:, 0::2, :] + w[:, 1::2, :])
    return 0.5*(a[0::2] + a[1::2])


def vol(ax, s, r1=3.0):
    rf = np.linspace(1.0, r1, 16*s + 1)
    tf = np.linspace(TH0, TH1, (16*s if ax == 'th' else 4) + 1)
    return (np.diff(rf**3)[None, None, :]*(-np.diff(np.cos(tf)))[None, :, None])


def main():
    act, tag, case = sys.argv[1:4]
    ss = [int(x) for x in sys.argv[4:]] or [1, 2, 4, 8]
    ax = cfg(case)[0]
    if act == 'run':
        b = os.path.join(so.W, 'bin', 'athena_%s_none_cpu' % tag)
        for s in ss:
            d, inp, npr = mk(tag, case, s)
            print(case, s, 'rc', so.run(b, inp, d, npr), flush=True)
        return
    E = {}
    for s in ss:
        d = mk(tag, case, s)[0]
        fs = sorted(f for f in os.listdir(os.path.join(d, 'bin')) if f.endswith('.bin'))
        E[s] = so.binall(os.path.join(d, 'bin', fs[-1]))
    errs = []
    for a, b in zip(ss[:-1], ss[1:]):
        ec = E[a]['m1_e']
        r1 = 5.0 if cfg(case)[1].get('atm') else 3.0
        ef = restrict(E[b]['m1_e'], ax, b, r1)
        w = vol(ax, a, r1)*np.ones_like(ec)
        mean = np.sum(w*ec)/np.sum(w)
        errs.append(np.sum(w*np.abs(ec - ef))/np.sum(w*np.abs(ec - mean)))
    print('%s %s t=%.6g' % (tag, case, E[ss[-1]]['time']))
    print('  L1 self-diff', ' '.join('%.3e' % e for e in errs))
    print('  orders      ', ' '.join('%.2f' % p for p in so.order(errs)))


if __name__ == '__main__':
    main()
