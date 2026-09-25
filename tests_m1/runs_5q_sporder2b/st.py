#!/usr/bin/env python3
"""runs_5q_sporder2b: COMBINED space-time convergence on the sp wedge.

dx and dt are refined together (implicit_cfl, or the hydro CFL, fixed), so the observed
order is min(space, time).  Self-convergence: e_n = L1_V(q_n - R q_2n)/L1_V(|q_n - <q_n>|)
per variable (E, F1 [, F2], and with gas T = (gamma-1) e/rho, the momentum and the total
energy), R the volume-weighted 2:1 restriction.

usage: st.py run  <case> <arm> [levels ...]     (runs the levels in parallel)
       st.py eval <case> <arm> [levels ...]
       st.py list
radial cases: levels n = 32 .. 512 (nx1 = n, 4 x 4 angular cells); lateral cases: levels
s = 1 2 4 8 (nx1 = n_lat = 16 s).  Arms: ARMS below (binary tag + <rad_m1> overrides).
"""
import os
import subprocess
import sys

import numpy as np

import solib as so

TH0, TH1, PH1 = 1.1707963267948966, 1.9707963267948966, 0.8
GM1 = 2.0/3.0

# arm: (binary tag, <rad_m1> overrides).  OLDK: every m1-sp-order2b key at its old value
OLDK = {'implicit_marshak_face': 'cell', 'vet_col_order2': 'false',
        'vet_col_fk_min': '0.3333333333333333', 'time2_vet_col': 'lag',
        'vet_col_reflect_top': 'false', 'vet_col_surface_face': 'false',
        'time2_vstage': 'false'}
ARMS = {
    'ref': ('ref', {}),                                   # fe12a518 defaults
    'new': ('new', {}),                                   # this branch's defaults
    # the fe12a518 behaviour with the new binary (the ref binary lacks gas_teq/gas_vth)
    'old': ('new', dict(OLDK)),
    'old_novimp': ('new', dict(OLDK, implicit_vimp='false')),
    'lag': ('new', {'time2_vet_col': 'lag'}),
    'rebuild': ('new', {'time2_vet_col': 'rebuild'}),
    'extrap': ('new', {'time2_vet_col': 'extrap'}),
    'clamp': ('new', {'vet_col_fk_min': '0.3333333333333333'}),
    'vacuum': ('new', {'vet_col_reflect_top': 'false'}),
    'sqcell': ('new', {'vet_col_surface_face': 'false'}),
    'vold': ('new', {'time2_vstage': 'false'}),
    'novimp': ('new', {'implicit_vimp': 'false'}),
    'be': ('new', {'time_scheme': 'be'}),
    'nc8': ('new', {'vet_col_ncore': 8}),
}


def base_rad(kap_s, kap_p=0.0, cfl=0.5, vc=True):
    return {'c_light': 1.0, 'chat_over_c': 1.0,
            'closure': 'vet_col' if vc else 'eddington', 'e_floor': '1.0e-30',
            'opacity': 'const', 'kappa_p': kap_p, 'kappa_f': kap_p, 'kappa_s': kap_s,
            'arad': 1.0, 'coupling': 'true', 'gas_feedback': 'false',
            'transport': 'implicit', 'time_scheme': 'hesdirk2',
            'implicit_cfl': '%.12g' % cfl, 'implicit_solver': 'bicgstab',
            'implicit_maxit': 200, 'implicit_precond': 'line',
            'implicit_tol': '1.0e-12', 'implicit_lin_tol': '1.0e-13',
            'implicit_bc_x1min': 'reflect', 'implicit_bc_x1max': 'reflect'} | (
                {'vet_col_ncore': 64} if vc else {})


def shell_prob(**kw):
    p = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_shell', 'gas_rho': 1.0,
         'gas_temp': '1.0e-16', 'e_out': '1.0e-2', 'e_amp': 1.0, 'r0': 2.0,
         'width': 0.2}
    p.update(kw)
    return p


def atm_prob(lat_amp=0.0, l2=0.0, l3=0.0):
    p = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_atm', 'atm_rho0': 1.0,
         'atm_rho_n': 2.0, 'atm_temp': '1.0e-3', 'atm_init': 'eddington',
         'e_out': '1.0e-6', 'user_srcs': 'true'}
    if lat_amp:
        p.update({'lat_amp': lat_amp, 'lat_l2': l2, 'lat_l3': l3})
    return p


ATM_RAD = {'c_light': 100.0, 'implicit_bc_x1min': 'flux', 'implicit_flux_x1min': 1.0,
           'implicit_bc_x1max': 'marshak', 'implicit_ebath_x1max': 0.0,
           'marshak_q': 0.5, 'kappa_s': 100.0, 'implicit_tol': '1.0e-11',
           'implicit_lin_tol': '1.0e-12'}


# case: geometry ('r' radial | 'th' | 'ph'), poly (stretched r), r1, tlim, hydro walls
# reflect (moving gas), gas output, rad, prob
NT = 128      # the fixed radial grid of the time-refinement cases (suffix _T)


def case_cfg(case):
    # suffix _T: TIME refinement only, nx1 = NT fixed, implicit_cfl = 8 cfl_case/l and
    # the hydro cfl_number = 0.8/l for the levels l = 1 2 4 8 16
    if case.endswith('_T'):
        c = case_cfg(case[:-2])
        c['tref'] = True
        return c
    c = {'geo': 'r', 'poly': case.endswith('_s'), 'r1': 3.0, 'gas': False,
         'hwall': 'outflow', 'tref': False}
    base = case[:-2] if case.endswith('_s') or case.endswith('_u') else case
    if base == 'pulse':
        # the vet_col radiating pulse (radial.py midm): rho kappa_s = 1, Marshak top
        c.update(T=1.0, rad=dict(base_rad(1.0), implicit_bc_x1max='marshak',
                                 marshak_q=0.5, implicit_ebath_x1max=0.0),
                 prob=shell_prob(e_out='1.0e-8'))
    elif base == 'pulse_edd':
        c.update(T=1.0, rad=dict(base_rad(1.0, vc=False), implicit_bc_x1max='marshak',
                                 marshak_q=0.5, implicit_ebath_x1max=0.0),
                 prob=shell_prob(e_out='1.0e-8'))
    elif base == 'fs':
        # free-streaming shell (rho kappa_s = 1e-3), vet_col, reflecting walls
        c.update(T=0.6, rad=base_rad(1.0e-3), prob=shell_prob())
    elif base == 'fsm':
        # the same with a Marshak top (the shell leaves through the top)
        c.update(T=1.5, rad=dict(base_rad(1.0e-3), implicit_bc_x1max='marshak',
                                 marshak_q=0.5, implicit_ebath_x1max=0.0),
                 prob=shell_prob(e_out='1.0e-8'))
    elif base == 'rtop':
        # rho kappa_s = 1 pulse between REFLECTING walls under vet_col (the mirror top)
        c.update(T=1.0, rad=base_rad(1.0), prob=shell_prob())
    elif base == 'atm':
        # the T-S4 vet_col atmosphere TRANSIENT (sph_atm, r = 1..5, c = 100, Eddington
        # start), t = 0.08; dt = 0.64 dx/c (fixed implicit_cfl)
        c.update(T=0.08, r1=5.0, rad=dict(base_rad(100.0), **ATM_RAD),
                 prob=atm_prob())
        c['rad']['implicit_cfl'] = '0.64'
    elif base.startswith('stiff'):
        # STIFF coupled gas + radiation: grey absorption rho kappa = K, the gas at
        # radiative equilibrium with the pulse (gas_teq), static (dbg_gas_force =
        # false), energy exchange on; c dt rho kappa = cfl dx K.  stiffK[_cC][_edd]
        tok = base.split('_')
        K = float(tok[0][5:])
        cfl = 0.5
        vc = True
        for t in tok[1:]:
            if t.startswith('c'):
                cfl = float(t[1:])
            if t == 'edd':
                vc = False
        # arad = 1e8, rho = 270: T ~ 0.01 (c_s ~ 0.1 c), gas and radiation heat
        # capacities comparable (1.5 rho ~ 4 a T^3)
        rad = dict(base_rad(0.0, kap_p=K/270.0, cfl=cfl, vc=vc), gas_feedback='true',
                   dbg_gas_force='false', arad='1.0e8')
        c.update(T=1.0, rad=rad, gas=True, hwall='reflect',
                 prob=shell_prob(e_out='0.1', gas_teq='true', gas_rho=270.0))
    elif base.startswith('mvr'):
        # moving gas, radial: v_r = 0.01 sin, rho kappa_s = 1000 + rho kappa_p = 10,
        # radiation force and heating on (vimp default on), t = 10 (mvrs: no kappa_p)
        # (radiation energy ~1e-12, the gas thermal 1.5e-3: the gas carries E)
        kp = 0.0 if 'mvrs' in base else 10.0
        rad = dict(base_rad(1000.0, kap_p=kp, cfl=50.0, vc=('edd' not in base)),
                   gas_feedback='true')
        c.update(T=10.0, rad=rad,
                 gas=True, hwall='reflect',
                 prob=shell_prob(gas_vr=0.01, gas_temp='1.0e-3', e_out='1.0e-12',
                                 e_amp='1.0e-12'))
    elif base.startswith('mvp'):
        # passive radiation (E ~ 1, gas_feedback = false) advected by the v_r = 0.01
        # sin gas flow through rho kappa_s = 1000 (vimp needs gas_feedback: off)
        c.update(T=10.0, rad=dict(base_rad(1000.0, cfl=50.0, vc=('edd' not in base))),
                 gas=True, hwall='reflect',
                 prob=shell_prob(gas_vr=0.01, gas_temp='1.0e-3', e_out='1.0'))
        if 'mvpt' in base:
            # the same at E ~ 1e-12 (the scale of mvr): tolerance / floor check
            c['prob'].update(e_out='1.0e-12', e_amp='1.0e-12')
    elif base.startswith('rsw'):
        # radiation-dominated moving slab: arad = 3e10, E = 0.03 (1 + 0.1 Gaussian) at
        # equilibrium with the gas (T ~ 1e-3): P_rad/P_gas ~ 10, E/(rho c^2) = 0.03; the
        # radiation-acoustic pulse (c_s ~ 0.11 c) moves the gas, rho kappa_s = 1000
        # (v tau/c ~ 2), rho kappa_p = 1, t = 2, dt = 0.5 dx/c
        c.update(T=2.0, rad=dict(base_rad(1000.0, kap_p=1.0, cfl=0.5,
                                          vc=('edd' not in base)),
                                 gas_feedback='true', arad='3.0e10'),
                 gas=True, hwall='reflect',
                 prob=shell_prob(e_out='0.03', e_amp='0.003', gas_teq='true'))
    elif base.startswith('lat'):
        # lateral modes (lateral.py): lat_<th|ph>_<atm|mid|mvt>
        _, ax, reg = base.split('_')
        c['geo'] = ax
        l2, l3 = (1.0, 0.0) if ax == 'th' else (0.0, 1.0)
        if reg == 'atm':
            c.update(T=0.08, r1=5.0, rad=dict(base_rad(100.0), **ATM_RAD),
                     prob=atm_prob(0.5, l2, l3), dtc=2.0e-4*16.0)
        elif reg == 'mid':
            c.update(T=1.0, rad=base_rad(1.0),
                     prob=shell_prob(e_out=1.0, e_amp=0.0, lat_amp=0.5, lat_l2=l2,
                                     lat_l3=l3), dtc=0.025)
        elif reg == 'mvt':
            # a lateral flow v_theta = 0.01 sin through a thick scattering shell with a
            # lateral E mode (E ~ 1e-12 = a T_gas^4: the gas carries the radiation, as
            # mvr), radiation force and heating on (vimp), t = 10
            c.update(T=10.0, rad=dict(base_rad(1000.0, kap_p=10.0),
                                      gas_feedback='true'),
                     gas=True, hwall='reflect', dtc=16.0,
                     prob=shell_prob(e_out='1.0e-12', e_amp=0.0, lat_amp=0.5,
                                     lat_l2=l2, lat_l3=l3, gas_vth=0.01,
                                     gas_temp='1.0e-3'))
    else:
        raise SystemExit('unknown case ' + case)
    return c


def levels(case):
    c = case_cfg(case)
    if c['tref']:
        return [1, 2, 4, 8, 16]
    return [1, 2, 4, 8] if c['geo'] != 'r' else [32, 64, 128, 256, 512]


def mk(case, arm, n):
    c = case_cfg(case)
    tag, over = ARMS[arm]
    geo = c['geo']
    d = os.path.join(so.W, 'st', case, '%s_%d' % (arm, n))
    os.makedirs(d, exist_ok=True)
    if c['tref']:
        n1, n2, n3, npr = NT, 4, 4, 1
    elif geo == 'r':
        n1, n2, n3, npr = n, 4, 4, 1
    else:
        n1 = 16*n
        n2 = 16*n if geo == 'th' else 4
        n3 = 16*n if geo == 'ph' else 4
        npr = 4 if n >= 8 else 1
    mesh = {'nghost': 2, 'nx1': n1, 'x1min': 1.0, 'x1max': c['r1'],
            'ix1_bc': c['hwall'], 'ox1_bc': c['hwall'],
            'nx2': n2, 'x2min': TH0, 'x2max': TH1, 'ix2_bc': 'reflect',
            'ox2_bc': 'reflect', 'nx3': n3, 'x3min': 0.0, 'x3max': PH1,
            'ix3_bc': 'periodic', 'ox3_bc': 'periodic', 'use_spherical_polar': 'true'}
    if c['poly']:
        mesh['use_grid_stretch_r_poly'] = 'true'
        mesh.update(so.POLY)
    mb = {'nx1': n1, 'nx2': n2//npr if geo == 'th' else n2,
          'nx3': n3//npr if geo == 'ph' else n3}
    time = {'evolution': 'dynamic', 'integrator': 'rk2', 'cfl_number': 0.3,
            'nlim': 1000000, 'tlim': c['T'], 'ndiag': 1000000}
    rad = dict(c['rad'])
    if geo != 'r':
        # dt = dtc/s: implicit_cfl from the smallest (lateral, r = 1) width
        dxmin = 0.8/n2 if geo == 'th' else 0.8*np.sin(TH0)/n3
        cl = float(rad['c_light'])
        rad['implicit_cfl'] = '%.12g' % (cl*c['dtc']/n/dxmin)
    if c['tref']:
        # both limiters scale with 1/l, so dt halves exactly from level to level
        rad['implicit_cfl'] = '%.12g' % (8.0*float(rad['implicit_cfl'])/n)
        time['cfl_number'] = '%.12g' % (0.8/n)
    rad.update(over)
    outs = []
    if geo == 'r':
        o = {'file_type': 'tab', 'variable': 'm1', 'data_format': '%24.16e',
             'slice_x2': 1.5207963267948967, 'slice_x3': 0.05, 'dt': c['T']}
        outs.append(o)
        if c['gas']:
            outs.append(dict(o, variable='hydro_w'))
    else:
        outs.append({'file_type': 'bin', 'variable': 'm1', 'dt': c['T']})
        if c['gas']:
            outs.append({'file_type': 'bin', 'variable': 'hydro_w', 'dt': c['T']})
    inp = os.path.join(d, 'in.athinput')
    so.mkinput(inp, mesh, mb, time, rad, c['prob'], outs)
    b = os.path.join(so.W, 'bin', 'athena_%s_none_cpu' % tag)
    return d, inp, b, npr


def launch(case, arm, ns):
    procs = []
    for n in ns:
        d, inp, b, npr = mk(case, arm, n)
        for f in os.listdir(d):
            if f in ('tab', 'bin'):
                subprocess.call(['rm', '-rf', os.path.join(d, f)])
        cmd = ['mpirun', '-np', str(npr), '--oversubscribe', '--bind-to', 'none', b,
               '-i', inp, '-d', d]
        lg = open(os.path.join(d, 'log'), 'w')
        env = dict(os.environ, OMP_NUM_THREADS='1')
        procs.append((n, subprocess.Popen(['nice'] + cmd, stdout=lg,
                                          stderr=subprocess.STDOUT, env=env)))
    for n, p in procs:
        print(case, arm, n, 'rc', p.wait(), flush=True)


def load(case, arm, n):
    c = case_cfg(case)
    d = mk(case, arm, n)[0]
    out = {}
    if c['geo'] == 'r':
        for v in (['m1'] + (['hydro_w'] if c['gas'] else [])):
            fs = sorted(f for f in os.listdir(os.path.join(d, 'tab'))
                        if f.endswith('.tab') and ('.%s.' % v) in f)
            t = so.tab(os.path.join(d, 'tab', fs[-1]))
            for k, a in t.items():
                out[k] = a
    else:
        for v in (['m1'] + (['hydro_w'] if c['gas'] else [])):
            fs = sorted(f for f in os.listdir(os.path.join(d, 'bin'))
                        if f.endswith('.bin') and ('.%s.' % v) in f)
            t = so.binall(os.path.join(d, 'bin', fs[-1]))
            for k, a in t.items():
                out[k] = a
    q = {'E': out['m1_e'], 'F1': out['m1_f1']}
    if c['geo'] != 'r':
        q['F%d' % (2 if c['geo'] == 'th' else 3)] = out['m1_f%d' % (
            2 if c['geo'] == 'th' else 3)]
    if c['gas']:
        rho = out['dens']
        q['T'] = GM1*out['eint']/rho
        q['M1'] = rho*out['velx']
        ek = 0.5*rho*(out['velx']**2 + out['vely']**2 + out['velz']**2)
        q['Eg'] = out['eint'] + ek
        if c['geo'] == 'th':
            q['M2'] = rho*out['vely']
    q['time'] = out['time']
    return q


def restrict(a, c, n):
    """2:1 volume restriction of level n+1 data to level n (r, and theta/phi)"""
    geo = c['geo']
    if c['tref']:
        return a
    nf = (2*n if geo == 'r' else 32*n)
    rf = so.rfaces(1.0, c['r1'], nf, c['poly'])
    a = so.restrict_r(a, rf)
    if geo == 'th':
        tf = np.linspace(TH0, TH1, 32*n + 1)
        w = -np.diff(np.cos(tf))[None, :, None]
        b = a*w
        return (b[:, 0::2, :] + b[:, 1::2, :])/(w[:, 0::2, :] + w[:, 1::2, :])
    if geo == 'ph':
        return 0.5*(a[0::2] + a[1::2])
    return a


def vol(c, n, shape):
    geo = c['geo']
    if c['tref']:
        n = NT
    rf = so.rfaces(1.0, c['r1'], n if geo == 'r' else 16*n, c['poly'])
    w = np.diff(rf**3)
    if geo == 'th':
        tf = np.linspace(TH0, TH1, 16*n + 1)
        w = w[None, None, :]*(-np.diff(np.cos(tf)))[None, :, None]
    return w*np.ones(shape)


def evaluate(case, arm, ns):
    c = case_cfg(case)
    Q = {n: load(case, arm, n) for n in ns}
    keys = [k for k in Q[ns[0]] if k != 'time']
    res = {}
    for k in keys:
        errs = []
        for a, b in zip(ns[:-1], ns[1:]):
            ec = Q[a][k]
            ef = restrict(Q[b][k], c, a)
            w = vol(c, a, ec.shape)
            mean = np.sum(w*ec)/np.sum(w)
            errs.append(np.sum(w*np.abs(ec - ef))/np.sum(w*np.abs(ec - mean)))
        res[k] = errs
    return res, Q[ns[-1]]['time']


def main():
    act = sys.argv[1]
    if act == 'list':
        print(' '.join(sorted(ARMS)))
        return
    case, arm = sys.argv[2:4]
    ns = [int(x) for x in sys.argv[4:]] or levels(case)
    if act == 'run':
        launch(case, arm, ns)
        return
    res, t = evaluate(case, arm, ns)
    for k, e in res.items():
        print('%-14s %-10s %-3s t=%-7.4g L1 %s | orders %s' % (
            case, arm, k, t, ' '.join('%.2e' % x for x in e),
            ' '.join('%5.2f' % p for p in so.order(e))))


if __name__ == '__main__':
    main()
