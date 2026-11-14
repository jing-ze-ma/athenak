#!/usr/bin/env python3
"""unify-1114: Cartesian (plane-parallel) twin of runs_5q_sporder2b/st.py.

The same rad_m1_beam sph_shell / sph_atm problems on a Cartesian slab x1 = 1..r1 (4 x 4
periodic transverse cells), uniform cells.  Self-convergence e_n = L1(q_n - R q_2n)/
L1(|q_n - <q_n>|), R the 2:1 average; _T cases refine dt only (nx1 = NT).
usage: cst.py run <case> <arm> [levels]  |  cst.py eval <case> <arm> [levels]
arms: BIN[tag] binary + <rad_m1> overrides (ARMS).
"""
import os
import subprocess
import sys

import numpy as np

import solib as so

BIN = {'u0': '/viper/ptmp2/jinma/builds/bin/athena_none_cpu_b2d852bb_u0'}
BIN['new'] = '/viper/ptmp2/jinma/builds/bin/athena_none_cpu_a49481fa_u2b'
BIN.update({k: v for k, v in (x.split('=') for x in os.environ.get('CST_BINS', '').split()
                              if '=' in x)})
OLDK = {'vet_col_order2': 'false', 'vet_col_fk_min': '0.3333333333333333',
        'time2_vet_col': 'lag'}
ARMS = {
    'old': ('u0', {}),                       # prod-1118 Cartesian defaults
    'new': ('new', {}),                      # the branch's defaults
    'newold': ('new', dict(OLDK)),           # new binary, old Cartesian values
    'lag': ('new', {'time2_vet_col': 'lag'}),
    'clamp': ('new', {'vet_col_fk_min': '0.3333333333333333'}),
    'o1': ('new', {'vet_col_order2': 'false'}),
    'mlin': ('new', {'implicit_marshak_face': 'linear'}),
    'mlsf': ('new', {'implicit_marshak_face': 'linear', 'vet_col_surface_face': 'true'}),
}
NT = 128


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


def case_cfg(case):
    if case.endswith('_T'):
        c = case_cfg(case[:-2])
        c['tref'] = True
        return c
    c = {'r1': 3.0, 'tref': False}
    mq = dict(implicit_bc_x1max='marshak', marshak_q=0.5, implicit_ebath_x1max=0.0)
    if case == 'pulse':      # radiating pulse, rho kappa_s = 1, Marshak top
        c.update(T=1.0, rad=dict(base_rad(1.0), **mq), prob=shell_prob(e_out='1.0e-8'))
    elif case == 'pulse_edd':
        c.update(T=1.0, rad=dict(base_rad(1.0, vc=False), **mq),
                 prob=shell_prob(e_out='1.0e-8'))
    elif case == 'fsm':      # free-streaming shell leaving through a Marshak top
        c.update(T=1.5, rad=dict(base_rad(1.0e-3), **mq), prob=shell_prob(e_out='1.0e-8'))
    elif case == 'thick':    # rho kappa_s = 20 diffusing pulse, Marshak top
        c.update(T=1.0, rad=dict(base_rad(20.0), **mq), prob=shell_prob(e_out='1.0e-8'))
    else:
        raise SystemExit('unknown case ' + case)
    return c


def levels(case):
    return [1, 2, 4, 8, 16] if case_cfg(case)['tref'] else [32, 64, 128, 256, 512]


def mk(case, arm, n):
    c = case_cfg(case)
    tag, over = ARMS[arm]
    d = os.path.join(so.W, 'st', case, '%s_%d' % (arm, n))
    os.makedirs(d, exist_ok=True)
    n1 = NT if c['tref'] else n
    mesh = {'nghost': 2, 'nx1': n1, 'x1min': 1.0, 'x1max': c['r1'],
            'ix1_bc': 'outflow', 'ox1_bc': 'outflow',
            'nx2': 4, 'x2min': -0.5, 'x2max': 0.5, 'ix2_bc': 'periodic',
            'ox2_bc': 'periodic', 'nx3': 4, 'x3min': -0.5, 'x3max': 0.5,
            'ix3_bc': 'periodic', 'ox3_bc': 'periodic'}
    mb = {'nx1': n1, 'nx2': 4, 'nx3': 4}
    time = {'evolution': 'dynamic', 'integrator': 'rk2', 'cfl_number': 0.3,
            'nlim': 1000000, 'tlim': c['T'], 'ndiag': 1000000}
    rad = dict(c['rad'])
    if c['tref']:
        rad['implicit_cfl'] = '%.12g' % (8.0*float(rad['implicit_cfl'])/n)
        time['cfl_number'] = '%.12g' % (0.8/n)
    rad.update(over)
    outs = [{'file_type': 'tab', 'variable': 'm1', 'data_format': '%24.16e',
             'slice_x2': 0.1, 'slice_x3': 0.1, 'dt': c['T']}]
    inp = os.path.join(d, 'in.athinput')
    so.mkinput(inp, mesh, mb, time, rad, c['prob'], outs)
    return d, inp, BIN[tag]


def launch(case, arm, ns):
    procs = []
    for n in ns:
        d, inp, b = mk(case, arm, n)
        cmd = ['nice', b, '-i', inp, '-d', d]
        lg = open(os.path.join(d, 'log'), 'w')
        env = dict(os.environ, OMP_NUM_THREADS='1')
        procs.append((n, d, subprocess.Popen(cmd, stdout=lg, stderr=subprocess.STDOUT,
                                             env=env)))
    for n, d, p in procs:
        rc = p.wait()
        lg = open(os.path.join(d, 'log')).read()
        print(case, arm, n, 'rc', rc, 'FATAL', lg.count('FATAL'),
              'NONCONV', lg.count('NON-CONVERGED'), flush=True)


def load(case, arm, n):
    d = mk(case, arm, n)[0]
    fs = sorted(f for f in os.listdir(os.path.join(d, 'tab'))
                if f.endswith('.tab') and '.m1.' in f)
    t = so.tab(os.path.join(d, 'tab', fs[-1]))
    return {'E': t['m1_e'], 'F1': t['m1_f1'], 'time': t['time']}


def evaluate(case, arm, ns):
    c = case_cfg(case)
    Q = {n: load(case, arm, n) for n in ns}
    res = {}
    for k in ('E', 'F1'):
        errs = []
        for a, b in zip(ns[:-1], ns[1:]):
            ec = Q[a][k]
            ef = Q[b][k] if c['tref'] else 0.5*(Q[b][k][0::2] + Q[b][k][1::2])
            errs.append(np.sum(np.abs(ec - ef))/np.sum(np.abs(ec - ec.mean())))
        res[k] = errs
    return res, Q[ns[-1]]['time']


def main():
    act = sys.argv[1]
    case, arm = sys.argv[2:4]
    ns = [int(x) for x in sys.argv[4:]] or levels(case)
    if act == 'run':
        launch(case, arm, ns)
        return
    res, t = evaluate(case, arm, ns)
    for k, e in res.items():
        print('%-10s %-8s %-3s t=%-7.4g L1 %s | orders %s' % (
            case, arm, k, t, ' '.join('%.2e' % x for x in e),
            ' '.join('%5.2f' % p for p in so.order(e))))


if __name__ == '__main__':
    main()
