#!/usr/bin/env python3
"""cart.py: Cartesian order suite for implicit M1 + vet_sc + halfrange (hr2prod-1110).

2-D isotropic Gaussian E pulse (rad_m1_beam m1_test = thick_pulse) in a static pure-scattering
medium, box [-1,1]^2 (x1 outflow = vet_sc sweep axis, x2 periodic), one MeshBlock, production
Cartesian rad_m1 keys of the beam deck cb_FIN3_r4 (blend_1022/p2), k1_redist 0.
  cases: pk1 (kappa_s 1.28: streaming), pk12 (12.8: transition, 0 < w < 1 band), pk128 (128: diffusive)
  modes: C (nx1 = nx2 = n, c dt/dx = NU_C), T (n = 128 fixed, dt = DT_T(case)/l)
usage: cart.py deck CASE ARM MODE LEV  -> prints the deck path
       cart.py eval CASE ARM MODE LEV...
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
H = os.environ.get('CART_H', os.path.join(os.getcwd(), 'cart_runs'))   # run tree
SRC = os.path.join(HERE, 'cb_FIN3_r4.athinput')   # rad_m1 keys of the Cartesian beam deck
sys.path.insert(0, HERE)
CASES = {'pk1': (1.28, 0.3), 'pk12': (12.8, 0.3), 'pk128': (128.0, 0.3)}
NU_C = 0.5
NT = 128
NU_T = 4.0   # c dt/dx of the coarsest T level at n = 128

NOCJ = {'implicit_thin_corr': 'none', 'implicit_thin_recon': 'dc'}
PLM = dict(NOCJ, implicit_hr_recon='plm',
           implicit_hr_recon_lim='smooth')
MOD = {
    'base': {},
    'cen': {'implicit_flux': 'central', 'implicit_flux_beam': 'closure',
            'implicit_thin_switch': 'none', 'implicit_thin_corr': 'none',
            'implicit_flux_faces': 'x1', 'implicit_blend': 'tau_f'},
    'nocj': NOCJ,
    'plm': PLM,
    'gmean': {'implicit_thin_gface': 'mean'},
    'cm': {'implicit_hr_cmem': 'true'},
    'ber': {'implicit_flux': 'berthon', 'implicit_thin_switch': 'none'},
    'nocorr': {'implicit_thin_corr': 'none'},
    'norecon': {'implicit_thin_recon': 'dc'},
    'frz0': {'implicit_thin_freeze': '0'},
    'xsc': {'time2_vet_sc': 'lag'},
    'edd': {'closure': 'eddington', 'implicit_hr_model': 'true'},
    'plmh': dict(PLM, implicit_hr_recon_q='h', implicit_hr_face='cell'),
    'hrcell': {'implicit_hr_face': 'cell'},
    'gjump': {'implicit_thin_gface': 'jump'},
    'dtc': {'implicit_blend_dtc': '1.0e-6'},
    'nofse': {'implicit_thin_switch': 'none'},
    'w1': {'implicit_thin_h0': '-2.0', 'implicit_thin_h1': '-1.0'},
    'wsx0': {'implicit_hr_wscale_x1': '0'},
    'wsl0': {'implicit_hr_wscale_lat': '0'},
    'vlag': {'time2_vet_sc': 'lag'},
}


def parse(path):
    blocks, order, b = {}, [], None
    for ln in open(path):
        ln = ln.split('#')[0].rstrip()
        if ln.startswith('<'):
            b = ln.strip('<> ')
            blocks[b] = {}
            order.append(b)
        elif '=' in ln and b:
            k, v = [x.strip() for x in ln.split('=', 1)]
            blocks[b][k] = v
    return blocks


def write(path, blocks):
    with open(path, 'w') as f:
        for b, d in blocks.items():
            f.write('<%s>\n' % b)
            for k, v in d.items():
                f.write('%-34s = %s\n' % (k, v))
            f.write('\n')


def rundir(case, arm, mode, lev):
    return os.path.join(H, 'runs', case, arm, '%s_%d' % (mode, lev))


def deck(case, arm, mode, lev):
    kap, tlim = CASES[case]
    B = parse(SRC)
    n = lev if mode == 'C' else NT
    dx = 2.0/n
    c = float(B['rad_m1']['c_light'])
    m = B['mesh']
    m.update({'nx1': n, 'x1min': -1.0, 'x1max': 1.0, 'ix1_bc': 'outflow', 'ox1_bc': 'outflow',
              'nx2': n, 'x2min': -1.0, 'x2max': 1.0, 'ix2_bc': 'periodic', 'ox2_bc': 'periodic'})
    B['meshblock'] = {'nx1': n, 'nx2': n, 'nx3': 1}
    if os.environ.get('CART_1D'):   # plane-parallel: pulse_1d, nx2 = 4 at every level
        m['nx2'] = 4
        B['meshblock']['nx2'] = 4
        if os.environ.get('CART_1D') == '1t':   # a TRUE 1-D mesh (nx2 = 1: the x1 hr model of RF0)
            m['nx2'] = 1
            B['meshblock']['nx2'] = 1
    if os.environ.get('CART_3D'):   # 3-D box, x3 periodic, refined with x1/x2
        m.update({'nx3': n, 'x3min': -1.0, 'x3max': 1.0, 'ix3_bc': 'periodic',
                  'ox3_bc': 'periodic'})
        B['meshblock']['nx3'] = n
    t = B['time']
    t.update({'dt_retry': 'false', 'nlim': '-1', 'tlim': '%.15g' % (tlim/c), 'ndiag': '100000',
              'cfl_number': '0.95'})
    for k in [k for k in t if k.startswith('dt_retry_')]:
        del t[k]
    dt = NU_C*dx/c if mode == 'C' else NU_T*dx/c/lev
    t['dt_fixed'] = '%.15g' % dt
    r = B['rad_m1']
    r.update({'kappa_p': '0.0', 'kappa_e': '0.0', 'kappa_f': '0.0', 'kappa_s': '%g' % kap,
              'time2_k1_redist': '0', 'implicit_cfl': '1.0e6', 'implicit_tol': '1.0e-10',
              'implicit_tol_rel': '0',
              # rotation off (user rule 11-10), the 512-ray vet_sc set (no 4-ray beam ring)
              'vet_gd_rotate_every': '0', 'vet_gd_rotate_dt': '0', 'vet_rotate_avg': '1',
              'vet_rotate_avg_ema': '0', 'vet_nmu': '8', 'vet_nphi': '32'})
    r.pop('vet_ring', None)
    for a in arm.split('~'):
        for k, v in MOD[a].items():
            r[k] = v
    # no boundary beams (the cb deck's implicit_beam*: two x1min beams, NOT part of the pulse test)
    for k in [k for k in r if k.startswith('implicit_beam')]:
        del r[k]
    if r.get('closure') == 'eddington':   # fixed tensor: no formal solution
        for k in [k for k in r if k.startswith('vet_')]:
            del r[k]
    B['problem'] = {'pgen_name': 'rad_m1_beam', 'm1_test': 'thick_pulse', 'gas_rho': '1.0',
                    'gas_eint': '0.01', 'pulse_amp': '1.0', 'pulse_width': os.environ.get('CART_W', '0.25'),
                    'pulse_x0': '0.0', 'pulse_y0': '0.0', 'pulse_bg': '1.0e-3'}
    if os.environ.get('CART_1D'):
        B['problem']['pulse_1d'] = 'true'
    for k in [k for k in B if k.startswith('output')]:
        del B[k]
    B['output1'] = {'file_type': 'bin', 'variable': 'm1', 'dt': '%.15g' % (tlim/c)}
    d = rundir(case, arm, mode, lev)
    os.makedirs(d, exist_ok=True)
    p = d + '/in.athinput'
    write(p, B)
    return p


def load(d):
    import bin_convert as bc
    fs = sorted(f for f in os.listdir(d + '/bin') if f.endswith('.bin'))
    b = bc.read_binary(d + '/bin/' + fs[-1])
    out = {}
    for v in ('m1_e', 'm1_f1', 'm1_f2'):
        a = np.asarray(b['mb_data'][v][0])
        out[v] = a.reshape(-1, a.shape[-2], a.shape[-1]) if a.ndim >= 3 and a.shape[-3] > 1 else a.reshape(a.shape[-2], a.shape[-1])
    return out, b['time']


def evaluate(case, arm, mode, levs):
    Q = {}
    for l in levs:
        Q[l], tt = load(rundir(case, arm, mode, l))
    for v in ('m1_e', 'm1_f1', 'm1_f2'):
        e = []
        for a, b in zip(levs[:-1], levs[1:]):
            qa, qb = Q[a][v], Q[b][v]
            if mode == 'C':
                if qb.ndim == 1 or (qb.ndim == 2 and qb.shape[0] == 1 and qa.ndim == 2):
                    qb = 0.5*(qb[..., 0::2] + qb[..., 1::2])
                elif qb.ndim == 2 and qb.shape[0] == qa.shape[0]:   # plane-parallel (CART_1D)
                    qb = 0.5*(qb[:, 0::2] + qb[:, 1::2])
                elif qb.ndim == 3:
                    qb = 0.125*sum(qb[a::2, b::2, c::2] for a in (0, 1) for b in (0, 1) for c in (0, 1))
                else:
                    qb = 0.25*(qb[0::2, 0::2] + qb[1::2, 0::2] + qb[0::2, 1::2] + qb[1::2, 1::2])
            e.append(np.abs(qa - qb).mean()/np.abs(qa - qa.mean()).mean())
        p = [np.log2(e[i]/e[i+1]) for i in range(len(e) - 1)]
        print('%-6s %-24s %s %-6s d %s | p %s' % (case, arm, mode, v, ' '.join('%.3e' % x for x in e),
                                                   ' '.join('%5.2f' % x for x in p)))


if __name__ == '__main__':
    act, case, arm, mode = sys.argv[1:5]
    levs = [int(x) for x in sys.argv[5:]]
    if act == 'deck':
        print(deck(case, arm, mode, levs[0]))
    else:
        evaluate(case, arm, mode, levs)
