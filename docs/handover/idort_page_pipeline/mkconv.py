#!/usr/bin/env python3
"""idort-1019: convergence table rows (conv_rows.json) from rwconv.json / spxconv.json / csconv.json / ladconv.json"""
import json, os
D = '/raven/ptmp/jinma/idort1020/page'
J = lambda f: json.load(open(os.path.join(D, f))) if os.path.exists(os.path.join(D, f)) else {}
rows = []
def st(p): return 'PASS' if p >= 1.8 else 'FLAG'
f2 = lambda v: ' '.join('%.2e' % x for x in v); o2 = lambda v: ' '.join('%.2f' % x for x in v)
for k, o in J('rwconv.json').items():
    arm, dr, P, t = k.split('|')
    if arm == 'TIGHT': continue
    for var in ('rho', 'E'):
        fit, fine, pairs = o[var]
        rows.append(dict(case='radwave %s (%s)' % (dr, 'PROD' if arm == 'PROD19' else 'v10 FIN3'), geom='Cartesian vet_sc', regime='P_rad/P_gas %s, tau %s' % ('3.3e-3' if float(P[1:]) < 1 else '3.33', t[1:]),
                         var=var, N='%s (dt ~ dx)' % ' '.join(map(str, o['N'])), err=f2(o['l1_' + var]) + ' (vs exact)', pairs=o2(pairs), fit='%.2f' % fit,
                         status=st(fine) if arm == 'PROD19' else 'INFO', arm=arm))
KD = {'X': 'space, fixed dt', 'C': 'dt ~ dx', 'L': 'lateral phi', 'T': 'dt / q at N 256'}
DESC = {'rsw': 'radiation-acoustic pulse, rho kappa_s 1000', 'rswt': 'radiation-acoustic pulse, rho kappa_s 10 (semi-thin)', 'rswv': 'rsw in a moving shell v_r 0.01 (suite deck)',
        'rswtv': 'rswt in a moving shell (suite deck)', 'rswv_F': 'rswv, user-approved deck (f0_init zero, outflow top)', 'rswtv_F': 'rswtv with the rswv deck fix (diagnostic)',
        'marsh_r1_u': 'spherical Marshak wave (front, not smooth)', 'patmT_L16_s': 'grey atmosphere relaxation (absorption)', 'rlat': 'lateral radiation-acoustic wave, kappa_s 1000',
        'rlatt': 'lateral wave, kappa_s 10', 'rlat_om5': 'rlat + solid-body rotation', 'rlatt_om5': 'rlatt + rotation', 'rswt_dt': 'rswt, dt refinement only'}
for k, r in J('spxconv.json').items():
    case, arm = k.split('|')
    if arm != 'P': continue
    base = case.replace('_C', '') if case.endswith('_C') else case
    for var in ('m1_e', 'dens'):
        if var not in r['err']: continue
        fine = r['order'][var][-1]
        s = st(fine)
        if case.startswith('marsh') or case in ('rswv', 'rswtv', 'rswv_C', 'rswtv_C'): s = 'INFO'
        if case == 'rswt_dt' and min(r['err'][var]) < 3e-8: s = 'PASS (floor)' if r['order'][var][0] >= 1.8 else s
        rows.append(dict(case='%s: %s' % (case, DESC.get(base, '')), geom='sp wedge vet_col + vet_gd', regime=KD[r['kind']], var='E' if var == 'm1_e' else 'rho',
                         N=' '.join(map(str, r['levs'])), err=f2(r['err'][var]) + ' (self-conv.)', pairs=o2(r['order'][var]), fit='%.2f' % r['fit'][var], status=s, arm='P'))
CL = {'gd': 'diffusion mode l=2, thick (kappa_s 1000), vet_col + vet_gd', 'gda': 'diffusion mode in an atmosphere (tau 20, Marshak top), vet_col + vet_gd', 'edd': 'diffusion mode, Eddington closure (harness keys only)'}
for k, r in J('csconv.json').items():
    mesh, case, tag, kind = k.split('|')
    if tag not in ('p19', 'r19'): continue
    if tag == 'r19' and case != 'edd': continue
    for cls in (('all', 'seam', 'vertex') if mesh == 'cs' else ('all',)):
        if cls not in r['err']: continue
        fine = r['order'][cls][-1]
        rows.append(dict(case='%s (%s)' % (CL[case], 'PROD' if tag == 'p19' else 'harness keys'), geom='cubed sphere' if mesh == 'cs' else 'sp wedge', regime='%s, %s cells' % ('space N, dt ~ dx' if kind == 'space' else 'time dt/q at N 16', cls),
                         var='E', N=' '.join(map(str, r['x'] + [r['x'][-1] * 2])) + (' (N panel)' if kind == 'space' else ' (q)'), err=f2(r['err'][cls]) + ' (self-conv.)', pairs=o2(r['order'][cls]),
                         fit='%.2f' % fine, status=st(fine) if tag == 'p19' else 'INFO', arm=tag))
for arm, d in (('p', J('ladconv.json').get('p', {})),):
    for k, r in d.items():
        case, a, kind, var = k.split('|')
        if var not in ('e', 'dens', 'f1') or a not in ('FIN', 'FINs'): continue
        fine = r['ord'][-1]
        rows.append(dict(case='ladder %s' % case, geom='Cartesian vet_sc', regime={'X': 'space, fixed dt', 'C': 'dt ~ dx', 'T': 'dt refinement'}[kind], var={'e': 'E', 'f1': 'F1', 'dens': 'rho'}[var],
                         N='-'.join(map(str, r['lev'])), err=f2(r['l1']) + ' (self-conv.)', pairs=o2(r['ord']), fit='%.2f' % fine, status=st(fine), arm='p'))
json.dump(rows, open(D + '/conv_rows.json', 'w'), indent=0)
print(len(rows), 'rows;', sum(r['status'] == 'FLAG' for r in rows), 'FLAG')
for r in rows:
    if r['status'] == 'FLAG': print('FLAG', r['geom'], '|', r['case'], '|', r['regime'], r['var'], 'pairs', r['pairs'])
