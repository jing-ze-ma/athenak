#!/usr/bin/env python3
"""raycheap17 (TASK-2026-11-10-orion-raycheap17.md): the 09-25 sp order suite (tests_m1/
runs_5q_sporder2b, st0925.py here) rerun with TODAY's production key set, the bisection
arms, and cubed-sphere twins.

usage (RC17_W = work dir, RC17_BIN = athena binary, RC17_MPI = mpirun command):
  st17.py probe                           S09 1-cycle runs -> $RC17_W/dt17.json (first dt)
  st17.py deck <case> <arm> <mode> <lev>  write one deck, print its path (Part 0 checks)
  st17.py run  <case> <arm> <mode> [lev]  run the levels (parallel, one mpirun each)
  st17.py eval <case> <arm> <mode> [lev]  errors and orders
  st17.py list
cases: pulse_u pulse_s pulse_edd_u fs_u lat_th_atm lat_ph_atm lat_th_mid mvr_u mvp_u (sp
  wedge, exactly the 09-25 case_cfg) and the cubed-sphere twins pulse_cs fs_cs
  pulse_edd_cs
  mvr_cs mvp_cs (radial; N = 8 per face) and latatm_cs latatmv_cs latmid_cs (P_2 mode,
  axis z or the cube vertex (1,1,1); N = 8 s, nx1 = 16 s).
modes: C = space-time (09-25: dx and dt together), X = space only (dt fixed at
  dt_C(finest)/4), T = time only (fixed grid: radial nx1 = 128, lateral s = 2; dt halves).
arms: S09 (09-25 keys, 09-25 dt recipe: C only), S09f (09-25 keys, dt_fixed ladder),
  SP (production keys), SP-<f> (SP minus feature f), S09+<f> (S09 plus feature f), and on
  the cs twins CS09 / CS09f / CSP / CSP-<f> (same key sets + the production cs mesh keys
  for CSP).
"""
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import solib as so      # noqa: E402
import st0925 as s9     # noqa: E402

W = so.W
BIN = os.environ.get('RC17_BIN', '/orion/ptmp/jinma/raycheap14_1109/athena_dea3e3732')
MPI = os.environ.get('RC17_MPI', 'mpirun').split()
NT = 128

SPCASES = ['pulse_u', 'pulse_s', 'pulse_edd_u', 'fs_u', 'lat_th_atm', 'lat_ph_atm',
           'lat_th_mid', 'mvr_u', 'mvp_u']
# cs twin -> (sp case whose rad/prob/T it copies, lateral axis or None)
CSCASES = {'pulse_cs': ('pulse_u', None), 'fs_cs': ('fs_u', None),
           'pulse_edd_cs': ('pulse_edd_u', None), 'mvr_cs': ('mvr_u', None),
           'mvp_cs': ('mvp_u', None), 'latatm_cs': ('lat_th_atm', 'z'),
           'latatmv_cs': ('lat_th_atm', 'vertex'), 'latmid_cs': ('lat_th_mid', 'z')}
# radiation Knudsen number Kn = 1/(rho kappa_t L) of the measured
# structure (L: pulse width
# 0.2, lateral half wavelength r dtheta ~ 0.8 r, atmosphere at its top r = 5): the pass
# rule (vet-acceptance-criteria): p >= 1.8 where Kn < 0.3, p > 0 (converging) where
# Kn >~ 0.3
KN = {'pulse_u': 5.0, 'pulse_s': 5.0, 'pulse_edd_u': 5.0, 'fs_u': 5.0e3,
      'lat_th_atm': 0.06, 'lat_ph_atm': 0.06, 'lat_th_mid': 0.6, 'mvr_u': 5.0e-3,
      'mvp_u': 5.0e-3}

# ---- key sets -----------------------------------------------------------------------


def read_resolved(case):
    """the 09-25 binary's resolved <rad_m1> and <time> of the case (ref0925/, from a
    2-cycle
    run of athena_new_none_cpu = m1-sp-order2b with the 09-25 deck)"""
    out = {'rad_m1': {}, 'time': {}, 'hydro': {}}
    blk = None
    for ln in open(os.path.join(HERE, 'ref0925', '%s.resolved.txt' % case)):
        ln = ln.strip()
        if ln.startswith('<'):
            blk = ln.strip('<>')
            continue
        if '=' in ln and blk in out:
            k, v = [x.strip() for x in ln.split('=', 1)]
            if v != '':
                out[blk][k] = v
    return out


def read_prod(cs):
    d = {}
    for ln in open(os.path.join(HERE, 'prod_keys_rad_m1.txt')):
        if ln.startswith('#') or not ln.strip():
            continue
        k, v, tag = ln.split()
        if tag == 'cs' and not cs:
            continue
        d[k] = v
    return d


# production cubed-sphere mesh keys (HG / RSG decks + run.cfg), named on the CSP arms
CS_MESH_PROD = {'cs_seam_flux': 'shared', 'cs_seam_blockend': 'extend',
                'cs_vertex_3cell': 'true', 'layout_invariant': 'true',
                'cs_seam_rho_guard': '4'}
# SP minus one feature (value None = do not name the key)
MINUS = {
    'col': {'vet_gd': 'false', 'vet_gd_replace': 'false'},
    'nofse': {'implicit_thin_switch': 'none'},
    'lag': {'implicit_recon_lag': 'picard', 'implicit_closure_lag': 'pass'},
    'central': {'implicit_flux': 'central', 'implicit_flux_beam': 'closure',
                'implicit_thin_switch': 'none', 'implicit_thin_corr': 'none'},
    'upw': {'implicit_enthalpy': 'upwind'},
    'nopos': {'implicit_pos_floor': 'false', 'implicit_pos_gas': 'false',
              'implicit_pos_floor_solve': 'false', 'implicit_pos_floor_cons': 'false',
              'implicit_realisable_coupling': 'false', 'time2_posfix': '0',
              'time2_auto_be': '0', 'time2_accept_fallback': 'false',
              'time2_k1_redist': '0', 'time2_pred_cap': '0'},
    'nofrz': {'implicit_thin_freeze': '0'},
    'tolrel': {'implicit_tol_rel': '0'},
    'tight': {'implicit_tol_rel': '0', 'implicit_tol': '1.0e-13',
              'implicit_lin_tol': '1.0e-14', 'implicit_lin_cnorm': '1.0e-14'},
    'pred': {'time2_vet_col': 'predict'},
    'norot': {'vet_gd_rotate_every': '0', 'vet_gd_rotate_dt': '0',
              'vet_rotate_avg': '1', 'vet_rotate_avg_ema': '0'},
    'nc64': {'vet_col_ncore': '64'},
}
# S09 plus one production feature: the production values of these keys (prefix* = all
# production keys with that prefix)
PLUS = {
    'gd': ['vet_gd*', 'vet_rotate*', 'vet_interp', 'vet_seg', 'vet_bc_*', 'vet_scatter*',
           'vet_path', 'vet_col_ncore'],
    'fse': ['implicit_flux', 'implicit_flux_faces', 'implicit_flux_beam',
            'implicit_blend*',
            'implicit_thin_*', 'implicit_hr_*'],
    'lag': ['implicit_recon_lag', 'implicit_closure_lag', 'implicit_recon_npass'],
    'pos': ['implicit_pos_*', 'implicit_realisable_coupling', 'time2_posfix',
            'time2_auto_be', 'time2_accept_fallback', 'time2_k1_*', 'time2_pred_cap'],
    'tol': ['implicit_tol*', 'implicit_lin_tol', 'implicit_lin_cnorm',
            'implicit_conv_est'],
    'x2': ['time2_vet_col', 'time2_vet_col_source', 'time2_start'],
    'frz': ['implicit_thin_freeze'],
}


def plus_keys(feat, prod):
    out = {}
    for pat in PLUS[feat]:
        for k, v in prod.items():
            if (pat.endswith('*') and k.startswith(pat[:-1])) or k == pat:
                out[k] = v
    return out


def arm_cfg(arm):
    """-> (base 's09' | 'sp', overrides, dt_fixed?, cs arm?)"""
    cs = arm.startswith('CS')
    a = arm[1:] if cs else arm         # CS09 -> S09, CSP -> SP, CSP-x -> SP-x
    if a == 'S09':
        return 's09', {}, False, cs
    if a == 'S09f':
        return 's09', {}, True, cs
    if a == 'SP':
        return 'sp', {}, True, cs
    if a.startswith('SP-'):
        return 'sp', dict(MINUS[a[3:]]), True, cs
    if a.startswith('S09+'):
        return 's09', {'__plus__': a[4:]}, True, cs
    raise SystemExit('unknown arm ' + arm)


# ---- decks --------------------------------------------------------------------------


def is_cs(case):
    return case in CSCASES


def sp_twin(case):
    return CSCASES[case][0] if is_cs(case) else case


def lateral(case):
    return s9.case_cfg(sp_twin(case))['geo'] != 'r'


def levels(case, mode):
    if mode == 'T':
        return [1, 2, 4, 8, 16]
    return [1, 2, 4, 8] if lateral(case) else [32, 64, 128, 256, 512]


def dtfix(case, mode, lev):
    """the dt_fixed ladder from the S09 probe (0.9 x the 09-25 first dt: the cap binds)"""
    tab = json.load(open(os.path.join(W, 'dt17.json')))
    c = sp_twin(case)
    if mode == 'C':
        return 0.9*tab['%s|C|%d' % (c, lev)]
    if mode == 'X':
        return 0.9*tab['%s|C|%d' % (c, levels(c, 'C')[-1])]/4.0
    # T: radial = the 09-25 _T level-1 dt (implicit_cfl x8, hydro cfl 0.8) / l; lateral =
    # 4x
    # the C dt at s = 2, / l
    if lateral(case):
        return 0.9*4.0*tab['%s|C|2' % c]/lev
    return 0.9*tab['%s|T|1' % c]/lev


def mk(case, arm, mode, lev, probe=False):
    base, over, fixed, csarm = arm_cfg(arm)
    if csarm != is_cs(case):
        raise SystemExit('arm %s does not match case %s (CS* arms on *_cs cases)'
                         % (arm, case))
    c9 = sp_twin(case)
    # the 09-25 deck of the sp twin at the matching level (C, or the _T variant)
    tcase = c9 + ('_T' if (mode == 'T' and not lateral(case)) else '')
    n9 = lev if mode != 'T' else (lev if not lateral(case) else 2)
    if mode == 'X':
        n9 = lev
    d9, inp9, _, npr = s9.mk(tcase, 'new', n9)
    blocks = parse_deck(inp9)
    cfg = s9.case_cfg(c9)
    rad = blocks['rad_m1']
    tim = blocks['time']
    res = read_resolved(c9)
    # S09: every resolved 09-25 key named at its 09-25 value (today's fresh-run defaults
    # cannot leak in); the deck's own values stay (they were the resolved ones).
    # SP starts
    # from the plain 09-25 case deck: keys production does not name take today's defaults
    if base == 's09':
        for k, v in res['rad_m1'].items():
            rad.setdefault(k, v)
    if base == 'sp':
        prod = read_prod(csarm)
        keep = {}
        if rad.get('closure') != 'vet_col':
            keep['closure'] = rad['closure']   # pulse_edd: Eddington stays
        rad.update(prod)
        rad.update(keep)
        if keep:
            for k in [k for k in rad
                      if k.startswith('vet_gd') or k.startswith('vet_rotate')]:
                del rad[k]
            rad['implicit_thin_switch'] = 'none'
        else:
            rad['vet_gd_rotate_dt'] = '%.12g' % (cfg['T']/16.0)
    if '__plus__' in over:
        rad.update(plus_keys(over.pop('__plus__'), read_prod(csarm)))
        if 'vet_gd_rotate_dt' not in rad and rad.get('vet_gd') == 'true':
            rad['vet_gd_rotate_dt'] = '%.12g' % (cfg['T']/16.0)
    for k, v in over.items():
        if v is None:
            rad.pop(k, None)
        else:
            rad[k] = v
    # time: dt_retry off everywhere (fresh-run default today is on); the fixed ladder
    tim['dt_retry'] = 'false'
    tim['nlim'] = '-1'
    if probe:
        tim['nlim'] = '1'
    elif fixed:
        tim['dt_fixed'] = '%.15g' % dtfix(case, mode, lev)
        tim['cfl_number'] = '0.95'
        rad['implicit_cfl'] = '1.0e6'
    elif mode != 'C':
        raise SystemExit('S09 / CS09 (09-25 dt recipe) is C mode only; use S09f for X/T')
    if is_cs(case):
        blocks, npr = to_cs(case, blocks, mode, lev, csarm and base == 'sp')
    elif lateral(case) and mode == 'T':
        pass   # s = 2 grid from s9.mk(n9 = 2)
    d = os.path.join(W, 'runs', case, arm, '%s_%d' % (mode, lev))
    if probe:
        d = os.path.join(W, 'probe', tcase, str(n9))
    os.makedirs(d, exist_ok=True)
    inp = os.path.join(d, 'in.athinput')
    write_deck(inp, blocks)
    return d, inp, npr


def to_cs(case, blocks, mode, lev, prodmesh):
    """the sp-wedge deck -> a cubed-sphere shell (full sphere, 6 panels)"""
    c9, ax = CSCASES[case]
    m = blocks['mesh']
    if ax is None:
        nf = 8
        n1 = int(m['nx1'])
    else:
        s = lev if mode != 'T' else 2
        nf, n1 = 8*s, 16*s
    for k in [k for k in m
              if k.startswith('use_grid_stretch') or k.startswith('f_stretch')]:
        del m[k]
    m.pop('use_spherical_polar', None)
    m.update({'nx1': n1, 'nx2': nf, 'nx3': nf, 'x2min': -1.0, 'x2max': 1.0, 'x3min': -1.0,
              'x3max': 1.0, 'ix2_bc': 'panel', 'ox2_bc': 'panel', 'ix3_bc': 'panel',
              'ox3_bc': 'panel', 'use_cubed_sphere': 'true'})
    if prodmesh:
        m.update(CS_MESH_PROD)
    blocks['meshblock'] = {'nx1': n1, 'nx2': 8, 'nx3': 8}
    p = blocks['problem']
    for k in ('lat_l2', 'lat_l3'):
        p.pop(k, None)
    if ax is not None:
        p.update({'lat_amp': 0.5, 'lat_l': 2, 'lat_axis': ax})
    p['shell_report'] = 'true'
    p['user_hist'] = 'true'
    # outputs: bin m1 (+ hydro_w) at T, and the user history (mode amplitude) at T
    outs = [k for k in blocks if k.startswith('output')]
    for k in outs:
        blocks[k]['file_type'] = 'bin'
        for q in ('slice_x2', 'slice_x3', 'data_format'):
            blocks[k].pop(q, None)
    T = blocks['time']['tlim']
    blocks['output%d' % (len(outs) + 1)] = {'file_type': 'hst', 'dt': T,
                                            'data_format': '%.17e'}
    nblk = 6*(nf//8)**2
    return blocks, min(nblk, 48)


def parse_deck(path):
    blocks, b = {}, None
    for ln in open(path):
        ln = ln.split('#')[0].rstrip()
        if ln.startswith('<'):
            b = ln.strip('<> ')
            blocks[b] = {}
        elif '=' in ln and b:
            k, v = [x.strip() for x in ln.split('=', 1)]
            blocks[b][k] = v
    return blocks


def write_deck(path, blocks):
    with open(path, 'w') as f:
        for b, d in blocks.items():
            f.write('<%s>\n' % b)
            for k, v in d.items():
                f.write('%-34s = %s\n' % (k, v))
            f.write('\n')


# ---- runs ---------------------------------------------------------------------------


def launch(jobs):
    procs = []
    for d, inp, npr in jobs:
        for f in os.listdir(d):
            if f in ('tab', 'bin'):
                subprocess.call(['rm', '-rf', os.path.join(d, f)])
        cmd = MPI + ['-np', str(npr), BIN, '-i', inp, '-d', d]
        lg = open(os.path.join(d, 'log'), 'w')
        env = dict(os.environ, OMP_NUM_THREADS='1')
        procs.append((d, subprocess.Popen(cmd, stdout=lg, stderr=subprocess.STDOUT,
                                          env=env)))
    for d, p in procs:
        print(d, 'rc', p.wait(), flush=True)


def probe():
    jobs = []
    for c in SPCASES:
        for n in levels(c, 'C'):
            jobs.append(mk(c, 'S09', 'C', n, probe=True))
        if not lateral(c):
            jobs.append(mk(c, 'S09', 'T', 1, probe=True))
    launch(jobs)
    tab = {}
    for c in SPCASES:
        for mode, ns in (('C', levels(c, 'C')), ('T', [1] if not lateral(c) else [])):
            for n in ns:
                tc = c + ('_T' if mode == 'T' else '')
                lg = os.path.join(W, 'probe', tc, str(n), 'log')
                for ln in open(lg):
                    if 'cycle=0 ' in ln and ' dt=' in ln:
                        tab['%s|%s|%d' % (c, mode, n)] = float(
                            ln.split('dt=')[1].split()[0])
                        break
    json.dump(tab, open(os.path.join(W, 'dt17.json'), 'w'), indent=1)
    print(json.dumps(tab, indent=1))


# ---- evaluation ---------------------------------------------------------------------


def shell(d):
    """per radial index Emean, F1mean (RadM1ShellReport, double) and the lateral spread"""
    r, e, f, sp = [], [], [], 0.0
    for ln in open(os.path.join(d, 'log')):
        t = ln.split()
        if len(t) == 10 and t[0] == 'SHELL' and t[1].isdigit():
            v = [float(x) for x in t[2:]]
            r.append(v[0])
            e.append(v[3])
            f.append(v[6])
            sp = max(sp, (v[4] - v[2])/max(abs(v[3]), 1e-300))
    return np.array(e), np.array(f), sp


def hst_amp(d):
    fs = sorted(x for x in os.listdir(d) if x.endswith('.user.hst') or x.endswith('.hst'))
    for fn in fs:
        rows = [ln.split() for ln in open(os.path.join(d, fn)) if not ln.startswith('#')]
        hdr = [ln for ln in open(os.path.join(d, fn))
               if ln.startswith('#') and 'VEP' in ln]
        if not rows or not hdr:
            continue
        names = [h.split('=')[-1] for h in hdr[-1].replace('#', '').split()]
        v = dict(zip(names, [float(x) for x in rows[-1]]))
        return ((v['VEP'] - v['VE']*v['VP']/v['V'])/(v['VPP'] - v['VP']**2/v['V']),
                v['time'] if 'time' in v else float(rows[-1][0]))
    return None, None


def binflat(d):
    """every bin variable flattened in file (block) order: same grid at every T level"""
    import bin_convert
    out = {}
    for fn in sorted(os.listdir(os.path.join(d, 'bin'))):
        if not fn.endswith('.bin'):
            continue
        b = bin_convert.read_binary(os.path.join(d, 'bin', fn))
        for v in b['var_names']:
            out[v] = np.concatenate([np.ravel(x) for x in b['mb_data'][v]])
    return out


def load_sp(case, arm, mode, n):
    """the sp wedge variables as st0925.load does, from this arm's directory"""
    d = os.path.join(W, 'runs', case, arm, '%s_%d' % (mode, n))
    c = s9.case_cfg(case)
    out = {}
    sub = 'tab' if c['geo'] == 'r' else 'bin'
    for v in (['m1'] + (['hydro_w'] if c['gas'] else [])):
        fs = sorted(f for f in os.listdir(os.path.join(d, sub))
                    if f.endswith('.' + sub) and ('.%s.' % v) in f)
        t = (so.tab if sub == 'tab' else so.binall)(os.path.join(d, sub, fs[-1]))
        out.update(t)
    q = {'E': out['m1_e'], 'F1': out['m1_f1']}
    if c['geo'] != 'r':
        a = 2 if c['geo'] == 'th' else 3
        q['F%d' % a] = out['m1_f%d' % a]
    if c['gas']:
        rho = out['dens']
        q['T'] = s9.GM1*out['eint']/rho
        q['M1'] = rho*out['velx']
        q['Eg'] = out['eint'] + 0.5*rho*(out['velx']**2 + out['vely']**2 + out['velz']**2)
    q['time'] = out['time']
    return q


def evaluate(case, arm, mode, ns):
    res = {}
    if is_cs(case):
        dd = {n: os.path.join(W, 'runs', case, arm, '%s_%d' % (mode, n)) for n in ns}
        if lateral(case):
            A = {n: hst_amp(dd[n])[0] for n in ns}
            res['Amode'] = [abs(A[a] - A[b])/abs(A[b]) for a, b in zip(ns[:-1], ns[1:])]
        if mode == 'T':
            Q = {n: binflat(dd[n]) for n in ns}
            for k in Q[ns[0]]:
                if k in ('m1_e', 'm1_f1', 'm1_f2', 'm1_f3', 'dens', 'eint', 'velx'):
                    res[k] = [np.abs(Q[a][k] - Q[b][k]).mean()
                              / np.abs(Q[b][k] - Q[b][k].mean()).mean()
                              for a, b in zip(ns[:-1], ns[1:])]
        else:
            S = {n: shell(dd[n]) for n in ns}
            r1 = s9.case_cfg(sp_twin(case))['r1']
            for j, k in ((0, 'Emean'), (1, 'F1mean')):
                errs = []
                for a, b in zip(ns[:-1], ns[1:]):
                    ec, eb = S[a][j], S[b][j]
                    ef = so.restrict_r(eb, so.rfaces(1.0, r1, len(eb))) \
                        if len(eb) == 2*len(ec) else eb
                    w = np.diff(so.rfaces(1.0, r1, len(ec))**3)
                    errs.append(np.sum(w*np.abs(ec - ef))
                                / np.sum(w*np.abs(ec - np.sum(w*ec)/np.sum(w))))
                res[k] = errs
            res['spread'] = [S[n][2] for n in ns]
        return res
    c = s9.case_cfg(case)
    if mode == 'T':
        c = dict(c, tref=True)
    Q = {n: load_sp(case, arm, mode, n) for n in ns}
    for k in [k for k in Q[ns[0]] if k != 'time']:
        errs = []
        for a, b in zip(ns[:-1], ns[1:]):
            ec = Q[a][k]
            ef = Q[b][k] if mode == 'T' else s9.restrict(Q[b][k], c, a)
            if mode == 'T' and lateral(case):
                w = s9.vol(dict(c, tref=False), 2, ec.shape)
            else:
                w = s9.vol(c, a, ec.shape)
            mean = np.sum(w*ec)/np.sum(w)
            errs.append(np.sum(w*np.abs(ec - ef))/np.sum(w*np.abs(ec - mean)))
        res[k] = errs
    return res


def verdict(case, p):
    kn = KN[sp_twin(case)]
    if kn < 0.3:
        return 'PASS' if p >= 1.8 else 'FAIL(thick p<1.8)'
    tag = '' if p >= 1.8 else ' [<1.8: loss vs 09-25 -> bisect]'
    return ('PASS' if p > 0.0 else 'FAIL(thin no conv)') + tag


def main():
    act = sys.argv[1]
    if act == 'list':
        print('sp cases', ' '.join(SPCASES))
        print('cs cases', ' '.join(CSCASES))
        print('arms S09 S09f SP', ' '.join('SP-' + k for k in MINUS),
              ' '.join('S09+' + k for k in PLUS), '(CS-prefixed on cs cases)')
        return
    if act == 'probe':
        probe()
        return
    case, arm, mode = sys.argv[2:5]
    ns = [int(x) for x in sys.argv[5:]] or levels(case, mode)
    if act == 'deck':
        print(mk(case, arm, mode, ns[0])[1])
        return
    if act == 'run':
        launch([mk(case, arm, mode, n) for n in ns])
        return
    res = evaluate(case, arm, mode, ns)
    for k, e in res.items():
        if k == 'spread':
            print('%-12s %-12s %s %-7s lateral spread max(Emax-Emin)/Emean per level %s'
                  % (case, arm, mode, k, ' '.join('%.2e' % x for x in e)))
            continue
        p = so.order(e)
        print('%-12s %-12s %s %-7s d %s | p %s | %s' % (
            case, arm, mode, k, ' '.join('%.3e' % x for x in e),
            ' '.join('%5.2f' % x for x in p), verdict(case, p[-1]) if len(p) else ''))


if __name__ == '__main__':
    main()
