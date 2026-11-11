"""sptests-1011: helpers for the production-path (vet_col + vet_gd, sp wedge) test suite.

mkinput(): an athinput from ordered dicts.  tab()/binall(): readers.  logstats(): Picard and
cost from a run log.  rfaces(): radial faces (uniform / poly as CoordSphericalPolar).
ARMS: the arm keys (NOTE-sptests.md 0.3).
"""
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, '/viper/ptmp2/jinma/blend_1022/wt/vis/python')
W = '/viper/ptmp2/jinma/blend_1022/p2/ana/sp'
# SPT_RUNS: the run tree (runs = the 82ef2e9c binary; runs_pc = prodcand)
RUNS = os.path.join(W, os.environ.get('SPT_RUNS', 'runs'))
BIN = {'new': os.environ.get('SPT_BIN', W + '/bin/athena_cpu_new')}
# the runs_5q poly stretch (refines toward the OUTER end by ~2.4x, as production refines
# toward the photosphere)
POLY = {'f_stretch_r_c1': '+0.134559', 'f_stretch_r_c2': '+4.942505',
        'f_stretch_r_c3': '-8.950106', 'f_stretch_r_c4': '+4.403042'}
TH0, TH1, PH1 = 1.1707963267948966, 1.9707963267948966, 1.0471975511965976

# ---- arms --------------------------------------------------------------------------------
GD = {'closure': 'vet_col', 'vet_gd': 'true', 'vet_gd_replace': 'true',
      'vet_col_lat_taucut': '30.0', 'vet_gd_gl_nmu': 10, 'vet_gd_gl_nphi': 12,
      'vet_gd_gl_stagger': 'true', 'vet_gd_rotate_every': 10,
      'vet_col_lat_offdiag': 'operator', 'implicit_op_stencil': 'true',
      'vet_gd_tangential': 'false', 'vet_gd_wall_interp': 'true', 'vet_gd_twin': 'true',
      'vet_gd_twin_fuse': 'true', 'vet_gd_twin_det': 'true', 'vet_source_noesrc': 'true',
      'vet_scatter': 'true', 'vet_scatter_form': 'ma', 'vet_scatter_j': 'relaxed',
      'implicit_face_opac_n': 'true', 'implicit_pos_floor': 'true',
      'implicit_pos_floor_solve': 'true', 'vet_gd_allow_clamp': 'true'}
HR = {'implicit_flux': 'blend', 'implicit_flux_faces': 'all',
      'implicit_flux_beam': 'halfrange', 'implicit_blend': 'idort_f',
      'implicit_blend_tau0': 1, 'implicit_blend_alpha': 1, 'implicit_blend_flo': 0.3,
      'implicit_blend_fhi': 0.6, 'implicit_blend_xthin': 30,
      'implicit_blend_xthin_mode': 'all', 'implicit_blend_xthin_wmin': 0.0,
      'implicit_blend_xthin_r0': 1.5, 'implicit_hr_recon': 'dc',
      'implicit_hr_recon_fresh': 2, 'implicit_hr_damp': 0.0, 'implicit_hr_pos': 'kill'}
FH = {'implicit_blend_xthin': 0, 'implicit_thin_switch': 'fse', 'implicit_thin_corr': 'hc',
      'implicit_thin_recon': 'dc', 'implicit_thin_bq': 'cell'}
H2 = {'time_scheme': 'hesdirk2', 'time2_vet_col': 'predict', 'implicit_vimp': 'true',
      'time2_vstage': 'true'}
CART = {'closure': 'vet_sc', 'vet_tensor': 'full', 'vet_nmu': 8, 'vet_nphi': 32}


def arm_keys(arm):
    """arm string: [H]{P|F|CP|CF}{c|f}[R0][A][D][N4]  (see NOTE 0.3)
    H hesdirk2+predict+vimp+vstage; P production hr (xthin 30); F fse+hc; CP / CF the
    Cartesian vet_sc analogues; c / f vet_bc_diff cell / face; R0 rotate_every 0; A taucut
    1e30 (vet_gd everywhere); D vet_gd_inner = diff (the old start); N4 gl_nmu 4"""
    k = {}
    # final_1012: 'ARM+key=val+key=val' -> explicit rad_m1 overrides applied last
    ovr = {}
    if '+' in arm:
        parts = arm.split('+')
        arm = parts[0]
        ovr = dict(x.split('~', 1) for x in parts[1:])
    a = arm
    h = a.startswith('H')
    if h:
        a = a[1:]
    cart = a.startswith('C')
    if cart:
        a = a[1:]
    base = a[0]
    a = a[1:]
    if cart:
        k.update(CART)
    else:
        k.update(GD)
    k.update(HR)
    if base == 'F':
        k.update(FH)
    elif base != 'P':
        raise SystemExit('bad arm ' + arm)
    k['time_scheme'] = 'be'
    if h:
        k.update(H2)
        if cart:
            k.pop('time2_vet_col')
            k['time2_vet_sc'] = 'predict'
    if a.startswith('c') or a.startswith('f'):
        k['vet_bc_diff'] = 'cell' if a[0] == 'c' else 'face'
        a = a[1:]
    while a:
        if a.startswith('R0'):
            k['vet_gd_rotate_every'] = 0
            a = a[2:]
        elif a.startswith('A'):
            k['vet_col_lat_taucut'] = '1.0e30'
            a = a[1:]
        elif a.startswith('D'):
            k['vet_gd_inner'] = 'diff'
            a = a[1:]
        elif a.startswith('G') and a[1:2].isdigit():
            # G<k>: k x the vet_gd directions per axis (gl_nmu 10k x gl_nphi 12k)
            g = int(a[1])
            k['vet_gd_gl_nmu'] = 10*g
            k['vet_gd_gl_nphi'] = 12*g
            a = a[2:]
        elif a.startswith('T0'):
            # vet_gd_twin off (no ray-noise subtraction)
            k['vet_gd_twin'] = 'false'
            k.pop('vet_gd_twin_fuse', None)
            k.pop('vet_gd_twin_det', None)
            a = a[2:]
        elif a.startswith('Q'):
            # 18 x 24 = 432 vet_gd directions (the largest set M1_VGD_NMAX = 432 allows)
            k['vet_gd_gl_nmu'] = 18
            k['vet_gd_gl_nphi'] = 24
            a = a[1:]
        elif a.startswith('X'):
            # sptests P2 limiter: implicit_hr_excess = iso
            k['implicit_hr_excess'] = 'iso'
            a = a[1:]
        elif a.startswith('Jc'):
            # vet_scatter_j = current: the scattering part of the relaxed source takes E,
            # not the relaxed E*
            k['vet_scatter_j'] = 'current'
            a = a[2:]
        elif a.startswith('St'):
            # hesdirk2: time2_vet_col_source = state (the source of the predicted state)
            k['time2_vet_col_source'] = 'state'
            a = a[2:]
        elif a.startswith('Sg'):
            # vet_col_source = gas (the default relaxed source -> B at huge c dt kappa_a)
            k['vet_col_source'] = 'gas'
            a = a[2:]
        elif a.startswith('M'):
            # sptestsC: the gas source named explicitly (matched = the default)
            k['implicit_gas_src'] = 'matched'
            a = a[1:]
        elif a.startswith('W'):
            # suiteF-1012: the FINAL production gas source and keys (all named in the deck):
            # conservative2 without the hydro shift, vimp at the Marshak end face, vet_seg emit,
            # the Lie start (time2_start = be, the default)
            k['implicit_gas_src'] = 'conservative2'
            k['implicit_cs2_shift'] = 'false'
            if os.environ.get('SPT_HLL'):
                # 4253af30: the conservative2 HLL / physical-face rule, named
                k['implicit_cs2_hll'] = os.environ['SPT_HLL']
            k['implicit_vimp_bc'] = 'true'
            k['vet_seg'] = 'emit'
            if k.get('time_scheme') == 'hesdirk2':
                k['time2_start'] = 'be'
            a = a[1:]
        elif a.startswith('V'):
            # suiteF-1012 PRODm: the production keys with implicit_gas_src = matched
            k['implicit_gas_src'] = 'matched'
            k['implicit_vimp_bc'] = 'true'
            k['vet_seg'] = 'emit'
            if k.get('time_scheme') == 'hesdirk2':
                k['time2_start'] = 'be'
            a = a[1:]
        elif a.startswith('Y'):
            # suiteF-1012: time2_start = slope (opt-in start, T ladders only)
            k['time2_start'] = 'slope'
            a = a[1:]
        elif a.startswith('K'):
            # sptestsC: the conservative gas source (consrc-1011, opt-in)
            k['implicit_gas_src'] = 'conservative'
            a = a[1:]
        elif a.startswith('Z'):
            # final_1012 (worker F): FIN = the combined candidate method (method-1012), every key
            # named.  vet_gd_rotate_dt is set per run in spt.mk (SPT_ROT: 'auto' = 20 x the run's
            # dt, T ladders 10 x the coarsest dt; '0' = rotate 0 control)
            k.update({'implicit_blend_xthin': 0, 'implicit_thin_switch': 'fse',
                      'implicit_thin_corr': 'hc', 'implicit_thin_hcsf': 'j',
                      'implicit_thin_recon': 'jse', 'implicit_thin_js_gd': 'x1',
                      'implicit_thin_bq': 'lin', 'implicit_thin_bq_gd': 'cell',
                      'implicit_thin_hcface': 'tau', 'implicit_thin_hcface_gd': 'x1',
                      'vet_seg': 'jump', 'vet_bc_diff': 'face', 'vet_gd_symref': 'col',
                      'time2_vet_col_source': 'blend', 'implicit_gas_src': 'conservative2',
                      'implicit_cs2_hll': 'drag', 'implicit_cs2_shift': 'false',
                      'implicit_vimp_bc': 'true', 'time_scheme': 'hesdirk2',
                      'time2_vet_sc': 'predict', 'time2_vet_col': 'predict',
                      'time2_vstage': 'true', 'implicit_vimp': 'true',
                      'time2_pred_cap': 0.5, 'time2_start': 'be',
                      'vet_gd_rotate_align': 'true', 'vet_gd_rotate_t0': 0.0,
                      'vet_gd_rotate_stretch': 0.0, 'vet_gd_rotate_slope': 'keep',
                      'vet_gd_rotate_dt': 'AUTO'})
            a = a[1:]
        elif a.startswith('F0'):
            # final_1012: implicit_f0_init = zero (radiation co-moving IC; moving-gas ICs, NOTE-sporder)
            k['implicit_f0_init'] = 'zero'
            a = a[2:]
        elif a.startswith('O'):
            # final_1012: OLD = the new keys at their old values (bitwise sanity vs 4253af30 HFfW drag)
            k.update({'implicit_thin_hcsf': 'j', 'implicit_thin_js_gd': 'none',
                      'implicit_thin_recon': 'dc', 'implicit_thin_bq': 'cell'})
            a = a[1:]
        elif a.startswith('N4'):
            k['vet_gd_gl_nmu'] = 4
            a = a[2:]
        else:
            raise SystemExit('bad arm suffix ' + arm)
    if cart:
        for x in [x for x in k if x.startswith('vet_gd') or x.startswith('vet_col')]:
            k.pop(x)
    k.update(ovr)
    return k


# ---- io ------------------------------------------------------------------------------------
def mkinput(path, blocks):
    """blocks: list of (name, dict)"""
    with open(path, 'w') as f:
        for b, d in blocks:
            f.write('<%s>\n' % b)
            for k, v in d.items():
                f.write('%-28s = %s\n' % (k, v))
            f.write('\n')


def run(binary, inp, rdir, np_=1, nthr=1, args=(), wait=True, rst=None):
    os.makedirs(rdir, exist_ok=True)
    src = ['-r', rst] if rst else ['-i', inp]
    cmd = ['mpirun', '-np', str(np_), '--oversubscribe', '--bind-to', 'none', binary
           ] + src + ['-d', rdir] + list(args)
    env = dict(os.environ, OMP_NUM_THREADS=str(nthr), OMP_PROC_BIND='false')
    lg = open(os.path.join(rdir, 'log'), 'w')
    p = subprocess.Popen(cmd, stdout=lg, stderr=subprocess.STDOUT, env=env)
    if wait:
        return p.wait()
    return p


def tab(fn):
    """a tab slice: dict of 1-D arrays + time"""
    with open(fn) as f:
        l0 = f.readline()
        l1 = f.readline()
    t = float(re.search(r'time=\s*([-+0-9.eE]+)', l0).group(1))
    cyc = int(re.search(r'cycle=\s*(\d+)', l0).group(1))
    names = l1[1:].split()
    a = np.loadtxt(fn, comments='#', ndmin=2)
    if len(names) > a.shape[1]:
        # a 2-D Cartesian slice: the header names k / x3v columns the data lacks
        names = names[:3] + names[len(names) - (a.shape[1] - 3):]
    out = {n: a[:, q] for q, n in enumerate(names)}
    out['time'] = t
    out['cycle'] = cyc
    return out


def tabs(d, var='m1'):
    fs = sorted(f for f in os.listdir(os.path.join(d, 'tab'))
                if f.endswith('.tab') and ('.%s.' % var) in f)
    return [os.path.join(d, 'tab', f) for f in fs]


def binall(fn):
    import bin_convert
    d = bin_convert.read_binary(fn)
    n1, n2, n3 = d['Nx1'], d['Nx2'], d['Nx3']
    b1, b2, b3 = d['nx1_out_mb'], d['nx2_out_mb'], d['nx3_out_mb']
    out = {}
    for v in d['var_names']:
        a = np.zeros((n3, n2, n1))
        for m, blk in enumerate(d['mb_data'][v]):
            lx1, lx2, lx3 = d['mb_logical'][m][:3]
            a[lx3*b3:(lx3+1)*b3, lx2*b2:(lx2+1)*b2, lx1*b1:(lx1+1)*b1] = blk
        out[v] = a
    out['time'] = d['time']
    return out


def logstats(d):
    """Picard mean/max, NON-CONVERGED, cycles, cpu time, clamps, FATAL/nan"""
    s = {}
    fn = os.path.join(d, 'log')
    if not os.path.exists(fn):
        return s
    txt = open(fn, errors='replace').read()
    m = re.search(r'Picard iterations mean=([-+0-9.eE]+) max=([-+0-9.eE]+) '
                  r'NON-CONVERGED=([-+0-9.eE]+)', txt)
    if m:
        s['pic'] = float(m.group(1))
        s['picmax'] = float(m.group(2))
        s['nonconv'] = float(m.group(3))
    m = re.search(r'cpu time used\s*=\s*([-+0-9.eE]+)', txt)
    if m:
        s['cpu'] = float(m.group(1))
    m = re.findall(r'cycle=(\d+)', txt)
    if m:
        s['cyc'] = int(m[-1])
    m = re.search(r'inner iterations mean=([-+0-9.eE]+)', txt)
    if m:
        s['inner'] = float(m.group(1))
    s['fatal'] = ('FATAL' in txt) or ('nan' in txt.lower().split('nanosec')[0] and
                                      ' nan' in txt.lower())
    s['clamp'] = len(re.findall(r'CLAMPED', txt))
    s['done'] = ('Terminating' in txt)
    return s


# ---- grids ---------------------------------------------------------------------------------
def rfaces(r0, r1, n, poly=False):
    xi = np.arange(n + 1)/n
    u = xi.copy()
    if poly:
        xik = xi.copy()
        for m in range(4):
            u += float(POLY['f_stretch_r_c%d' % (m + 1)])*xik*(1.0 - xi)
            xik = xik*xi
    return r0 + (r1 - r0)*u


def rcent(rf, sp=True):
    """the sp volume centroid (x1v): 3/4 (r_r^4 - r_l^4)/(r_r^3 - r_l^3)"""
    if not sp:
        return 0.5*(rf[1:] + rf[:-1])
    a, b = rf[:-1], rf[1:]
    return 0.75*(b**4 - a**4)/(b**3 - a**3)


def restrict_r(a, rf, sp=True):
    w = np.diff(rf**3) if sp else np.diff(rf)
    s = a*w
    return (s[..., 0::2] + s[..., 1::2])/(w[0::2] + w[1::2])


def order(e):
    e = np.asarray(e, dtype=float)
    return np.log2(e[:-1]/e[1:])


def fmt(x):
    return ' '.join('%.3e' % v for v in x)


def fmto(x):
    return ' '.join('%5.2f' % v for v in x)
