#!/usr/bin/env python3
"""order_1009: space / time / combined convergence of the implicit VET face flux
(half-range blend vs central) on the Cartesian vet_sc slab (nx2 = 4, one block in x1).

usage: od.py run  <jobfile> [P]       run every line of the job file (xargs-like pool)
       od.py jobs <study> > jobfile    print the job lines of a study
       od.py eval <study>               print the convergence table of a study
A run = (case, arm, scheme, mode, level).  mode: X space (dt fixed), T time (grid fixed),
C combined (implicit_cfl fixed).  Self-convergence: e = L1(q_l - R q_l+1)/L1(|q_l - <q_l>|)
(R = 2:1 restriction for X/C; identity for T), per variable E, F1.
"""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

W = os.environ['XTF_RUN']                      # run tree
B = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # bundle dir
TPL = {'pulse': B + '/inp/pulse.athinput', 'atm': B + '/inp/atm.athinput'}

HR = {'implicit_flux': 'blend', 'implicit_flux_faces': 'all',
      'implicit_flux_beam': 'halfrange', 'implicit_blend': 'idort_f',
      'implicit_blend_tau0': '1.0', 'implicit_blend_flo': '0.3',
      'implicit_blend_fhi': '0.6', 'implicit_blend_xthin': '30', 'implicit_blend_alpha': '1.0'}
ARMS = {'cen': {'implicit_flux': 'central', 'implicit_flux_faces': 'x1'},
        'hr': HR,
        'hr1': dict(HR, implicit_flux='berthon'),     # pure half-range (w = 1)
        'hrx0': dict(HR, implicit_blend_xthin='0.0'),  # no transparency override
        }
BASE = {'transport': 'implicit', 'closure': 'vet_sc', 'implicit_solver': 'bicgstab',
        'implicit_tol': '1.0e-12', 'implicit_lin_tol': '1.0e-13', 'implicit_maxit': '200',
        'implicit_lin_maxit': '400', 'implicit_vimp': 'false'}
ARMS['hrb'] = dict(HR, implicit_blend_xthin_mode='beam')  # xthinfix-1009
ARMS['hrs'] = dict(HR, implicit_blend_xthin_mode='steep')      # steep, X0 30
ARMS['hrs15'] = dict(ARMS['hrs'], implicit_blend_xthin='15')
ARMS['hrs60'] = dict(ARMS['hrs'], implicit_blend_xthin='60')
ARMS['hrsp'] = dict(ARMS['hrs'], implicit_hr_recon='plm')      # steep + plm
ARMS['hrp'] = dict(HR, implicit_hr_recon='plm')                # all + plm
ARMS['cenv'] = dict(ARMS['cen'], implicit_vimp='true')
TPL.update({'rw_t10': B + '/inp/rw_t10.athinput', 'rw_t1000': B + '/inp/rw_t1000.athinput'})


def case_cfg(case):
    """case: pulse_k<kappa>  (periodic Gaussian E pulse, pure scattering)"""
    if case.startswith('pulse_k'):
        k = float(case[7:])
        # thick: T = diffusion time of the width (sigma^2 doubles), thin: 0.2 light units
        T = max(0.2, 0.00375*k)
        # periodic formal solution: npass = ceil(14*0.95/(kappa L)) + 1, capped at 16
        npass = min(int(np.ceil(14*0.95/k)) + 1, 16)
        return {'tpl': 'pulse', 'T': T, 'L': 1.0,
                'rad': {'kappa_s': '%g' % k, 'vet_x1_periodic': 'true',
                        'vet_x1_npass': str(npass)}}
    if case == 'atm':
        # grey plane-parallel atmosphere (G1v), tau_total ~ 100, steady state vs Hopf
        return {'tpl': 'atm', 'T': 2000.0, 'rad': {}, 'L': 1.0, 'ic': 1.0e4}
    if case.startswith('rw_t'):
        # linear radiation-acoustic wave, P_rad/P_gas = 100, tau_lambda = 10 | 1e3, c = 1e3,
        # exact eigenmode IC (runs_3r_radwave), one period; hydro cfl 0.8 must not limit
        import re
        tl = float(re.search(r'tlim = ([0-9.e+-]+)', open(TPL[case]).read()).group(1))
        # hydro-limited dt (cfl 0.8): 251 (tau 10) / 25 (tau 1e3) steps at N = 32
        st = {'rw_t10': (4200, 1100, 270), 'rw_t1000': (500, 128, 40)}[case]
        return {'tpl': case, 'T': tl, 'rad': {}, 'L': 1.0, 'c': 1.0e3, 'nX': st[0],
                'nT': st[1], 'nC': st[2], 'NT': 128, 'hcfl': '0.8', 'gas': True}
    raise SystemExit('unknown case ' + case)


# levels: X/C: nx1 = 32..512; T: nx1 = NT, dt = dt0/l, l = 1..16
NX = [32, 64, 128, 256, 512]
NT = 256
TL = [1, 2, 4, 8, 16]


def rundir(case, arm, sch, mode, lev, tag=''):
    return os.path.join(W, 'runs', case, '%s_%s_%s%s_%d' % (arm, sch, mode, tag, lev))


def write_input(path, tpl, mesh, rad, time):
    s = open(TPL[tpl]).read().split('\n')
    out, blk = [], None
    sets = {'mesh': dict(mesh), 'meshblock': {'nx1': mesh['nx1']}, 'time': dict(time),
            'rad_m1': dict(rad)}
    for ln in s:
        st = ln.strip()
        if st.startswith('<') and st.endswith('>'):
            if blk in sets:
                for k, v in sets[blk].items():
                    out.append('%-24s = %s' % (k, v))
                sets[blk] = {}
            blk = st[1:-1]
            out.append(ln)
            continue
        if blk in sets and '=' in st and not st.startswith('#'):
            k = st.split('=')[0].strip()
            if k in sets[blk]:
                out.append('%-24s = %s' % (k, sets[blk].pop(k)))
                continue
        out.append(ln)
    if blk in sets:
        for k, v in sets[blk].items():
            out.append('%-24s = %s' % (k, v))
    open(path, 'w').write('\n'.join(out) + '\n')


def mk(case, arm, sch, mode, lev, cfl, tag=''):
    """cfl: X: the implicit_cfl at nx1 = 512 / 4 (dt fixed for all levels);
    T: the implicit_cfl at l = 1 on nx1 = NT; C: the implicit_cfl of every level"""
    c = case_cfg(case)
    d = rundir(case, arm, sch, mode, lev, tag)
    os.makedirs(d, exist_ok=True)
    T = c['T']
    cl = c.get('c', 1.0)
    nX, nT, nC, NT_ = c.get('nX', 300), c.get('nT', 40), c.get('nC', 40), c.get('NT', NT)
    if mode == 'X':      # dt = T/nX on every grid
        n = lev
        ic = cl*T/nX*n
    elif mode == 'T':    # nx1 = NT, dt = T/(nT l)
        n = NT_
        ic = cl*T/(nT*lev)*NT_
    elif mode == 'S':    # steady state at implicit_cfl = case ic
        n = lev
        ic = c['ic']
    else:                # c dt/dx fixed: dt = T/nC at nx1 = 32
        n = lev
        ic = cl*T/nC*32
    rad = dict(BASE, time_scheme=sch, implicit_cfl='%.12g' % ic)
    rad.update(c['rad'])
    rad.update(ARMS[arm])
    time = {'tlim': '%.17g' % c['T'], 'cfl_number': c.get('hcfl', '1.0e6'), 'nlim': '-1'}
    if c['tpl'] == 'atm':
        time['tlim'] = '%.12g' % c['T']
    mesh = {'nx1': n}
    inp = os.path.join(d, 'in.athinput')
    write_input(inp, c['tpl'], mesh, rad, time)
    return d, inp


def jobs(study):
    """study = case:arm1,arm2:sch1,sch2:modes:cfl  (modes e.g. XTC)"""
    case, arms, schs, modes, cfl = study.split(':')
    cfl = float(cfl)
    lines = []
    for arm in arms.split(','):
        for sch in schs.split(','):
            for mode in modes:
                for lev in (TL if mode == 'T' else NX):
                    d, inp = mk(case, arm, sch, mode, lev, cfl)
                    lines.append('%s %s' % (d, inp))
    return lines


def run_one(line):
    d, inp = line.split()
    if os.path.exists(os.path.join(d, 'wall.txt')):
        return d + ' skip'
    subprocess.call(['bash', B + '/one.sh', d, inp])
    return d + ' ' + open(os.path.join(d, 'wall.txt')).read().strip()


def load(d):
    tdir = os.path.join(d, 'tab')
    fs = sorted(f for f in os.listdir(tdir) if f.endswith('.tab'))
    a = np.loadtxt(os.path.join(tdir, fs[-1]))
    hdr = open(os.path.join(tdir, fs[-1])).readline()
    # columns: i x1v E F1 F2 F3 (athenak tab: gid i x1v vars...)
    return a, hdr


def cols(d, var='m1'):
    tdir = os.path.join(d, 'tab')
    fs = sorted(f for f in os.listdir(tdir) if f.endswith('.tab') and '.%s.' % var in f)
    fn = os.path.join(tdir, fs[-1])
    names = None
    with open(fn) as f:
        for ln in f:
            if ln.startswith('#') and 'x1v' in ln:
                names = ln[1:].split()
            if not ln.startswith('#'):
                break
    a = np.loadtxt(fn)
    t = float(open(fn).readline().split('time=')[1].split()[0])
    names = [x for x in names if x not in ('k', 'x3v')]
    out = {k: a[:, names.index(k)] for k in names}
    out['time'] = t
    return out


def evaluate(case, arm, sch, mode, tag='', levs=None):
    levs = levs or (TL if mode == 'T' else NX)
    Q = {}
    for lev in levs:
        d = rundir(case, arm, sch, mode, lev, tag)
        try:
            Q[lev] = cols(d)
            if Q[lev]['time'] < 0.999999*case_cfg(case)['T']:
                del Q[lev]
                break
            if case.startswith('rw'):
                Q[lev].update(cols(d, 'hydro_w'))
        except Exception:
            break
    levs = [lv for lv in levs if lv in Q]
    res = {}
    for k in ('m1_e', 'm1_f1', 'dens', 'velx', 'eint'):
        if k not in Q[levs[0]]:
            continue
        errs = []
        for a, b in zip(levs[:-1], levs[1:]):
            qa, qb = Q[a][k], Q[b][k]
            if mode != 'T':
                qb = 0.5*(qb[0::2] + qb[1::2])
            mean = qa.mean()
            errs.append(np.abs(qa - qb).sum()/np.abs(qa - mean).sum())
        res[k] = errs
    return res, levs


def main():
    act = sys.argv[1]
    if act == 'jobs':
        print('\n'.join(jobs(sys.argv[2])))
    elif act == 'run':
        P = int(sys.argv[3]) if len(sys.argv) > 3 else 12
        lines = [ln for ln in open(sys.argv[2]).read().split('\n') if ln.strip()]
        # longest (finest) first
        with ThreadPoolExecutor(P) as ex:
            for r in ex.map(run_one, lines[::-1]):
                print(r, flush=True)
    elif act == 'eval':
        case, arms, schs, modes, cfl = sys.argv[2].split(':')
        for arm in arms.split(','):
            for sch in schs.split(','):
                for mode in modes:
                    lv = None
                    if len(sys.argv) > 3:
                        for tok in sys.argv[3].split(';'):
                            md, ls = tok.split('=')
                            if md == mode:
                                lv = [int(x) for x in ls.split(',')]
                    res, levs = evaluate(case, arm, sch, mode, levs=lv)
                    for k, e in res.items():
                        e = np.array(e)
                        o = np.log2(e[:-1]/e[1:]) if len(e) > 1 else []
                        print('%-12s %-4s %-8s %s %-6s lev %s  L1 %s | ord %s' % (
                            case, arm, sch, mode, k.replace('m1_', ''), '-'.join(map(str, levs)),
                            ' '.join('%.2e' % x for x in e),
                            ' '.join('%5.2f' % x for x in o)))


if __name__ == '__main__':
    main()
