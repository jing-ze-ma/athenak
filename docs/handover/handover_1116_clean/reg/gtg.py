#!/usr/bin/env python3
"""vet_gd closure tests of the clean cell layout on a 3-D spherical-polar wedge.
usage: gtg.py run <bin> <case> <tag> <N> <n2> <n3> <mode> [key=val ...]
  case: thin | dif | pulse0 (pulse, E_bath 0) | lat ; mode: sph | gd
  phi period pi/2 (HEALPix permutation-exact), theta as gt.py; extra keys -> <rad_m1>
       gtg.py score <case> <tag> <N>       (shell-averaged L1 vs exact, lateral spread)
       gtg.py clos <tag_gd> <tag_sph> <case> <N>  (closure dumps: per-component errors)
"""
import glob
import math
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.environ['REG'] + '/geo')
import gt  # noqa: E402

W = os.environ['REG'] + '/runs'
PH1 = 0.5*math.pi
FULL = bool(os.environ.get('FULL'))
TH0, TH1 = (0.0, math.pi) if FULL else (gt.TH0, gt.TH1)
if FULL:
    PH1 = 2.0*math.pi


def base(case):
    return 'pulse' if case == 'pulse0' else case


def inp(case, n, n2, n3, mode, extra, path):
    bc = base(case)
    c0 = gt.cfg
    def cfg(cs):
        c = c0(cs)
        c['n2'] = lambda nn: n2
        if case == 'pulse0':
            c['eb'] = (0.0, 0.0)
        return c
    gt.cfg = cfg
    gt.inp(bc, n, path)
    gt.cfg = c0
    s = open(path).read()
    s = s.replace('nx3 = 1\n', f'nx3 = {n3}\n')
    if FULL:
        s = s.replace(f'x2min = {gt.TH0!r}', 'x2min = 0.0')
        s = s.replace(f'x2max = {gt.TH1!r}', f'x2max = {math.pi!r}')
        s = s.replace('ix2_bc = reflect\n', 'ix2_bc = polar\n')
        s = s.replace('ox2_bc = reflect\n', 'ox2_bc = polar\n')
        s = s.replace('use_spherical_polar = true\n',
                      'use_spherical_polar = true\nuse_polar_boundary = true\n')
    if 'MB3' in os.environ:
        i = s.index('<meshblock>')
        s = s[:i] + s[i:].replace(f'nx3 = {n3}\n', f"nx3 = {os.environ['MB3']}\n", 1)
    if 'MB2' in os.environ:
        i = s.index('<meshblock>')
        s = s[:i] + s[i:].replace(f'nx2 = {n2}\n', f"nx2 = {os.environ['MB2']}\n", 1)
    if os.environ.get('THPER'):
        s = s.replace('ix2_bc = reflect\n', 'ix2_bc = periodic\n')
        s = s.replace('ox2_bc = reflect\n', 'ox2_bc = periodic\n')
    if 'EOUT' in os.environ:
        s = s.replace('e_out = ', f"e_out = {os.environ['EOUT']}\n# e_out = ")
    for kv in os.environ.get('SUBS', '').split(';'):
        if '=>' in kv:
            s = s.replace(kv.split('=>')[0] + '\n', kv.split('=>')[1] + '\n')
    if 'NLIM' in os.environ:
        s = s.replace('nlim = 400\n', f"nlim = {os.environ['NLIM']}\n")
    s = s.replace(f'x3max = {gt.PH1!r}', f'x3max = {PH1!r}')
    if mode == 'gd':
        s = s.replace('closure = eddington\n', 'closure = vet_col\n')
        s = s.replace('clean_closure = sphere\n', 'clean_closure = extrap1\n')
        s = s.replace('scheme = clean\n', 'scheme = clean\nvet_col_lat = true\nvet_gd = true\n'
                      'vet_gd_twin = false\nvet_gd_rotate_dt = 0\n')
        # physical theta walls by default (THPER=1: the artificial theta periodicity)
        if not (os.environ.get('THPER') or FULL or any('theta_wall' in kv for kv in extra)):
            s = s.replace('scheme = clean\n', 'scheme = clean\nvet_gd_theta_wall = vacuum\n')
    add = ''.join(f'{kv.split("=", 1)[0]} = {kv.split("=", 1)[1]}\n' for kv in extra)
    s = s.replace('scheme = clean\n', 'scheme = clean\n' + add)
    open(path, 'w').write(s)


def run(binary, case, tag, n, n2, n3, mode, extra):
    d = os.path.join(W, case, f'{tag}_{n}')
    os.makedirs(d, exist_ok=True)
    inp(case, n, n2, n3, mode, extra, os.path.join(d, 'sp.athinput'))
    env = dict(os.environ, CLN_DUMP=os.path.join(d, 'cln.dump'))
    with open(os.path.join(d, 'run.log'), 'w') as lg:
        subprocess.run([binary, '-i', 'sp.athinput'], cwd=d, stdout=lg,
                       stderr=subprocess.STDOUT, env=env)


def binlast(d, first=False):
    sys.path.insert(0, os.environ['ATHENA_VIS'])
    import bin_convert
    fs = sorted(glob.glob(os.path.join(d, 'bin', '*.bin')))
    b = bin_convert.read_binary(fs[0 if first else -1])
    md = b['mb_data']
    return (b['time'], np.asarray(md['m1_e'][0]), np.asarray(md['m1_f1'][0]),
            np.asarray(md['m1_f2'][0]), np.asarray(md['m1_f3'][0]))


def geom(case, n, n2, n3):
    c = gt.cfg(base(case))
    rf = np.linspace(c['r'][0], c['r'][1], n + 1)
    th = np.linspace(TH0, TH1, n2 + 1)
    rv = 0.75*(rf[1:]**4 - rf[:-1]**4)/(rf[1:]**3 - rf[:-1]**3)
    v = ((np.cos(th[:-1]) - np.cos(th[1:]))[None, :, None]*
         ((rf[1:]**3 - rf[:-1]**3)/3.0)[None, None, :]*np.ones((n3, 1, 1))*PH1/n3)
    return c, rv, v


def score(case, tag, n):
    d = os.path.join(W, case, f'{tag}_{n}')
    t, e, f1, f2, f3 = binlast(d)
    n3, n2, n1 = e.shape
    c, rv, v = geom(case, n1, n2, n3)
    bc = base(case)
    ee, ye = gt.exact(bc, rv, t)
    ee, ye = ee[None, None, :], ye[None, None, :]
    w = v.copy()
    sub = 0.0
    if 'win' in c:
        w = w*((rv > c['win'][0]) & (rv < c['win'][1]))[None, None, :]
        sub = c['eout']*(1.0 + (c['cl']*t/rv)**2)[None, None, :]
    elif bc in ('dif', 'lat'):
        sub = c['eout']
    l1e = np.sum(np.abs(e - ee)*w)/np.sum(np.abs(ee - sub)*w)
    l1y = np.sum(np.abs(f1/c['cl'] - ye)*w)/np.sum(np.abs(ye)*w)
    spr = np.abs(e - e.mean(axis=(0, 1))[None, None, :]).max()/np.abs(e).max()
    lat = max(np.abs(f2).max(), np.abs(f3).max())/np.abs(f1).max()
    e0 = np.sum(binlast(d, True)[1]*v)
    print('%s %s N %d (%dx%d) t %.6g L1(E) %.3e L1(F) %.3e lat-spread(E) %.2e '
          'max|F_lat|/max|F_r| %.2e  sumEV rel %.3e' % (case, tag, n1, n2, n3, t, l1e, l1y,
                                                      spr, lat, np.sum(e*v)/e0 - 1.0))


def dump(d):
    a = np.fromfile(os.path.join(d, 'cln.dump'), dtype=np.int64, count=5)
    x = np.fromfile(os.path.join(d, 'cln.dump'), dtype=np.float64, offset=40)
    x = x.reshape(a)[0]
    g = 2
    return x[:, g:-g, g:-g, g:-g]


NAMES = (['D_rr', 'D_tt', 'D_pp', 'D_rt', 'D_rp', 'D_tp', 'h+_r', 'h+_t', 'h+_p', 'h-_r',
          'h-_t', 'h-_p'] + ['k+%d%d' % (dd, b) for dd in range(3) for b in range(3)] +
         ['C1', 'J'])


def sphere(rv, rs=1.0):
    """analytic sphere closure (CleanClosureSphere) at radii rv, CLN_* layout"""
    q = rs/rv
    mc = np.sqrt(np.maximum(1.0 - q*q, 0.0))
    c = np.zeros((23, rv.size))
    drr = (1.0 + mc + mc*mc)/3.0
    ht = 0.5*(np.arcsin(q) - mc*q)*(1.0 + mc)/(np.pi*q*q)
    kr = q*(1.0 + mc)/(3.0*np.pi)
    c[0], c[1], c[2] = drr, 0.5*(1.0 - drr), 0.5*(1.0 - drr)
    c[6], c[7], c[8] = 0.5*(1.0 + mc), ht, ht
    c[10], c[11] = -ht, -ht
    c[12], c[15], c[16], c[18], c[20] = drr, kr, 0.25*(1.0 - drr), kr, 0.25*(1.0 - drr)
    c[21] = 1.0
    c[22] = 2.0*(1.0 - mc)
    return c


def clos(tgd, case, n, rmax=1e30):
    a = dump(os.path.join(W, case, f'{tgd}_{n}'))
    n3, n2, n1 = a.shape[1:]
    c, rv, v = geom(case, n1, n2, n3)
    b = sphere(rv)[:, None, None, :]*np.ones((1, n3, n2, 1))
    sel = (rv <= rmax)[None, None, :]*np.ones((n3, n2, 1), bool)
    jin = np.zeros((n3, n2, n1), bool)
    jin[:, n2//4:n2 - n2//4, :] = True
    out = {}
    for q, nm in enumerate(NAMES):
        if nm == 'C1':
            continue
        da = a[q] - b[q]
        if nm == 'J':
            da = da/b[q]
        out[nm] = (np.abs(da[sel]).max(), np.sum((np.abs(da)*v)[sel])/np.sum(v[sel]),
                   np.abs(da[sel & jin]).max())
    return out


def clostab(tags, case, n, rmax, comps):
    print(f'{case} N {n} r <= {rmax}: L1 (volume mean |err|), max; J relative')
    print('%-14s' % 'tag' + ''.join('%-20s' % c for c in comps))
    for t in tags:
        o = clos(t, case, n, rmax)
        print('%-14s' % t + ''.join('%.2e/%.2e     ' % (o[c][1], o[c][0]) for c in comps))


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6]),
            int(sys.argv[7]), sys.argv[8], sys.argv[9:])
    elif sys.argv[1] == 'score':
        score(sys.argv[2], sys.argv[3], int(sys.argv[4]))
    elif sys.argv[1] == 'clos':
        clostab(sys.argv[5].split(','), sys.argv[2], int(sys.argv[3]), float(sys.argv[4]),
                ['J', 'D_rr', 'D_rt', 'D_rp', 'h+_r', 'h+_t', 'h-_t', 'k+10', 'k+11'])
