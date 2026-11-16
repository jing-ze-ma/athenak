#!/usr/bin/env python3
"""Gas-feedback tests of the clean cell layout on spherical polar (R v2 momentum form).
usage: fb.py run <bin> <case> <tag> <N> <n2> <n3> <mode> [key=val ...]   (mode sph | gd)
  case fthin: the T6a radiating sphere (bath 4 at r = 1, thin) through a uniform
              scattering gas (problem keys via env RHO, KS, TG; NLIM steps): the radiation
              force on the gas is chi F/c with r^2 F = L/(4 pi) constant (div F = 0 in a
              steady pure-scattering field); clean_init = steady starts it steady
  case flz:   latc (diffusion, reflecting r walls, theta walls, phi periodic; hydro r walls
              reflecting) with the gas in solid-body rotation gas_rot = OM: a closed box,
              the drag hands L_z to the radiation; total L_z (gas + radiation) conserved
  case flzm:  the same with a vacuum outer face (lat, e_out = 0, no bath; vet_gd)
       fb.py force <case> <tag> <N>   gas momentum vs t chi F/c (cell and analytic)
       fb.py hst <case> <tag> <N>     closed-box budgets from the am_hist history
"""
import glob
import subprocess
import os
import sys

import numpy as np

T = os.environ['REG'] + '/fb'
sys.path.insert(0, T)
import gt  # noqa: E402
import gtg  # noqa: E402

W = os.path.join(T, 'runsf')
gtg.W = W
BASE = {'fthin': 'thin', 'flz': 'latc', 'flzm': 'lat'}


def run(binary, case, tag, n, n2, n3, mode, extra):
    d = os.path.join(W, case, f'{tag}_{n}')
    os.makedirs(d, exist_ok=True)
    for f in glob.glob(os.path.join(d, '*.hst')):
        os.remove(f)
    p = os.path.join(d, 'sp.athinput')
    gtg.inp(BASE[case], n, n2, n3, mode, extra, p)
    s = open(p).read()
    if os.environ.get('FB', '1') == '1':
        s = s.replace('gas_feedback = false', 'gas_feedback = true')
    ks = os.environ.get('KS')
    if ks:
        i = s.index('kappa_s = ')
        s = s[:i] + f'kappa_s = {ks}' + s[s.index('\n', i):]
    s = s.replace('gas_rho = 1.0\n', f"gas_rho = {os.environ.get('RHO', '1.0')}\n")
    s = s.replace('gas_temp = 1.0e-16\n', f"gas_temp = {os.environ.get('TG', '1e-16')}\n")
    if 'NLIM' in os.environ:
        i = s.index('nlim = ')
        s = s[:i] + f"nlim = {os.environ['NLIM']}" + s[s.index('\n', i):]
        s = s.replace('tlim = 1e+30', 'tlim = 1e30')
    add = 'user_hist = true\nam_hist = true\n'
    if case == 'flzm':
        # vet_gd has no reflecting outer x1: outer vacuum, no radiation near that face
        i = s.index('implicit_ebath_x1max = ')
        s = s[:i] + 'implicit_ebath_x1max = 0.0' + s[s.index('\n', i):]
        i = s.index('e_out = ')
        s = s[:i] + 'e_out = 1.0e-30' + s[s.index('\n', i):]
    if case in ('flz', 'flzm'):
        add += f"gas_rot = {os.environ.get('OM', '2.0e-3')}\n"
        s = s.replace('ix1_bc = outflow', 'ix1_bc = reflect')
        s = s.replace('ox1_bc = outflow', 'ox1_bc = reflect')
    s = s.replace('width = ', add + 'width = ')
    dto = 1.0e30 if case == 'fthin' else gt.cfg('lat')['T']*(1.0 - 1.0e-9)
    s += f"""
<output7>
file_type = hst
data_format = %24.16e
dt = 1.0e-12

<output8>
file_type = bin
bin_double = true
variable = hydro_u
dt = {dto!r}
"""
    if 'TIMEKEY' in os.environ:
        s = s.replace('<time>\n', '<time>\n' + os.environ['TIMEKEY'] + '\n')
    if 'HYDKEY' in os.environ:
        s = s.replace('<hydro>\n', '<hydro>\n' + os.environ['HYDKEY'] + '\n')
    if 'OUTEVERY' in os.environ:
        import re
        s = re.sub(r'(<output[18]>[^<]*?)\ndt = [^\n]*', r'\1\ndt = 1.0e-12', s)
    open(p, 'w').write(s)
    with open(os.path.join(d, 'run.log'), 'w') as lg:
        pre = ['mpirun', '-np', os.environ['NP'], '--bind-to', 'none'] if 'NP' in os.environ else []
        subprocess.run(pre + [binary, '-i', 'sp.athinput'], cwd=d, stdout=lg,
                       stderr=subprocess.STDOUT)


def glue(b):
    """the MeshBlocks of a bin output joined into (n3, n2, n1) arrays"""
    lg = np.asarray(b['mb_logical'])[:, :3]
    nb = lg.max(axis=0) + 1
    out = {}
    for k, v in b['mb_data'].items():
        a = [np.asarray(x) for x in v]
        n3, n2, n1 = a[0].shape
        g = np.zeros((nb[2]*n3, nb[1]*n2, nb[0]*n1))
        for q, x in enumerate(a):
            i, j, kk = lg[q]
            g[kk*n3:(kk+1)*n3, j*n2:(j+1)*n2, i*n1:(i+1)*n1] = x
        out[k] = g
    return out


def hydro(d, first=False):
    sys.path.insert(0, os.environ['ATHENA_VIS'])
    import bin_convert
    fs = sorted(glob.glob(os.path.join(d, 'bin', '*.hydro_u.*.bin')))
    b = bin_convert.read_binary(fs[0 if first else -1])
    return b['time'], glue(b)


def m1(d):
    sys.path.insert(0, os.environ['ATHENA_VIS'])
    import bin_convert
    fs = sorted(glob.glob(os.path.join(d, 'bin', '*.m1.*.bin')))
    b = bin_convert.read_binary(fs[-1])
    return b['time'], glue(b)


def force(case, tag, n):
    """fthin: gas radial momentum per unit volume after t vs t chi F/c: against the
    cells' own final y_r (local force law) and against L/(4 pi r^2 c), L from the
    outermost cells' F r^2 and from the analytic F(r_in) = 1"""
    d = os.path.join(W, case, f'{tag}_{n}')
    t, h = hydro(d)
    tm, r = m1(d)
    p = open(os.path.join(d, 'sp.athinput')).read()
    val = lambda k: float(p.split(f'\n{k} = ')[1].split()[0])
    chi = val('kappa_s')*val('gas_rho')
    cl = val('c_light')
    n3, n2, n1 = h['mom1'].shape
    c, rv, v = gtg.geom(BASE[case], n1, n2, n3)
    mr = h['mom1']
    fy = r['m1_f1']
    loc = mr/(t*chi*fy/cl) - 1.0
    r2f = (fy*rv[None, None, :]**2).mean(axis=(0, 1))
    an = mr/(t*chi/(cl*rv[None, None, :]**2)) - 1.0
    tot = np.sum(mr*v)/np.sum(t*chi/(cl*rv[None, None, :]**2)*v) - 1.0
    lat = max(np.abs(h['mom2']).max(), np.abs(h['mom3']).max())/np.abs(mr).max()
    vmax = np.abs(mr/h['dens']).max()
    print('%s %s N %d (%dx%d) t %.4g chi %.1e: m_r/(t chi y_r)-1 max %.2e mean %.2e | '
          'vs 1/r^2: L1 %.2e total %.2e | r^2F/c in,out %.4f %.4f | lat %.1e vmax/c %.1e'
          % (case, tag, n1, n2, n3, t, chi, np.abs(loc).max(), loc.mean(),
             np.sum(np.abs(an)*v)/np.sum(v), tot, r2f[0]/cl, r2f[-1]/cl, lat, vmax/cl))


def hst(case, tag, n):
    d = os.path.join(W, case, f'{tag}_{n}')
    f = os.path.join(d, 'sp.user.hst')
    hd = open(f).read().split('\n')
    lab = [x.split('=')[1] for x in hd[1].replace('#', '').split()]
    a = np.atleast_2d(np.loadtxt(f))
    a = a[np.nonzero(a[:, 0] == 0.0)[0][-1]:]   # the last run in the file
    col = {k: a[:, q] for q, k in enumerate(lab)}
    lz = col['Lgz'] + col['Lrz']
    en = col['Eg'] + col['Er']
    sc = max(np.abs(col['Lrz']).max(), np.abs(col['Lgz']).max(), 1e-300)
    lx = col['Lgx'] + col['Lrx']
    print('  abs drifts: Lz %.2e (eps sum|L| %.1e)  E %.2e (eps E %.1e)' % (
        np.abs(lz - lz[0]).max(), 2.2e-16*np.abs(col['Lgz']).max(),
        np.abs(en - en[0]).max(), 2.2e-16*en[0]))
    print('%s %s N %d: steps %d t %.4g | Lz(gas+rad) drift max %.2e (of max|Lrz| %.2e, '
          'final Lgz %.3e Lrz %.3e) | Lx drift %.2e (of max|Lrx|) | E(gas+rad) drift %.2e '
          '(of max|dEr| %.2e)' %
          (case, tag, n, len(lz) - 1, col['time'][-1], np.abs(lz - lz[0]).max()/sc, sc,
           col['Lgz'][-1], col['Lrz'][-1],
           np.abs(lx - lx[0]).max()/max(np.abs(col['Lrx']).max(), 1e-300),
           np.abs(en - en[0]).max()/np.abs(col['Er'] - col['Er'][0]).max(),
           np.abs(col['Er'] - col['Er'][0]).max()))


def pr(case, tag, n):
    """radial momentum budget: change of sum V (m_r + y_r/c) per step (bin outputs every
    step) against the clean write-back's sum V (-h div_r/c) (DBGPR lines, one solve per
    step: time_scheme = be); exch = sum V (dm_r + dy_r/c + h div_r/c)"""
    sys.path.insert(0, os.environ['ATHENA_VIS'])
    import bin_convert
    d = os.path.join(W, case, f'{tag}_{n}')
    cl = float(open(os.path.join(d, 'sp.athinput')).read().split('\nc_light = ')[1].split()[0])
    hs = sorted(glob.glob(os.path.join(d, 'bin', '*.hydro_u.*.bin')))
    ms = sorted(glob.glob(os.path.join(d, 'bin', '*.m1.*.bin')))
    p = []
    for fh, fm in zip(hs, ms):
        bh, bm = bin_convert.read_binary(fh), bin_convert.read_binary(fm)
        gh, gm = glue(bh), glue(bm)
        n3, n2, n1 = gh['dens'].shape
        c, rv, v = gtg.geom(BASE[case], n1, n2, n3)
        p.append((bh['cycle'], np.sum(v*(gh['mom1'] + gm['m1_f1']/cl**2)),
                  np.sum(v*(np.abs(gh['mom1']) + np.abs(gm['m1_f1'])/cl**2))))
    lg = open(os.path.join(d, 'run.log')).read()
    lg = lg[lg.find('clean_init steady:') + 1:] if 'clean_init steady:' in lg else lg
    src, cur = {}, {}
    ex = []
    for ln in lg.split('\n'):
        if ln.startswith('DBGPR'):
            w = ln.split()
            cy = int(w[2])
            src[cy] = float(w[8])   # be: the last (accepted) solve of the cycle
            cur[cy] = abs(float(w[10]))
            ex.append(abs(float(w[6]) - float(w[8]))/max(float(w[12]), 1e-300))
    wsc = 2.0*np.pi/(gtg.PH1)   # bins hold the wedge, DBGPR the code volumes: same
    res, ab = [], []
    for q in range(1, len(p)):
        cy0, cy1 = p[q-1][0], p[q][0]
        s_ = sum(src.get(cy, 0.0) for cy in range(cy0, cy1))
        dp = p[q][1] - p[q-1][1]
        res.append((dp - s_)/max(abs(s_), 1e-300))
        ab.append(abs(dp - s_)/(2.2e-16*(p[q][2] + sum(2.0*cur.get(cy, 0.0) for cy in range(cy0, cy1))) + 1e-300))
    res = np.array(res)
    print('  max |dPr - src| / (eps (sum V(|m_r| + |y_r|/c) + 2|curvature|)) = %.2e' % max(ab))
    if os.environ.get("PRV"): print(np.array(ab))
    print('%s %s N %d: steps %d | max |dPr - src|/|src| %.2e (median %.2e) | exch/abs max '
          '%.2e | Pr %.4e, sum|m_r| %.4e' % (case, tag, n, len(res), np.abs(res).max(),
                                             np.median(np.abs(res)), max(ex), p[-1][1], p[-1][2]))


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6]),
            int(sys.argv[7]), sys.argv[8], sys.argv[9:])
    else:
        globals()[sys.argv[1]](sys.argv[2], sys.argv[3], int(sys.argv[4]))
