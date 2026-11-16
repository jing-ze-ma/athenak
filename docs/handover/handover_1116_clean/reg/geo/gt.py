#!/usr/bin/env python3
"""Spherical-polar tests of the clean cell layout.
usage: gt.py run <bin> <case> <tag> <N...>   (runs on the login node, nice)
       gt.py score <case> <tag> <N...>
cases: pulse   1-D free streaming (sphere closure R = 0, D_rr = 1): E0 = exp(-((r-6)/0.4)^2),
               F = 0; exact r^2 E = (G(r-ct) + G(r+ct))/2, r^2 y = (G(r-ct) - G(r+ct))/2
       pulse2  the same on a 2-D r-theta wedge (nx2 = N/8)
       dif     T4 spherical diffusion (rho kappa_s = 4e4, c = 10, D = c/(3 chi)):
               E0 = e_out + exp(-40 (r-3)^2), u = r E obeys u_t = D u_rr
       thin    T6a black-body sphere R = 1 (flux 1, bath 4 at r = 1), thin: steady
               E = 2 (1 - sqrt(1 - 1/r^2)), F = 1/r^2
       lat     dif with a lateral mode (1 + 0.5 cos(pi (theta - th0)/dth)), 2-D (nx2 = N/4),
               self-convergence
"""
import glob
import math
import os
import subprocess
import sys

import numpy as np

W = os.environ['REG'] + '/geo/runs'
TH0, TH1, PH1 = 1.1707963267948966, 1.9707963267948966, 1.0471975511965976


def cfg(case):
    c = dict(n2=lambda n: 1, eout=1.0e-12, eamp=1.0, ks=1.0e-10, cl=1.0, rs=0.0,
             bc1=('reflect', 'marshak'), eb=(0.0, 0.0), extra={})
    if case in ('pulse', 'pulse2'):
        # a background e_out (positivity of the tails); scored on the domain of
        # dependence [r_in + cT, r_out - cT] of the open ends
        c.update(r=(1.0, 24.0), r0=10.0, w=0.3, T=4.0, cfl=0.5, eout=0.1, eb=(0.0, 0.4),
                 win=(5.0, 20.0))
        if case == 'pulse2':
            c['n2'] = lambda n: n//8
    elif case in ('dif', 'lat', 'difc', 'latc'):
        D = 10.0/(3.0*4.0e4)
        c.update(r=(1.0, 5.0), r0=3.0, w=1.0/math.sqrt(40.0), T=0.22**2/D, ks=4.0e4,
                 cl=10.0, rs=100.0, eout=math.exp(-10.0), eb=(0.0, math.exp(-10.0)),
                 dt64=0.22**2/D/512)
        if case in ('difc', 'latc'):
            c['bc1'] = ('reflect', 'reflect')
        if case in ('lat', 'latc'):
            c['n2'] = lambda n: n//4
            c['extra'] = {'lat_amp': 0.5, 'lat_l2': 1.0}
    elif case == 'thin':
        c.update(r=(1.0, 10.0), r0=0.0, w=1.0, eamp=0.0, eout=1.0e-12, rs=1.0,
                 bc1=('marshak', 'marshak'), eb=(4.0, 0.0), steady=True)
    return c


def inp(case, n, path):
    c = cfg(case)
    n2 = c['n2'](n)
    dr = (c['r'][1] - c['r'][0])/n
    if c.get('steady'):
        dt = 4.0*dr/c['cl']
        tlim, nlim = 1.0e30, 400
    elif 'dt64' in c:
        dt = c['dt64']*64.0/n
        tlim, nlim = c['T'], -1
    else:
        dt = c['cfl']*dr/c['cl']
        tlim, nlim = c['T'], -1
    # implicit_cfl is relative to the smallest cell width (r_in dtheta on a 2-D wedge)
    dxmin = dr if n2 == 1 else min(dr, c['r'][0]*(TH1 - TH0)/n2)
    s = f"""<job>
basename = sp

<mesh>
nghost = 2
nx1 = {n}
x1min = {c['r'][0]}
x1max = {c['r'][1]}
ix1_bc = outflow
ox1_bc = outflow
nx2 = {n2}
x2min = {TH0!r}
x2max = {TH1!r}
ix2_bc = reflect
ox2_bc = reflect
nx3 = 1
x3min = 0.0
x3max = {PH1!r}
ix3_bc = periodic
ox3_bc = periodic
use_spherical_polar = true

<meshblock>
nx1 = {n}
nx2 = {n2}
nx3 = 1

<time>
evolution = dynamic
integrator = rk2
cfl_number = 0.3
nlim = {nlim}
tlim = {tlim!r}
ndiag = 1000000

<hydro>
eos = ideal
reconstruct = plm
rsolver = hllc
gamma = 1.666666666666667
rad_signal_speed = false

<rad_m1>
c_light = {c['cl']}
chat_over_c = 1.0
closure = eddington
e_floor = 1.0e-30
opacity = const
kappa_p = 0.0
kappa_e = 0.0
kappa_f = 0.0
kappa_s = {c['ks']!r}
arad = 1.0
coupling = true
gas_feedback = false
transport = implicit
implicit_cfl = {dt*c['cl']/dxmin!r}
implicit_solver = bicgstab
implicit_precond = line
implicit_recon = dc
implicit_bc_x1min = {c['bc1'][0]}
implicit_bc_x1max = {c['bc1'][1]}
implicit_ebath_x1min = {c['eb'][0]!r}
implicit_ebath_x1max = {c['eb'][1]!r}
time_scheme = hesdirk2
ic_rad = none
scheme = clean
clean_layout = cell
clean_closure = sphere
clean_sphere_radius = {c['rs']!r}
clean_tol = 1.0e-14
clean_tol_nl = 1.0e-13
clean_maxnl = 60
{os.environ.get('GT_RAD', '').replace(',', chr(10))}

<problem>
pgen_name = rad_m1_beam
m1_test = sph_shell
gas_rho = 1.0
gas_temp = 1.0e-16
e_out = {c['eout']!r}
e_amp = {c['eamp']!r}
r0 = {c['r0']!r}
width = {c['w']!r}
"""
    for k, v in c['extra'].items():
        s += f"{k} = {v}\n"
    s += f"""
<output1>
file_type = {'tab' if n2 == 1 else 'bin'}
variable = m1
data_format = %24.16e
dt = {1.0e30 if c.get('steady') else c['T']*(1.0 - 1.0e-9)!r}
"""
    s += f"""
<output2>
file_type = bin
bin_double = true
variable = m1
dt = {1.0e30 if c.get('steady') else c['T']*(1.0 - 1.0e-9)!r}
"""
    if n2 > 1:
        s = s[:s.index('<output2>')]
        s = s.replace('variable = m1\n', 'bin_double = true\nvariable = m1\n', 1)
    open(path, 'w').write(s)


def run(binary, case, tag, ns):
    for n in ns:
        d = os.path.join(W, case, f'{tag}_{n}')
        os.makedirs(d, exist_ok=True)
        inp(case, n, os.path.join(d, 'sp.athinput'))
        with open(os.path.join(d, 'run.log'), 'w') as lg:
            subprocess.run(['nice', '-n', '10', binary, '-i', 'sp.athinput'], cwd=d,
                           stdout=lg, stderr=subprocess.STDOUT)


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(sys.argv[2], sys.argv[3], sys.argv[4], [int(x) for x in sys.argv[5:]])


def exact(case, r, t):
    """exact E, F/c at radius r and time t"""
    c = cfg(case)
    if case in ('pulse', 'pulse2'):
        def g(x):
            return x*x*(c['eout'] + c['eamp']*np.exp(-((x - c['r0'])/c['w'])**2))
        a, b = g(r - c['cl']*t), g(r + c['cl']*t)
        return 0.5*(a + b)/r**2, 0.5*(a - b)/r**2
    if case in ('dif', 'lat', 'difc', 'latc'):
        D = c['cl']/(3.0*c['ks'])
        al = 1.0/c['w']**2
        be = 1.0 + 4.0*al*D*t
        s = r - c['r0']
        u = np.exp(-al*s*s/be)*(c['r0']/np.sqrt(be) + s/be**1.5)
        # F = -(c/(3 chi)) dE/dr, as y = F/c
        du = np.exp(-al*s*s/be)*(1.0/be**1.5 - 2.0*al*s/be*(c['r0']/np.sqrt(be) +
                                                             s/be**1.5))
        dE = du/r - u/r**2
        return c['eout'] + u/r, -(1.0/(3.0*c['ks']))*dE
    if case == 'thin':
        mc = np.sqrt(np.maximum(1.0 - 1.0/r**2, 0.0))
        return 2.0*(1.0 - mc), 1.0/r**2


def faces(case, n):
    c = cfg(case)
    return np.linspace(c['r'][0], c['r'][1], n + 1)


def tab(d):
    f = sorted(glob.glob(os.path.join(d, 'tab', '*.tab')))[-1]
    t = float(open(f).readline().split('time=')[1].split()[0])
    a = np.loadtxt(f)
    return t, {'x1v': a[:, 2], 'E': a[:, -4], 'F1': a[:, -3], 'F2': a[:, -2]}


def score(case, tag, ns):
    c = cfg(case)
    rows = []
    for n in ns:
        d = os.path.join(W, case, f'{tag}_{n}')
        t, v = tab(d)
        rf = faces(case, n)
        vol = (rf[1:]**3 - rf[:-1]**3)/3.0
        r = v['x1v']
        e = v['E'] if 'E' in v else v['r00']
        y = v['F1']/c['cl'] if 'F1' in v else v['r01']/c['cl']
        ee, ye = exact(case, r, t)
        sub = c['eout'] if case in ('dif', 'lat') else 0.0
        if 'win' in c:
            vol = vol*((r > c['win'][0]) & (r < c['win'][1]))
            sub = exact(case, r, 0.0)[0]*0.0 + c['eout']*(1.0 + (c['cl']*t/r)**2)
        l1e = np.sum(np.abs(e - ee)*vol)/np.sum(np.abs(ee - sub)*vol)
        l1y = np.sum(np.abs(y - ye)*vol)/np.sum(np.abs(ye)*vol)
        rows.append((n, t, l1e, l1y))
    print(f'{case} {tag}: N  t  L1(E)  L1(F)  orders')
    for q, (n, t, a, b) in enumerate(rows):
        o = ''
        if q > 0:
            o = '%.2f %.2f' % (math.log2(rows[q-1][2]/a), math.log2(rows[q-1][3]/b))
        print('%5d %.6g %.3e %.3e %s' % (n, t, a, b, o))


if __name__ == '__main__' and sys.argv[1] == 'score':
    score(sys.argv[2], sys.argv[3], [int(x) for x in sys.argv[4:]])


def energy(case, tag, n):
    """total E V (wedge factor dropped) of the first and last tab outputs"""
    d = os.path.join(W, case, f'{tag}_{n}')
    fs = sorted(glob.glob(os.path.join(d, 'tab', '*.tab')))
    c = cfg(case)
    n2 = c['n2'](n)
    r = np.linspace(c['r'][0], c['r'][1], n + 1)
    th = np.linspace(TH0, TH1, n2 + 1)
    vol = np.outer(np.cos(th[:-1]) - np.cos(th[1:]), (r[1:]**3 - r[:-1]**3)/3.0).ravel()
    out = []
    for f in (fs[0], fs[-1]):
        a = np.loadtxt(f)
        out.append(np.sum(a[:, -4]*vol))
    print(case, tag, n, 'sum E V: t0 %.17e t1 %.17e rel %.3e' % (out[0], out[1],
          out[1]/out[0] - 1.0))


if __name__ == '__main__' and sys.argv[1] == 'energy':
    energy(sys.argv[2], sys.argv[3], int(sys.argv[4]))


def binall(d):
    """(time, x1v, x2v, E, y1, y2) of every bin output, double precision"""
    sys.path.insert(0, os.environ['ATHENA_VIS'])
    import bin_convert
    out = []
    for f in sorted(glob.glob(os.path.join(d, 'bin', '*.bin'))):
        b = bin_convert.read_binary(f)
        md = b['mb_data']
        out.append((b['time'], np.asarray(md['m1_e'][0][0]),
                    np.asarray(md['m1_f1'][0][0]), np.asarray(md['m1_f2'][0][0])))
    return out


def vol2(case, n):
    c = cfg(case)
    r = np.linspace(c['r'][0], c['r'][1], n + 1)
    th = np.linspace(TH0, TH1, c['n2'](n) + 1)
    return np.outer(np.cos(th[:-1]) - np.cos(th[1:]), (r[1:]**3 - r[:-1]**3)/3.0)


def energy2(case, tag, n):
    a = binall(os.path.join(W, case, f'{tag}_{n}'))
    v = vol2(case, n)
    e0, e1 = np.sum(a[0][1]*v), np.sum(a[-1][1]*v)
    print(case, tag, n, 'sum E V: t0 %.17e t1 %.17e rel %.3e' % (e0, e1, e1/e0 - 1.0))


def selfconv(case, tag, ns):
    """L1 of E against the finest level, restricted by 2x2 volume averaging"""
    c = cfg(case)
    fin = binall(os.path.join(W, case, f'{tag}_{ns[-1]}'))[-1]
    vf = vol2(case, ns[-1])
    rows = []
    for n in ns[:-1]:
        a = binall(os.path.join(W, case, f'{tag}_{n}'))[-1]
        f = ns[-1]//n
        ef = fin[1]*vf
        n2, n1 = a[1].shape
        ec = ef.reshape(n2, f, n1, f).sum(axis=(1, 3))
        vc = vf.reshape(n2, f, n1, f).sum(axis=(1, 3))
        ec = ec/vc
        l1 = np.sum(np.abs(a[1] - ec)*vc)/np.sum(np.abs(ec - c['eout'])*vc)
        rows.append((n, l1))
    for q, (n, l1) in enumerate(rows):
        o = '' if q == 0 else '%.2f' % math.log2(rows[q-1][1]/l1)
        print('%s %s %5d vs %d: L1(E) %.3e %s' % (case, tag, n, ns[-1], l1, o))


def score2(case, tag, n):
    """2-D radial case: L1 vs the exact solution and max |y_theta|/max|y_r|"""
    c = cfg(case)
    a = binall(os.path.join(W, case, f'{tag}_{n}'))[-1]
    rf = np.linspace(c['r'][0], c['r'][1], n + 1)
    rv = 0.75*(rf[1:]**4 - rf[:-1]**4)/(rf[1:]**3 - rf[:-1]**3)
    ee, ye = exact(case, rv, a[0])
    v = vol2(case, n)*((rv > c['win'][0]) & (rv < c['win'][1]))
    sub = c['eout']*(1.0 + (c['cl']*a[0]/rv)**2)
    l1e = np.sum(np.abs(a[1] - ee)*v)/np.sum(np.abs(ee - sub)*v)
    l1y = np.sum(np.abs(a[2]/c['cl'] - ye)*v)/np.sum(np.abs(ye)*v)
    print('%s %s %d t %.6g L1(E) %.3e L1(F) %.3e max|F_th|/max|F_r| %.2e '
          'max_j spread of E %.2e' % (case, tag, n, a[0], l1e, l1y,
                                      np.abs(a[3]).max()/np.abs(a[2]).max(),
                                      np.abs(a[1] - a[1][:1]).max()/np.abs(a[1]).max()))


if __name__ == '__main__' and sys.argv[1] in ('energy2', 'selfconv', 'score2'):
    if sys.argv[1] == 'selfconv':
        selfconv(sys.argv[2], sys.argv[3], [int(x) for x in sys.argv[4:]])
    else:
        globals()[sys.argv[1]](sys.argv[2], sys.argv[3], int(sys.argv[4]))
