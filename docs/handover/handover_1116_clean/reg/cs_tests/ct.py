#!/usr/bin/env python3
"""Cubed-sphere tests of the clean cell layout (S1).
usage: ct.py deck <case> <dir> <nr> <nlat> <mblat> [key=val ...]
       ct.py score <case> <dir>      (pulse, dif: L1 vs exact; y_tan)
       ct.py cons <dir>              (sum E V, sum y V first/last output)
       ct.py uni <dir>               (well balance)
       ct.py selfconv <case> <dir_coarse...> <dir_fine>
cases: uni pulse dif lat (lat_axis via key), closed (reflect both, asymmetric pulse)
"""
import glob
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.environ['ATHENA_VIS'])

FRAMES = {0: ((0, 1, 0), (0, 0, 1), (1, 0, 0)), 1: ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
          2: ((0, -1, 0), (0, 0, 1), (-1, 0, 0)), 3: ((0, 1, 0), (-1, 0, 0), (0, 0, 1)),
          4: ((1, 0, 0), (0, 0, 1), (0, -1, 0)), 5: ((0, 1, 0), (1, 0, 0), (0, 0, -1))}


def cfg(case):
    c = dict(r=(1.0, 24.0), r0=10.0, w=0.3, eout=0.1, eamp=1.0, ks=1.0e-10, cl=1.0,
             clos='sphere', rs=0.0, bc=('reflect', 'marshak'), eb=(0.0, 0.0), extra={},
             T=4.0, cfl=0.5)
    if case == 'uni':
        c.update(r=(1.0, 3.0), eout=1.0, eamp=0.0, clos='eddington', bc=('reflect',
                 'reflect'), T=None, nlim=10, cfl=2.0)
    elif case == 'pulse':
        c.update(eb=(0.0, 0.4), win=(5.0, 20.0))
    elif case in ('dif', 'lat'):
        D = 10.0/(3.0*4.0e4)
        c.update(r=(1.0, 5.0), r0=3.0, w=1.0/math.sqrt(40.0), T=0.22**2/D,
                 ks=4.0e4, cl=10.0, rs=100.0, eout=math.exp(-10.0),
                 eb=(0.0, math.exp(-10.0)), steps8=64)
        if case == 'lat':
            c['extra'] = {'lat_amp': 0.5, 'lat_l': 2, 'lat_axis': 'z'}
            c.update(T=0.25*c['T'], steps8=16)
    elif case == 'thin':
        # T6a: black-body sphere R = 1 (bath 4 at r = 1, flux 1) in vacuum, steady
        c.update(r=(1.0, 10.0), r0=0.0, w=1.0, eamp=0.0, eout=1.0e-12, rs=1.0,
                 bc=('marshak', 'marshak'), eb=(4.0, 0.0), T=None, nlim=100, cfl=4.0)
    elif case == 'feed':
        # closed box (hydro and radiation reflect on both shells), gas feedback on: a
        # rotating, radially moving, laterally structured absorbing gas and an asymmetric
        # radiation field exchange energy, L_z and radial momentum
        c.update(r=(1.0, 3.0), r0=2.0, w=0.3, eout=1.0, eamp=1.0, ks=2.0, cl=10.0,
                 clos='eddington', bc=('reflect', 'reflect'), T=None, nlim=10, cfl=4.0,
                 hbc='reflect', fb='true', gtemp=0.5,
                 extra={'lat_amp': 0.5, 'lat_l': 1, 'lat_axis': 'vertex',
                        'gas_rot': 0.3, 'gas_vr': 0.2, 'rho_lat_amp': 0.4})
    elif case == 'closed':
        # asymmetric pulse far from the walls in a closed box (reflect both): energy and
        # Cartesian momentum change only by the walls, where E ~ 0
        c.update(r=(1.0, 20.0), r0=10.5, w=0.4, eout=0.0, clos='eddington', bc=('reflect',
                 'reflect'), T=1.0, extra={'lat_amp': 0.9, 'lat_l': 1,
                                           'lat_axis': 'vertex'})
    return c


def deck(case, d, nr, nlat, mblat, kv):
    c = cfg(case)
    os.makedirs(d, exist_ok=True)
    dr = (c['r'][1] - c['r'][0])/nr
    dxl = c['r'][0]*0.5*math.pi/nlat
    if 'steps8' in c:
        nst = c['steps8']*nr//8
        dt = c['T']/nst
    elif case == 'thin':
        dt = c['cfl']*dr/c['cl']
    elif c['T'] is None:
        dt = c['cfl']*min(dr, dxl)/c['cl']
    else:
        dt = c['cfl']*dr/c['cl']
    tlim = c['T'] if c['T'] is not None else 1.0e30
    nlim = c.get('nlim', -1)
    mbr = int(kv.pop('mbr', nr))
    extra = dict(c['extra'])
    rad = {}
    for k, v in kv.items():
        if k.startswith('rad.'):
            rad[k[4:]] = v
        else:
            extra[k] = v
    s = f"""<job>
basename = cs

<mesh>
nghost = 2
nx1 = {nr}
x1min = {c['r'][0]}
x1max = {c['r'][1]}
ix1_bc = {c.get('hbc', 'outflow')}
ox1_bc = {c.get('hbc', 'outflow')}
nx2 = {nlat}
x2min = -1.0
x2max = 1.0
ix2_bc = panel
ox2_bc = panel
nx3 = {nlat}
x3min = -1.0
x3max = 1.0
ix3_bc = panel
ox3_bc = panel
use_cubed_sphere = true

<meshblock>
nx1 = {mbr}
nx2 = {mblat}
nx3 = {mblat}

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
gas_feedback = {c.get('fb', 'false')}
transport = implicit
implicit_cfl = {dt*c['cl']/min(dr, dxl)!r}
implicit_solver = bicgstab
implicit_precond = mg
implicit_krylov_fuse = 1
implicit_recon = dc
implicit_bc_x1min = {c['bc'][0]}
implicit_bc_x1max = {c['bc'][1]}
implicit_ebath_x1min = {c['eb'][0]!r}
implicit_ebath_x1max = {c['eb'][1]!r}
time_scheme = hesdirk2
ic_rad = none
scheme = clean
clean_layout = cell
clean_closure = {c['clos']}
clean_sphere_radius = {c['rs']!r}
clean_tol = 1.0e-14
clean_tol_nl = 1.0e-13
clean_maxnl = 60
"""
    for k, v in rad.items():
        import re
        if re.search(r'(?m)^%s = ' % k, s):
            s = re.sub(r'(?m)^%s = .*$' % k, f"{k} = {v}", s)
        else:
            s += f"{k} = {v}\n"
    s += f"""
<problem>
pgen_name = rad_m1_beam
m1_test = sph_shell
gas_rho = 1.0
gas_temp = {c.get('gtemp', 1.0e-16)!r}
e_out = {c['eout']!r}
e_amp = {c['eamp']!r}
r0 = {c['r0']!r}
width = {c['w']!r}
"""
    for k, v in extra.items():
        s += f"{k} = {v}\n"
    s += f"""
<output1>
file_type = bin
bin_double = true
variable = m1
dt = {(c['T']*(1.0 - 1.0e-9) if c['T'] else 1.0e30)!r}

<output2>
file_type = bin
bin_double = true
variable = hydro_u
dt = {(c['T']*(1.0 - 1.0e-9) if c['T'] else 1.0e30)!r}
"""
    open(os.path.join(d, 'cs.athinput'), 'w').write(s)


def solid(x1, x2, y1, y2):
    f = lambda x, y: np.arctan(x*y/np.sqrt(1.0 + x*x + y*y))
    return f(x2, y2) - f(x1, y2) - f(x2, y1) + f(x1, y1)


def lat_vec(p, t, u0, u1, d):
    """angle * unit normal (+index) of a lateral face (as the code)"""
    a, b, n = (np.array(v, float) for v in FRAMES[p])
    def pt(x, y):
        return a[None, :]*x[:, None] + b[None, :]*y[:, None] + n[None, :]
    if d == 0:
        P, Q = pt(t, u1), pt(t, u0)
    else:
        P, Q = pt(u0, t), pt(u1, t)
    nv = np.cross(P, Q)
    s = np.linalg.norm(nv, axis=1)
    th = np.arctan2(s, np.sum(P*Q, axis=1))
    return (th/s)[:, None]*nv


def load(f, nbp=None):
    """per block: time, E, y (Cartesian, F/c), rhat, r (centroid radius), volume"""
    import bin_convert
    b = bin_convert.read_binary(f)
    md = b['mb_data']
    nmb = b['n_mbs']
    nbp = nbp or nmb//6
    g = b['mb_geometry']
    nx1, nx2, nx3 = b['nx1_mb'], b['nx2_mb'], b['nx3_mb']
    out = []
    for m in range(nmb):
        p = m//nbp
        x1a, x1b, x2a, x2b, x3a, x3b = g[m]
        rf = np.linspace(x1a, x1b, nx1 + 1)
        x2f = np.linspace(x2a, x2b, nx2 + 1)
        x3f = np.linspace(x3a, x3b, nx3 + 1)
        tx, ty = np.tan(0.25*np.pi*x2f), np.tan(0.25*np.pi*x3f)
        S = np.zeros((nx3, nx2, 3))
        om = np.zeros((nx3, nx2))
        for k in range(nx3):
            for j in range(nx2):
                one = np.ones(1)
                v = (lat_vec(p, tx[j+1]*one, ty[k]*one, ty[k+1]*one, 0)
                     - lat_vec(p, tx[j]*one, ty[k]*one, ty[k+1]*one, 0)
                     + lat_vec(p, ty[k+1]*one, tx[j]*one, tx[j+1]*one, 1)
                     - lat_vec(p, ty[k]*one, tx[j]*one, tx[j+1]*one, 1))
                S[k, j] = -0.5*v[0]
                om[k, j] = solid(tx[j], tx[j+1], ty[k], ty[k+1])
        rh = S/np.linalg.norm(S, axis=2)[:, :, None]
        vr = (rf[1:]**3 - rf[:-1]**3)/3.0
        rc = 0.75*(rf[1:]**4 - rf[:-1]**4)/(rf[1:]**3 - rf[:-1]**3)
        vol = om[:, :, None]*vr[None, None, :]
        E = np.asarray(md['m1_e'][m])
        F = np.stack([np.asarray(md['m1_f%d' % q][m]) for q in (1, 2, 3)], axis=-1)
        out.append(dict(E=E, F=F, rh=np.broadcast_to(rh[:, :, None, :], F.shape),
                        r=np.broadcast_to(rc, E.shape), vol=vol))
    return b['time'], out


def outs(d, v='m1'):
    return sorted(glob.glob(os.path.join(d, 'bin', '*.%s.*.bin' % v)))


def cat(bl, key):
    return np.concatenate([b[key].reshape(-1, *b[key].shape[3:]) for b in bl])


def exact(case, r, t):
    c = cfg(case)
    if case == 'thin':
        mc = np.sqrt(np.maximum(1.0 - 1.0/r**2, 0.0))
        return 2.0*(1.0 - mc), 1.0/r**2
    if case == 'pulse':
        def g(x):
            return x*x*(c['eout'] + c['eamp']*np.exp(-((x - c['r0'])/c['w'])**2))
        a, b = g(r - c['cl']*t), g(r + c['cl']*t)
        return 0.5*(a + b)/r**2, 0.5*(a - b)/r**2
    D = c['cl']/(3.0*c['ks'])
    al = 1.0/c['w']**2
    be = 1.0 + 4.0*al*D*t
    s = r - c['r0']
    u = np.exp(-al*s*s/be)*(c['r0']/np.sqrt(be) + s/be**1.5)
    du = np.exp(-al*s*s/be)*(1.0/be**1.5 - 2.0*al*s/be*(c['r0']/np.sqrt(be) +
                                                         s/be**1.5))
    dE = du/r - u/r**2
    return c['eout'] + u/r, -(1.0/(3.0*c['ks']))*dE


def score(case, d):
    c = cfg(case)
    t, bl = load(outs(d)[-1])
    E, F, rh, r, v = (cat(bl, k) for k in ('E', 'F', 'rh', 'r', 'vol'))
    y = F/c['cl']   # (y_r, y_theta, y_phi) in the frame about the hydro's r_hat
    yr = y[:, 0]
    yt = np.linalg.norm(y[:, 1:], axis=1)
    ee, ye = exact(case, r, t)
    sub = c['eout']
    if 'win' in c:
        v = v*((r > c['win'][0]) & (r < c['win'][1]))
        sub = c['eout']*(1.0 + (c['cl']*t/r)**2)
    l1e = np.sum(np.abs(E - ee)*v)/np.sum(np.abs(ee - sub)*v)
    l1y = np.sum(np.abs(yr - ye)*v)/np.sum(np.abs(ye)*v)
    # lateral spread of E per shell
    rs = np.unique(np.round(r, 12))
    spr = max(np.ptp(E[np.isclose(r, q)]) for q in rs)/np.abs(E).max()
    print('%s %s t %.6g L1(E) %.3e L1(y_r) %.3e max|y_t|/max|y_r| %.2e max shell spread '
          'E %.2e' % (case, d, t, l1e, l1y, yt.max()/np.abs(yr).max(), spr))
    return l1e, l1y


def cons(d):
    fs = outs(d)
    res = []
    for f in (fs[0], fs[-1]):
        t, bl = load(f)
        E, F, v = (cat(bl, k) for k in ('E', 'F', 'vol'))
        res.append((t, np.sum(E*v), np.sum(F*v[:, None], axis=0), np.sum(np.abs(F)*v[:, None])))
    (t0, e0, p0, a0), (t1, e1, p1, a1) = res
    print('%s t %.4g -> %.4g: sum E V %.17e -> %.17e rel %.2e | sum F V %s -> %s, '
          '|dP|/sum|F|V %.2e' % (d, t0, t1, e0, e1, e1/e0 - 1.0, p0, p1,
                                  np.abs(p1 - p0).max()/a1))


def frame_hydro(p, x2v, x3v):
    """rhat, e_xi, e_eta at the cell centres (as cubed_sphere.hpp)"""
    a, b, n = (np.array(v, float) for v in FRAMES[p])
    x, y = np.meshgrid(np.tan(0.25*np.pi*x2v), np.tan(0.25*np.pi*x3v))  # (k, j)
    dl = np.sqrt(1.0 + x*x + y*y)
    rh = (a*x[..., None] + b*y[..., None] + n)/dl[..., None]
    e1 = (a*dl[..., None] - rh*x[..., None])/np.sqrt(1.0 + y*y)[..., None]
    e2 = (b*dl[..., None] - rh*y[..., None])/np.sqrt(1.0 + x*x)[..., None]
    return rh, e1, e2


def budget(d, cl):
    """total energy, L_z (gas + radiation), radial momentum, first/last output"""
    import bin_convert
    res = []
    fm, fh = outs(d), outs(d, 'hydro_u')
    for q in (0, -1):
        bm, bh = bin_convert.read_binary(fm[q]), bin_convert.read_binary(fh[q])
        nmb = bm['n_mbs']
        tot = np.zeros(6)
        for m in range(nmb):
            p = m//(nmb//6)
            g = bm['mb_geometry'][m]
            nx1, nx2, nx3 = bm['nx1_mb'], bm['nx2_mb'], bm['nx3_mb']
            rf = np.linspace(g[0], g[1], nx1 + 1)
            x2f, x3f = np.linspace(g[2], g[3], nx2 + 1), np.linspace(g[4], g[5], nx3 + 1)
            x2v, x3v = 0.5*(x2f[1:] + x2f[:-1]), 0.5*(x3f[1:] + x3f[:-1])
            tx, ty = np.tan(0.25*np.pi*x2f), np.tan(0.25*np.pi*x3f)
            om = solid(tx[None, :-1], tx[None, 1:], ty[:-1, None], ty[1:, None])
            V = om[:, :, None]*((rf[1:]**3 - rf[:-1]**3)/3.0)[None, None, :]
            ql = rf[:-1]/rf[1:]
            r = 0.25*(ql*ql + 1.0)/((ql*ql + ql + 1.0)/3.0)*(rf[1:] + rf[:-1])
            rh, e1, e2 = frame_hydro(p, x2v, x3v)
            cc = np.sum(e1*e2, axis=-1)[..., None]
            md, mh = bm['mb_data'], bh['mb_data']
            E = np.asarray(md['m1_e'][m])
            Fr, Ft, Fp = (np.asarray(md['m1_f%d' % k][m]) for k in (1, 2, 3))
            m1, m2, m3 = (np.asarray(mh[k][m]) for k in ('mom1', 'mom2', 'mom3'))
            en = np.asarray(mh['ener'][m])
            u1, u2 = (m2 - cc*m3)/(1 - cc*cc), (m3 - cc*m2)/(1 - cc*cc)
            pc = m1[..., None]*rh[:, :, None, :] + u1[..., None]*e1[:, :, None, :] + \
                u2[..., None]*e2[:, :, None, :]
            R = r[None, None, :]*np.sqrt(rh[..., 0]**2 + rh[..., 1]**2)[..., None]
            lg = r[None, None, :]*(rh[:, :, None, 0]*pc[..., 1] - rh[:, :, None, 1]*pc[..., 0])
            tot += [np.sum(E*V), np.sum(en*V), np.sum(lg*V), np.sum(R*Fp*V)/cl**2,
                    np.sum(m1*V), np.sum(Fr*V)/cl**2]
        res.append(tot)
    a, b = res
    print('%s: E_rad+E_gas %.15e -> %.15e rel %.2e (E_rad %.4e -> %.4e)' %
          (d, a[0]+a[1], b[0]+b[1], (b[0]+b[1])/(a[0]+a[1]) - 1, a[0], b[0]))
    sl = abs(a[2]) + abs(b[2]) + abs(a[3]) + abs(b[3])
    print('   L_z gas %.6e -> %.6e, rad %.6e -> %.6e, total change/scale %.2e' %
          (a[2], b[2], a[3], b[3], ((b[2]+b[3]) - (a[2]+a[3]))/sl))
    print('   gas alone: E_gas rel %.2e, L_z gas rel %.2e' % (b[1]/a[1] - 1, b[2]/a[2] - 1))
    print('   radial momentum gas %.6e -> %.6e, rad %.6e -> %.6e' % (a[4], b[4], a[5], b[5]))


def uni(d):
    fs = outs(d)
    t, bl = load(fs[-1])
    E, F = cat(bl, 'E'), cat(bl, 'F')
    print('%s t %.4g max|E-1| %.2e max|F|/c %.2e' % (d, t, np.abs(E - 1.0).max(),
                                                   np.abs(F).max()))


def assemble(f):
    """E and V on the global (panel, k, j, i) grid (panel = block // blocks per panel)"""
    import bin_convert
    t, bl = load(f)
    b = bin_convert.read_binary(f)
    nmb, g = b['n_mbs'], b['mb_geometry']
    n1, n2, n3 = b['Nx1'], b['Nx2'], b['Nx3']
    x1a, x1b = b['x1min'], b['x1max']
    E = np.zeros((6, n3, n2, n1))
    V = np.zeros_like(E)
    for m in range(nmb):
        p = m//(nmb//6)
        e = bl[m]['E']
        k0 = int(round((g[m][4] + 1.0)/2.0*n3))
        j0 = int(round((g[m][2] + 1.0)/2.0*n2))
        i0 = int(round((g[m][0] - x1a)/(x1b - x1a)*n1))
        sl = (p, slice(k0, k0 + e.shape[0]), slice(j0, j0 + e.shape[1]),
              slice(i0, i0 + e.shape[2]))
        E[sl] = e
        V[sl] = bl[m]['vol']
    return E, V


def selfconv(case, ds):
    ef, vf = assemble(outs(ds[-1])[-1])
    res = []
    for d in ds[:-1]:
        ec, vc = assemble(outs(d)[-1])
        f = ef.shape[1]//ec.shape[1]
        fr = ef.shape[3]//ec.shape[3]
        _, n3, n2, n1 = ec.shape
        er = (ef*vf).reshape(6, n3, f, n2, f, n1, fr).sum(axis=(2, 4, 6))
        vr = vf.reshape(6, n3, f, n2, f, n1, fr).sum(axis=(2, 4, 6))
        er = er/vr
        res.append(np.sum(np.abs(ec - er)*vc)/np.sum(np.abs(er - cfg(case)['eout'])*vc))
        o = '' if len(res) == 1 else ' order %.2f' % math.log2(res[-2]/res[-1])
        print('%s vs %s: L1(E) %.3e%s' % (d, ds[-1], res[-1], o))

if __name__ == '__main__':
    a = sys.argv
    if a[1] == 'deck':
        kv = dict(x.split('=', 1) for x in a[7:])
        deck(a[2], a[3], int(a[4]), int(a[5]), int(a[6]), kv)
    elif a[1] == 'score':
        score(a[2], a[3])
    elif a[1] == 'cons':
        cons(a[2])
    elif a[1] == 'budget':
        budget(a[2], float(a[3]))
    elif a[1] == 'uni':
        uni(a[2])
    elif a[1] == 'selfconv':
        selfconv(a[2], a[3:])
