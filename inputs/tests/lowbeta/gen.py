#!/usr/bin/env python3
"""lowbeta-tests-1102: generate decks + task lists for the low-beta Cartesian MHD study.
usage: gen.py <test>   (cpaw2d cpaw1d lw vortex rd rotor blast)
Writes runs/<test>/<case>/run.athinput + meta.json and appends to tasks_<test>.txt."""
import json, math, os, sys

B = '/resnick/groups/carnegie_poc/jingze/lowbeta_1102'
BIP = B + '/athena_bip'
LBT = B + '/athena_lbt'

# scheme label -> <mhd> keys (all existing keys; nothing new)
SCHEMES = {
    'hlle_gs':  dict(rsolver='hlle'),
    'hlle_sm':  dict(rsolver='hlle', emf_smooth='1.0e-4'),
    'hlld_gs':  dict(rsolver='hlld'),
    'hlld_sm':  dict(rsolver='hlld', emf_smooth='1.0e-4'),
    'hlld_uct': dict(rsolver='hlld', emf_scheme='uct_hlld'),
    'lhlld_gs': dict(rsolver='lhlld', lhlld_cb_cont='true'),
    'lhlld_sm': dict(rsolver='lhlld', lhlld_cb_cont='true', emf_smooth='1.0e-4'),
    'lhlld_uct': dict(rsolver='lhlld', lhlld_cb_cont='true', emf_scheme='uct_hlld'),
}
S1D = ['hlle_gs', 'hlld_gs', 'lhlld_gs']   # no corner EMF in 1-D


def deck(mesh, time, mhd, prob, outs):
    s = '<job>\nbasename = run\n\n<mesh>\nnghost = 4\n'
    for k, v in mesh.items():
        s += f'{k} = {v}\n'
    s += '\n<meshblock>\n'
    for d in ('nx1', 'nx2', 'nx3'):
        s += f'{d} = {mesh.get(d, 1)}\n'
    s += '\n<time>\nevolution = dynamic\nintegrator = rk2\nnlim = -1\nndiag = 100000\n'
    for k, v in time.items():
        s += f'{k} = {v}\n'
    s += '\n<mhd>\neos = ideal\n'
    for k, v in mhd.items():
        s += f'{k} = {v}\n'
    s += '\n<problem>\n'
    for k, v in prob.items():
        s += f'{k} = {v}\n'
    for n, o in enumerate(outs, 1):
        s += f'\n<output{n}>\n'
        for k, v in o.items():
            s += f'{k} = {v}\n'
    return s


def mesh2d(n1, n2, x1, x2, bc='periodic', n3=1, x3=(-0.5, 0.5)):
    m = dict(nx1=n1, x1min=x1[0], x1max=x1[1], ix1_bc=bc, ox1_bc=bc,
             nx2=n2, x2min=x2[0], x2max=x2[1], ix2_bc=bc, ox2_bc=bc,
             nx3=n3, x3min=x3[0], x3max=x3[1], ix3_bc=bc, ox3_bc=bc)
    return m


def mesh1d(n1, x1, bc='periodic'):
    return dict(nx1=n1, x1min=x1[0], x1max=x1[1], ix1_bc=bc, ox1_bc=bc,
                nx2=1, x2min=-0.5, x2max=0.5, ix2_bc='periodic', ox2_bc='periodic',
                nx3=1, x3min=-0.5, x3max=0.5, ix3_bc='periodic', ox3_bc='periodic')


HST = lambda dt: dict(file_type='hst', dt=dt)
LOG = dict(file_type='log', dcycle=1)
def BIN(var, dt):
    return dict(file_type='bin', variable=var, dt=dt, bin_double='true')

tasks = []


def add(test, case, binary, text, meta, tmax_s=3600):
    d = f'{B}/runs/{test}/{case}'
    os.makedirs(d, exist_ok=True)
    open(f'{d}/run.athinput', 'w').write(text)
    meta.update(test=test, case=case, binary=binary)
    json.dump(meta, open(f'{d}/meta.json', 'w'), indent=1)
    tasks.append(f'{d} {binary} {tmax_s}')


def cpaw(dim):
    betas = [0.2, 1e-2, 1e-3, 1e-4, 1e-5]
    Ns = [32, 64, 128, 256]
    schemes = list(SCHEMES) if dim == 2 else S1D
    for beta in betas:
        pres = beta*0.5*(1.0 + 0.01)       # beta = p/(|B|^2/2), |B|^2 = 1 + 0.1^2
        for rec in ('plm', 'ppm4'):
            for sc in schemes:
                for N in Ns:
                    if dim == 2:
                        mesh = mesh2d(2*N, N, (0, 2), (0, 1))
                        T = 2/math.sqrt(5)
                    else:
                        mesh = mesh1d(N, (0, 1))
                        T = 1.0
                    mhd = dict(reconstruct=rec, gamma='1.6666666666666667', **SCHEMES[sc])
                    prob = dict(iprob=3, b_par=1.0, b_perp=0.1,
                                pres=repr(pres), right_polar='true')
                    outs = [HST(T/50), LOG, BIN('mhd_w_bcc', T),
                            BIN('mhd_w_e', T/40)]
                    case = f'b{beta:.0e}_{rec}_{sc}_N{N}'
                    add(f'cpaw{dim}d', case, LBT,
                        deck(mesh, dict(cfl_number=0.4, tlim=repr(T)), mhd, prob, outs),
                        dict(beta=beta, rec=rec, scheme=sc, N=N, pres=pres, T=T, dim=dim))


def lw():
    # background of the AthenaK deck (bx0, by0, bz0) = (1, sqrt 2, 0.5), |B|^2 = 3.25
    betas = [0.369, 1e-2, 1e-3, 1e-4]
    Ns = [32, 64, 128, 256]
    modes = {'fast': 0, 'alf': 1, 'slow': 2}
    for beta in betas:
        pg = 0.6 if beta == 0.369 else beta*0.5*3.25
        for mn, wf in modes.items():
            for sc in SCHEMES:
                for N in Ns:
                    mesh = mesh2d(2*N, N, (0, 3), (0, 1.5))
                    mhd = dict(reconstruct='plm', gamma='1.6666666666666667', **SCHEMES[sc])
                    prob = dict(pgen_name='linear_wave', wave_flag=wf, amp='1.0e-6', dens=1.0,
                                pgas=repr(pg), vx0=0.0, vy0=0.0, vz0=0.0, bx0=1.0,
                                by0='1.4142135623730951', bz0=0.5)
                    outs = [LOG]
                    case = f'b{beta:.0e}_{mn}_{sc}_N{N}'
                    add('lw', case, BIP,
                        deck(mesh, dict(cfl_number=0.4, tlim=1.0), mhd, prob, outs),
                        dict(beta=beta, mode=mn, scheme=sc, N=N, pgas=pg), tmax_s=7200)


def vortex():
    mu = 5.389489439
    kap = math.sqrt(2.0)*mu
    vtil = kap/(2*math.pi)
    tend = 0.05
    ncross = tend*vtil/(20.0*math.sqrt(2.0))
    for sc in SCHEMES:
        for N in [32, 64, 128, 256, 512, 1024]:
            mesh = mesh2d(N, N, (-10, 10), (-10, 10))
            mhd = dict(reconstruct='plm', gamma='1.6666666666666667', **SCHEMES[sc])
            prob = dict(pgen_name='balsara_vortex', vtil=repr(vtil), beta_k=0.5,
                        ncross=repr(ncross))
            outs = [HST(tend/100), LOG, BIN('mhd_w_bcc', tend),
                    BIN('mhd_w_e', tend/50)]
            add('vortex', f'{sc}_N{N}', BIP,
                deck(mesh, dict(cfl_number=0.4, tlim=tend), mhd, prob, outs),
                dict(scheme=sc, N=N, vtil=vtil, btil=math.sqrt(0.5)*vtil, tend=tend))


def rd():
    for sc in S1D:
        for N in [200, 800, 3200, 25600]:
            mesh = mesh1d(N, (0, 1), bc='outflow')
            mhd = dict(reconstruct='plm', gamma=2.0, **SCHEMES[sc])
            prob = dict(pgen_name='shock_tube', shock_dir=1, xshock=0.5,
                        dl=0.1, ul=0, vl=0, wl=0, pl=0.1, bxl=5.0, byl=0.1, bzl=0,
                        dr=0.1, ur=0, vr=0, wr=0, pr=0.1, bxr=5.0, byr=-0.1, bzr=0)
            outs = [HST(0.0002), LOG, BIN('mhd_w_bcc', 0.01)]
            add('rd', f'{sc}_N{N}', BIP,
                deck(mesh, dict(cfl_number=0.4, tlim=0.01), mhd, prob, outs),
                dict(scheme=sc, N=N))


def rotor():
    arms = dict(SCHEMES)
    arms['hlld_gs_fofc'] = dict(rsolver='hlld', fofc='true')
    arms['hlle_gs_fofc'] = dict(rsolver='hlle', fofc='true')
    for sc, keys in arms.items():
        for N in [200, 400, 800]:
            mesh = mesh2d(N, N, (0, 1), (0, 1), bc='outflow')
            mhd = dict(reconstruct='plm', gamma='1.6666666666666667', **keys)
            prob = dict(iprob=1, user_hist='true')
            outs = [HST(0.0006), LOG, BIN('mhd_w_bcc', 0.12)]
            add('rotor', f'{sc}_N{N}', LBT,
                deck(mesh, dict(cfl_number=0.4, tlim=0.12), mhd, prob, outs),
                dict(scheme=sc, N=N, tend=0.12))


def blast():
    arms = dict(SCHEMES)
    arms['hlld_gs_fofc'] = dict(rsolver='hlld', fofc='true')
    arms['hlle_gs_fofc'] = dict(rsolver='hlle', fofc='true')
    s2 = math.sqrt(2.0)
    setups = {
        # Tremblin 2024 8.2.3 (= GS05 blast): beta = 0.1/(2 pi) = 0.016
        'b2e-2': dict(dom=(0, 1), bc='periodic', c=0.5, gamma='1.6666666666666667',
                      p_in=10.0, p_amb=0.1, b0=2*math.sqrt(math.pi), tend=0.2),
        # Balsara & Spicer 1999 / MD21 6.4: beta = 2.5e-4
        'b2.5e-4': dict(dom=(-0.5, 0.5), bc='outflow', c=0.0, gamma=1.4,
                        p_in=1000.0, p_amb=0.1, b0=100/math.sqrt(4*math.pi), tend=0.01),
        # Tremblin 2024 8.2.4: B = 250/sqrt2 (1,1), beta = 3.2e-6
        'b3e-6': dict(dom=(0, 1), bc='periodic', c=0.5, gamma=1.4,
                      p_in=1000.0, p_amb=0.1, b0=250.0, tend=0.02),
    }
    for lab, st in setups.items():
        dims = [(2, 256)] + ([(3, 64), (3, 128)] if lab == 'b2.5e-4' else [])
        for dim, N in dims:
            for sc, keys in arms.items():
                if dim == 3 and sc not in ('hlle_gs', 'hlld_gs', 'hlld_sm', 'hlld_uct',
                                           'lhlld_gs', 'lhlld_uct', 'hlld_gs_fofc'):
                    continue
                d0, d1 = st['dom']
                if dim == 2:
                    mesh = mesh2d(N, N, (d0, d1), (d0, d1), bc=st['bc'])
                else:
                    mesh = mesh2d(N, N, (d0, d1), (d0, d1), bc=st['bc'], n3=N, x3=(d0, d1))
                mhd = dict(reconstruct='plm', gamma=st['gamma'], **keys)
                prob = dict(iprob=2, user_hist='true', p_in=st['p_in'], p_amb=st['p_amb'],
                            b0=repr(st['b0']), theta=repr(math.pi/2), phi=repr(math.pi/4),
                            xc=st['c'], yc=st['c'], zc=st['c'], r_blast=0.1)
                te = st['tend']
                outs = [HST(te/200), LOG, BIN('mhd_w_bcc', te)]
                add('blast', f'{lab}_{dim}d_N{N}_{sc}', LBT,
                    deck(mesh, dict(cfl_number=0.4, tlim=te), mhd, prob, outs),
                    dict(scheme=sc, N=N, dim=dim, setup=lab, tend=te,
                         beta=st['p_amb']/(0.5*st['b0']**2)), tmax_s=5400)


if __name__ == '__main__':
    t = sys.argv[1]
    {'cpaw2d': lambda: cpaw(2), 'cpaw1d': lambda: cpaw(1), 'lw': lw, 'vortex': vortex,
     'rd': rd, 'rotor': rotor, 'blast': blast}[t]()
    open(f'{B}/tasks_{t}.txt', 'w').write('\n'.join(tasks) + '\n')
    print(t, len(tasks), 'tasks')
