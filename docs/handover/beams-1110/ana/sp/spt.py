#!/usr/bin/env python3
"""sptests-1011 driver: spt.py run <case> <arm> [levels...] | spt.py dir <case> <arm> <lev>

Cases (NOTE-sptests.md 0.4).  A level is nx1 (radial cases: 4 x 4 lateral cells; 3-D cases:
the lateral grid of the case), or for a T ladder (case suffix _T) the dt divisor l (dt =
dt0/l on the fixed grid NT).  Suffix _s: poly-stretched r.  Arms: spl.arm_keys().
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spl  # noqa: E402

TH0, TH1, PH1 = spl.TH0, spl.TH1, spl.PH1
CG = 1.0e64        # arad of the held-gas cases: B = a T^4 with a tiny T (no hydro dt limit)


def rad_base(c=1.0, kp=0.0, ks=0.0, cfl=1.0e3):
    return {'c_light': c, 'chat_over_c': 1.0, 'e_floor': '1.0e-30', 'opacity': 'const',
            'kappa_p': kp, 'kappa_e': kp, 'kappa_f': kp, 'kappa_s': ks, 'arad': 1.0,
            'coupling': 'true', 'gas_feedback': 'false', 'transport': 'implicit',
            'implicit_cfl': cfl, 'implicit_solver': 'bicgstab', 'implicit_precond': 'line',
            'implicit_maxit': 200, 'implicit_tol': '1.0e-10', 'implicit_lin_tol': '1.0e-12',
            'implicit_closure_lag': 'step', 'marshak_q': 0.5,
            # prodcand e66ee669 runs ic_rad = fs_relax on fresh vet_gd starts by default: the
            # tests set their own initial state (time-dependent cases need it exactly)
            'ic_rad': 'none'}


def shell(**kw):
    p = {'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_shell', 'gas_rho': 1.0,
         'gas_temp': '1.0e-16', 'e_out': '1.0e-8', 'e_amp': 0.0, 'r0': 2.0, 'width': 0.2}
    p.update(kw)
    return p


def case_cfg(case):
    """geometry + physics of a case"""
    c = {'poly': False, 'geo': 'r', 'nt': None, 'gas': False, 'hwall': 'outflow',
         'nlat': (4, 4), 'th': (TH0, TH1), 'ph': (0.0, PH1), 'th_bc': 'reflect',
         'outs': 'tab', 'nout': 1, 'np': 1, 'nthr': 1}
    name = case
    if name.endswith('_T'):
        c['tref'] = True
        name = name[:-2]
    if name.endswith('_s'):
        c['poly'] = True
        name = name[:-2]
    elif name.endswith('_u'):
        name = name[:-2]
    cart = name.endswith('_cart')
    if cart:
        name = name[:-5]
        c['geo'] = 'cart'
    tok = name.split('_')
    t = tok[0]
    rk = None
    nlx = None
    for x in tok[1:]:
        if x.startswith('r') and x[1:].isdigit():
            rk = int(x[1:])
        if x.startswith('n') and x[1:].isdigit():
            nlx = int(x[1:])
    c['test'] = t
    if t == 'thin':
        # T6a: black-body inner sphere (flux 1) in a thin medium, steady
        c.update(r=(1.0, 10.0), steady=True, nlim=40,
                 rad=dict(rad_base(1.0, ks=1.0e-6, cfl=1.0e3), implicit_bc_x1min='flux',
                          implicit_flux_x1min=1.0, implicit_bc_x1max='marshak',
                          implicit_ebath_x1max=0.0, vet_gd_inner='emit'),
                 prob=shell())
    elif t == 'core':
        # T6b: IDORT 4.1.2 homogeneous sphere rho kappa_a = 10, B = 1 inside r = 1
        tc = (1.0/CG)**0.25
        c.update(kscale=1.0e-30, r=(0.2, 7.0), steady=True, nlim=60,
                 rad=dict(rad_base(1.0, kp=1.0, cfl=1.0e3), arad=CG,
                          implicit_bc_x1min='reflect', implicit_bc_x1max='marshak',
                          implicit_ebath_x1max=0.0),
                 prob=shell(gas_rho='1.0e-7', gas_temp='1.0e-30', core_r=1.0,
                            core_rho=10.0, core_temp='%.17g' % tc, gas_hold='true',
                            user_srcs='true'))
    elif t == 'satm':
        # T10: spherical extended grey atmosphere, rho = 100 r^-2, pure scattering,
        # flux 1 in, Marshak top, steady (sph_atm, Eddington start)
        c.update(r=(1.0, 5.0), steady=True, nlim=60,
                 rad=dict(rad_base(100.0, ks=100.0, cfl=1.0e4), implicit_bc_x1min='flux',
                          implicit_flux_x1min=1.0, implicit_bc_x1max='marshak',
                          implicit_ebath_x1max=0.0),
                 prob={'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_atm', 'atm_rho0': 1.0,
                       'atm_rho_n': 2.0, 'atm_temp': '1.0e-16', 'atm_init': 'eddington',
                       'e_out': '1.0e-6', 'user_srcs': 'true'})
    elif t in ('patm', 'patmT'):
        # PRIORITY (photosphere): rho kappa = 896 r^-8 on r in [1, 5]: tau = 128 at the
        # base, tau = 1 at r = 2 (R_ph), thin envelope to 2.5 R_ph; flux 1 in, vacuum top;
        # patm pure scattering (S = J), patmT absorption with the gas in radiative
        # equilibrium (gas_feedback, the density held); Eddington start; steady at a
        # production-like c dt/dx_min = 1e3
        rad = dict(rad_base(1.0, ks=896.0, cfl=1.0e3), implicit_bc_x1min='flux',
                   implicit_flux_x1min=1.0, implicit_bc_x1max='marshak',
                   implicit_ebath_x1max=0.0)
        if t == 'patmT':
            rad.update(kappa_p=896.0, kappa_e=896.0, kappa_f=896.0, kappa_s=0.0,
                       arad=CG, gas_feedback='true')
            # gas heat capacity ~ the radiation energy at the top (with a = 1e64, T ~
            # 1e-16): rho x 1e20, kappa x 1e-20 (unscaled, rho c_v T / E ~ 1e-20 at the top
            # and the exchange Newton returns T = NaN there at c dt/dx >= 1e3)
            c['kscale'] = 1.0e-20
            for x in tok[1:]:
                if x.startswith('e') and x[1:].isdigit():
                    # coupled-gas scattering split: kappa_a = eps chi, kappa_s = (1 - eps) chi
                    # (steady radiative equilibrium S = J: the same exact reference)
                    eps = 10.0**(-int(x[1:]))
                    rad.update(kappa_p=896.0*eps, kappa_e=896.0*eps, kappa_f=896.0*eps,
                               kappa_s=896.0*(1.0 - eps))
        c.update(r=(1.0, 5.0), steady=True, nlim=nlx or 150, two=(80, nlx or 300, '1000.0'),
                 rad=rad, gas=(t == 'patmT'),
                 prob={'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_atm', 'atm_rho0': 1.0,
                       'atm_rho_n': 8.0, 'atm_temp': '1.0e-30', 'atm_init': 'eddington',
                       'e_out': '1.0e-6', 'user_srcs': 'true'})
    elif t == 'pinc':
        # PRIORITY (thermal beam, collimated): hot opaque core r < 1.5 (B = 1) under a cold
        # opaque shell 2 < r < 4 pierced by a radial cone channel of half-angle 0.06 on the
        # wedge axis (the vet_gd analogue of the Cartesian channel beams xb20/ba0)
        tc = (1.0/CG)**0.25
        cfl = 300.0
        for x in tok[1:]:
            if x.startswith('c') and x[1:].isdigit():
                cfl = float(x[1:])
        c.update(kscale=1.0e-30, geo='3d', r=(1.0, 5.0), steady=False, T=40.0,
                 dt0=None, cfl=cfl, nlx=nlx,
                 th=(0.5*math.pi - 0.35, 0.5*math.pi + 0.35),
                 ph=(0.0, math.pi/3.0), nlat=(48, 48), th_bc='periodic', outs='bin',
                 nthr=8, nout=4,
                 rad=dict(rad_base(1.0, kp=1.0), arad=CG, implicit_bc_x1min='reflect',
                          implicit_bc_x1max='marshak', implicit_ebath_x1max=0.0),
                 prob=shell(gas_rho='1.0e-6', gas_temp='1.0e-70', core_r=1.5,
                            core_rho=50.0, core_temp='%.17g' % tc, screen_r0=2.0,
                            screen_r1=4.0, screen_rho=500.0, screen_temp='1.0e-70',
                            hole_th='%.17g' % (0.5*math.pi), hole_ph='%.17g' %
                            (math.pi/6.0), hole_a=0.06, gas_hold='true',
                            user_srcs='true'))
    elif t == 'pin':
        # PRIORITY (thermal beam on vet_gd): a hot opaque core (r < 1.5, rho kappa_a 50,
        # B = 1) behind a cold opaque screen 1.6 < r < 1.8 (rho kappa_a 500) with a hole
        # of angular radius 0.1 on the wedge axis: a pinhole thermal beam; thin elsewhere
        tc = (1.0/CG)**0.25
        c.update(kscale=1.0e-30, geo='3d', r=(1.0, 5.0), steady=False, T=40.0,
                 dt0=None, cfl=300.0, nlx=nlx,
                 th=(0.5*math.pi - 0.35, 0.5*math.pi + 0.35),
                 ph=(0.0, math.pi/3.0), nlat=(48, 48), th_bc='periodic', outs='bin',
                 nthr=8, nout=4,
                 rad=dict(rad_base(1.0, kp=1.0), arad=CG, implicit_bc_x1min='reflect',
                          implicit_bc_x1max='marshak', implicit_ebath_x1max=0.0),
                 prob=shell(gas_rho='1.0e-6', gas_temp='1.0e-70', core_r=1.5,
                            core_rho=50.0, core_temp='%.17g' % tc, screen_r0=2.0,
                            screen_r1=2.2, screen_rho=500.0, screen_temp='1.0e-70',
                            hole_th='%.17g' % (0.5*math.pi), hole_ph='%.17g' %
                            (math.pi/6.0), hole_a=0.12, gas_hold='true',
                            user_srcs='true'))
    elif t.startswith('nlte'):
        # T2: IDORT 4.1.3 non-LTE scattering atmosphere, eps = 10^-le, on a shell of
        # width 20 at r0 = 20 rk (rk = 1 -> r in [10, 30]); bath 1 below (thick), vacuum top
        le = int(t[4:])
        eps = 10.0**(-le)
        r0 = 20.0*(rk or 1)
        tg = (1.0/CG)**0.25
        c.update(kscale=1.0e-30, r=(r0 - 10.0, r0 + 10.0), steady=True, nlim=nlx or 80,
                 rad=dict(rad_base(1.0, kp=eps, ks=1.0 - eps, cfl=1.0e6), arad=CG,
                          implicit_bc_x1min='marshak', implicit_ebath_x1min=1.0,
                          implicit_bc_x1max='marshak', implicit_ebath_x1max=0.0),
                 prob={'pgen_name': 'rad_m1_beam', 'm1_test': 'sph_atm',
                       'atm_rho0': '%.17g' % (1.0e-3*math.exp(20.0)), 'atm_scale_h': 1.0,
                       'atm_temp': '%.17g' % tg, 'atm_init': 'vacuum', 'e_out': '1.0e-6',
                       'user_srcs': 'true'})
        if cart:
            c['r'] = (10.0, 30.0)
        if 'tb' in tok[1:]:
            # two-stage: A 80 steps at implicit_cfl 1e6, B 300 steps at c dt/dx 1e3
            c['two'] = (nlx or 80, 300, '1000.0')
    elif t == 'pulse':
        # T3: Gaussian pulse rho kappa_s = 1 between vacuum ends, width 2 at r0 = 2 rk
        r0 = 2.0*(rk or 1)
        c.update(r=(r0 - 1.0, r0 + 1.0), steady=False, T=1.0, dt0=1.0/64,
                 rad=dict(rad_base(1.0, ks=1.0), implicit_bc_x1min='marshak',
                          implicit_ebath_x1min=0.0, implicit_bc_x1max='marshak',
                          implicit_ebath_x1max=0.0, vet_gd_inner='emit',
                          implicit_tol='1.0e-12', implicit_lin_tol='1.0e-13'),
                 prob=shell(e_out='1.0e-8', e_amp=1.0, r0=r0, width=0.2))
        if cart:
            c['r'] = (1.0, 3.0)
            c['prob']['r0'] = 2.0
    elif t == 'dif':
        # T4: J21 static diffusion: E = exp(-40 (r - r0)^2), rho kappa_s = 4e4, c = 10;
        # r0 = 3 rk on [r0 - 2, r0 + 2]; outputs at sqrt(Dt) = 0.13, 0.18, 0.22
        r0 = 3.0*(rk or 1)
        D = 10.0/(3.0*4.0e4)
        c.update(r=(r0 - 2.0, r0 + 2.0), steady=False, T=0.22**2/D, dt0=0.22**2/D/512,
                 outt=[0.13**2/D, 0.18**2/D, 0.22**2/D],
                 rad=dict(rad_base(10.0, ks=4.0e4), implicit_bc_x1min='reflect',
                          implicit_bc_x1max='marshak',
                          implicit_ebath_x1max='%.17g' % math.exp(-10.0),
                          implicit_tol='1.0e-12', implicit_lin_tol='1.0e-13'),
                 prob=shell(e_out='%.17g' % math.exp(-10.0), e_amp=1.0, r0=r0,
                            width='%.17g' % (1.0/math.sqrt(40.0))))
        if cart:
            c['r'] = (1.0, 5.0)
            c['prob']['r0'] = 3.0
    elif t.startswith('relax'):
        # T1: thermal relaxation (JSD12 sec 4 / IDORT 4.3.1): rho kappa_a = 100, c = 100,
        # a = 1; relax1 (T, E) = (1, 100), relax2 (100, 1); reflect in, Marshak bath E_eq
        T0, E0 = (1.0, 100.0) if t == 'relax1' else (100.0, 1.0)
        eq = relax_eq(T0, E0)
        c.update(r=(1.0, 2.0), steady=False, T=1.0e-3, dt0=1.0e-3/8,
                 outt=[2.5e-4, 5.0e-4, 1.0e-3], gas=True, hwall='reflect',
                 rad=dict(rad_base(100.0, kp=100.0), gas_feedback='true',
                          implicit_bc_x1min='reflect', implicit_bc_x1max='marshak',
                          implicit_ebath_x1max='%.17g' % eq[1],
                          implicit_gas_newton='true',
                          implicit_tol='1.0e-12', implicit_lin_tol='1.0e-13'),
                 prob=shell(gas_rho=1.0, gas_temp=T0, e_out=E0, e_amp=0.0))
        if 'early' in tok[1:]:
            # the transient itself (relaxation time 1/(c kappa) = 1e-4): outputs at 0.5, 1,
            # 2 relaxation times, dt0 = 2.5e-5 (the T ladder resolves it from dt0/1)
            c.update(T=2.0e-4, dt0=2.5e-5, outt=[5.0e-5, 1.0e-4, 2.0e-4])
    elif t in ('rsw', 'rswt', 'rswv', 'rswtv'):
        # T9: vcpred rswm: radiation-dominated moving shell (P_rad/P_gas ~ 10), the
        # radiation-acoustic pulse, rho kappa_s = 1000 (rswt: 10), kappa_p = 1, t = 2;
        # sptestsC: rswv / rswtv add a radial flow v_r = 0.01 sin(pi (r - 1)/2) (moving shell)
        ks = 10.0 if t.startswith('rswt') else 1000.0
        c.update(r=(1.0, 3.0), steady=False, T=2.0, dt0=1.0/64, gas=True,
                 hwall='reflect',
                 rad=dict(rad_base(1.0, kp=1.0, ks=ks), gas_feedback='true',
                          arad='3.0e10', implicit_bc_x1min='reflect',
                          implicit_bc_x1max='marshak', implicit_ebath_x1max='0.03',
                          implicit_tol='1.0e-12', implicit_lin_tol='1.0e-13'),
                 prob=shell(e_out='0.03', e_amp='0.003', gas_teq='true'))
        c['rad']['kappa_s'] = ks
        if t.endswith('v'):
            c['prob']['gas_vr'] = '0.01'
        if 'q' in tok[1:]:
            # quiet: no pulse (the uniform equilibrium must stay put)
            c['prob']['e_amp'] = 0.0
        if rk:
            # sptestsC large-radius control: r in [2 rk - 1, 2 rk + 1]
            c['r'] = (2.0*rk - 1.0, 2.0*rk + 1.0)
            c['prob']['r0'] = 2.0*rk
        if cart:
            c['r'] = (1.0, 3.0)
            c['prob']['r0'] = 2.0
    elif t in ('rlat', 'rlatt'):
        # sptestsC lateral radiation-acoustic wave (L_z conservation): the rsw medium (c 1,
        # arad 3e10, E 0.03, gas at radiative equilibrium, P_rad/P_gas ~ 10), rho kappa_s
        # 1000 (rlatt: 10), rho kappa_p 1, E x (1 + 0.1 cos(2 pi phi/(pi/3))), at rest;
        # r in [1, 2], 16 radial x 4 theta x lev phi cells, t = 4; phi periodic, theta and
        # radial walls reflect: total L_z = 0 is conserved exactly
        ks = 10.0 if t == 'rlatt' else 1000.0
        c.update(r=(1.0, 2.0), steady=False, T=4.0, dt0=4.0/256, xdt='scalephi',
                 outt=[0.5], gas=True, hwall='reflect', geo='3d', nlat=(4, 64),
                 outs='bin', nthr=8,
                 rad=dict(rad_base(1.0, kp=1.0, ks=ks), gas_feedback='true',
                          arad='3.0e10', implicit_bc_x1min='reflect',
                          implicit_bc_x1max='marshak', implicit_ebath_x1max='0.03',
                          implicit_tol='1.0e-12', implicit_lin_tol='1.0e-13'),
                 prob=shell(e_out='0.03', e_amp=0.0, gas_teq='true', lat_amp=0.1,
                            lat_l2=0.0, lat_l3=1.0))
        for x in tok[1:]:
            if x == 'rt':
                # closed box: reflecting radiation top (vet_col_reflect_top mirrors the FS)
                c['rad']['implicit_bc_x1max'] = 'reflect'
                c['rad'].pop('implicit_ebath_x1max')
            if x.startswith('om') and x[2:].isdigit():
                # solid-body rotation Omega = k 1e-3 (net L_z; radiation dragged along)
                c['prob']['gas_rot'] = '%.17g' % (int(x[2:])*1.0e-3)
    elif t.startswith('marsh'):
        # T5 (sptestsC): spherical Marshak wave: cold gas (rho 1, T 1e-3, gamma 5/3), rho
        # kappa_a = 10, c = 100, a = 1; Marshak bath E = 1 (T = 1) at r_min, vacuum-ish
        # top; width 1 at r_min = rk (r1: [1, 2], r1000: [1000, 1001]); t = 0.05
        r0 = 1.0*(rk or 1)
        c.update(r=(r0, r0 + 1.0), steady=False, T=0.05, dt0=0.05/256,
                 outt=[0.0125], gas=True, hwall='reflect',
                 rad=dict(rad_base(100.0, kp=10.0), gas_feedback='true',
                          implicit_bc_x1min='marshak', implicit_ebath_x1min=1.0,
                          implicit_bc_x1max='marshak', implicit_ebath_x1max='1.0e-12',
                          implicit_gas_newton='true', vet_gd_inner='emit',
                          implicit_tol='1.0e-12', implicit_lin_tol='1.0e-13'),
                 prob=shell(gas_rho=1.0, gas_temp='1.0e-3', e_out='1.0e-12', e_amp=0.0))
        if cart:
            c['r'] = (1.0, 2.0)
    elif t == 'drot':
        # dynamic diffusion on the wedge (sptestsC): heavy gas (rho 1e6, kappa x 1e-6) in
        # solid-body rotation Omega = 0.03 about z, rho kappa_s = 1000, c = 10; E = 1 +
        # 0.1 cos(2 pi phi/(pi/3)) advected by the gas (v tau/c ~ 7) and diffused
        # (D k^2 t ~ 0.5); r in [1, 2], 16 radial x 4 theta x lev phi cells; t = 10
        c.update(kscale=1.0e-6, r=(1.0, 2.0), steady=False, T=10.0, dt0=10.0/400,
                 xdt='scalephi', outt=[1.0], gas=True, hwall='reflect', geo='3d',
                 nlat=(4, 64), outs='bin', nthr=8,
                 rad=dict(rad_base(10.0, ks=1000.0), gas_feedback='true', arad=1.0,
                          implicit_bc_x1min='reflect', implicit_bc_x1max='marshak',
                          implicit_ebath_x1max=1.0,
                          implicit_tol='1.0e-12', implicit_lin_tol='1.0e-13'),
                 prob=shell(gas_rho=1.0, gas_temp='1.0e-6', e_out=1.0, e_amp=0.0,
                            lat_amp=0.1, lat_l2=0.0, lat_l3=1.0, gas_rot=0.03))
    elif t in ('shd', 'beam', 'beam1'):
        # T7 / T8: 3-D wedge (theta about the equator, phi a 30 deg period)
        c.update(kscale=1.0e-30, geo='3d', r=(1.0, 6.0 if t == 'shd' else 5.0), steady=True, nlim=40,
                 th=(0.5*math.pi - 0.35, 0.5*math.pi + 0.35),
                 ph=(0.0, math.pi/3.0), nlat=(48, 48), th_bc='periodic',
                 outs='bin', nthr=16,
                 rad=dict(rad_base(1.0, ks=1.0e-6, cfl=1.0e3), implicit_bc_x1min='flux',
                          implicit_flux_x1min=1.0, implicit_bc_x1max='marshak',
                          implicit_ebath_x1max=0.0, vet_gd_inner='emit'),
                 prob=shell())
        for x in tok[1:]:
            if x.startswith('c') and x[1:].isdigit():
                # c dt/dx_min of the steady run (default 1e3)
                c['rad']['implicit_cfl'] = float(x[1:])
        if t == 'shd':
            c['rad'].update(kappa_p=1.0, kappa_e=1.0, kappa_f=1.0, kappa_s=0.0)
            c['prob'].update(gas_rho='1.0e-6', gas_temp='1.0e-30', clump_r=3.0,
                             clump_th='%.17g' % (0.5*math.pi), clump_ph='%.17g' %
                             (math.pi/6.0), clump_a=0.4, clump_rho=1.0e3,
                             clump_temp='1.0e-30', gas_hold='true', user_srcs='true')
        else:
            # two beams crossing at r ~ 2.3 (beam: both; beam1: the lower one only)
            c['beams'] = [((0.5*math.pi + 0.15, 0.5*math.pi + 0.25), +1)]
            if t == 'beam':
                c['beams'].append(((0.5*math.pi - 0.25, 0.5*math.pi - 0.15), -1))
    else:
        raise SystemExit('unknown case ' + case)
    for x in tok[1:]:
        if x.startswith('h') and x[1:].isdigit() and not c['steady']:
            # sptestsC: k equally spaced outputs (conservation histories)
            c['outt'] = [c['T']/int(x[1:])]
    if 'one' in tok[1:]:
        c.pop('two', None)
    for x in tok[1:]:
        if x.startswith('A') and x[1:].isdigit() and c.get('two'):
            c['acfl'] = '1.0e%s' % x[1:]
    # lateral-grid tokens: L<n> = n x n cells over the default wedge; Lp = production-like
    # angular cells (0.008 rad): 8 theta cells over 0.064 rad (periodic), 64 phi over pi/6
    for x in tok[1:]:
        if x == 'Lq':
            # production-like angular cells (0.008 rad) and a band wide enough for the
            # near-tangent reach (64 x 64 lateral cells: no clamped reads)
            c.update(nlat=(64, 64), th=(0.5*math.pi - 0.256, 0.5*math.pi + 0.256),
                     ph=(0.0, math.pi/6.0), th_bc='periodic')
            c['nthr'] = 16
        if x == 'Lp':
            c.update(nlat=(8, 64), th=(0.5*math.pi - 0.032, 0.5*math.pi + 0.032),
                     ph=(0.0, math.pi/6.0), th_bc='periodic')
            c['nthr'] = max(c['nthr'], 4)
        elif x.startswith('L') and x[1:].isdigit():
            c['nlat'] = (int(x[1:]), int(x[1:]))
            c['nthr'] = max(c['nthr'], 4)
    return c


def relax_eq(T0, E0, rho=1.0, gm1=2.0/3.0, a=1.0):
    etot = rho*T0/gm1 + E0
    lo, hi = 0.0, (etot/a)**0.25 + etot
    for _ in range(200):
        m = 0.5*(lo + hi)
        if rho*m/gm1 + a*m**4 > etot:
            hi = m
        else:
            lo = m
    return lo, a*lo**4


def dxmin(c, n):
    rf = spl.rfaces(c['r'][0], c['r'][1], n, c['poly'])
    dr = np.diff(rf).min()
    if c['geo'] == 'cart':
        return dr
    n2, n3 = c['nlat']
    dth = (c['th'][1] - c['th'][0])/n2
    dph = (c['ph'][1] - c['ph'][0])/n3
    smin = min(math.sin(c['th'][0]), math.sin(c['th'][1]), 1.0)
    return min(dr, rf[0]*dth, rf[0]*smin*dph)


def mk(case, arm, lev):
    c = case_cfg(case)
    d = os.path.join(spl.RUNS, case, '%s_%d' % (arm, lev))
    os.makedirs(d, exist_ok=True)
    tref = c.get('tref', False)
    n1 = 128 if tref else lev
    if tref and c['test'] == 'relax':
        n1 = 64
    n2, n3 = c['nlat']
    if c['test'] in ('drot', 'rlat', 'rlatt'):
        n1 = 16
        if not tref:
            n3 = lev
        c['nlat'] = (n2, n3)
    sp = c['geo'] != 'cart'
    mesh = {'nghost': 2, 'nx1': n1, 'x1min': c['r'][0], 'x1max': c['r'][1],
            'ix1_bc': c['hwall'], 'ox1_bc': c['hwall']}
    if sp:
        mesh.update({'nx2': n2, 'x2min': '%.17g' % c['th'][0],
                     'x2max': '%.17g' % c['th'][1], 'ix2_bc': c['th_bc'],
                     'ox2_bc': c['th_bc'], 'nx3': n3, 'x3min': '%.17g' % c['ph'][0],
                     'x3max': '%.17g' % c['ph'][1], 'ix3_bc': 'periodic',
                     'ox3_bc': 'periodic', 'use_spherical_polar': 'true'})
        mb = {'nx1': n1, 'nx2': n2, 'nx3': n3}
    else:
        # vet_sc needs a multi-D mesh: 4 periodic x2 cells
        mesh.update({'nx2': 4, 'x2min': -0.5, 'x2max': 0.5, 'ix2_bc': 'periodic',
                     'ox2_bc': 'periodic', 'nx3': 1, 'x3min': -0.5, 'x3max': 0.5,
                     'ix3_bc': 'periodic', 'ox3_bc': 'periodic'})
        mb = {'nx1': n1, 'nx2': 4, 'nx3': 1}
    if c['poly']:
        mesh['use_grid_stretch_r_poly'] = 'true'
        mesh.update(spl.POLY)
    rad = dict(c['rad'])
    prob = dict(c['prob'])
    ksc = c.get('kscale')
    if ksc:
        # held-gas emitters: the opacity per mass x K, every density / K (rho kappa and
        # B = a T^4 unchanged): the gas heat capacity per step >> the exchange, so a held
        # temperature really emits, and T ~ 1e-16 keeps the hydro dt out of the way
        for k in ('kappa_p', 'kappa_e', 'kappa_f', 'kappa_s'):
            rad[k] = '%.17g' % (float(rad[k])*ksc)
        for k in ('gas_rho', 'core_rho', 'clump_rho', 'screen_rho', 'atm_rho0'):
            if k in prob:
                prob[k] = '%.17g' % (float(prob[k])/ksc)
    rad.update(spl.arm_keys(arm))
    if not sp:
        rad.pop('vet_gd_inner', None)
    if 'vet_gd_inner' in spl.arm_keys(arm):
        rad['vet_gd_inner'] = spl.arm_keys(arm)['vet_gd_inner']
    if c['geo'] == '3d' and 'beams' in c:
        # vet_gd_beam (one beam per input key set: the lower beam; the second via a
        # second run is not possible -> the crossing case uses a symmetric pair)
        pass
    if str(rad.get('gas_feedback')) == 'false':
        rad.pop('implicit_vimp', None)
        rad.pop('time2_vstage', None)
    hyd = {'eos': 'ideal', 'reconstruct': 'plm', 'rsolver': 'hllc',
           'gamma': '1.666666666666667'}
    if not c['gas'] or c['test'] in ('patm',):
        # held / static gas: no radiation-pressure signal speed in the hydro dt
        hyd['rad_signal_speed'] = 'false'
    if c['test'] == 'patm' or case.startswith('patmT'):
        hyd['rad_signal_speed'] = 'false'
    if sp and c['th_bc'] == 'periodic':
        hyd['sp_x2_periodic_image'] = 'true'
    time = {'evolution': 'dynamic', 'integrator': 'rk2', 'cfl_number': 0.3,
            'nlim': 1000000, 'tlim': 1.0e30, 'ndiag': 1000000}
    outs = []
    if c['steady']:
        time['nlim'] = c['nlim']
        dtout = 1.0e30
    elif c.get('cfl'):
        rad['implicit_cfl'] = c['cfl']
        time['nlim'] = c.get('nlx') or 100
        # an output every 100 cycles
        dtout = 100.0*float(c['cfl'])*dxmin(c, n1)/float(rad['c_light'])*(1.0 - 1e-9)
    else:
        dt = c['dt0']/(lev if tref else 1.0)
        if not tref and c.get('xdt') == 'scale':
            dt = c['dt0']*(n1/128.0)
        if not tref and c.get('xdt') == 'scalephi':
            dt = c['dt0']*(64.0/n3)
        rad['implicit_cfl'] = '%.17g' % (dt*float(rad['c_light'])/dxmin(c, n1))
        time['tlim'] = '%.17g' % c['T']
        dtout = c['T']
    if rad.get('vet_gd_rotate_dt') == 'AUTO':
        # final_1012 FIN: physical-time rotation.  T ladders (tref): 10 x the coarsest dt
        # (dt0), the same rotation times on every level; else 20 x the run's radiation dt.
        # SPT_ROT=0 -> rotate 0 control (vet_gd_rotate_dt 0, vet_gd_rotate_every 0)
        if os.environ.get('SPT_ROT') == '0' or not sp or str(rad.get('vet_gd_rotate_every')) == '0':
            rad['vet_gd_rotate_dt'] = 0.0
            rad['vet_gd_rotate_every'] = 0
        else:
            if tref and not c['steady'] and not c.get('cfl'):
                rdt = 10.0*c['dt0']*(1.0 - 2.0e-6)
            else:
                rdt = 20.0*float(rad['implicit_cfl'])*dxmin(c, n1)/float(rad['c_light'])
            rad['vet_gd_rotate_dt'] = '%.17g' % rdt
    ot = c.get('outt')
    vars_ = ['m1'] + (['hydro_w', 'hydro_u'] if c['gas'] else [])
    if c['test'] in ('pinc', 'pin'):
        # P2 diagnosis: the vet_gd formal-solution moments and the transported face fluxes
        vars_ += ['m1_fs', 'm1_face']
    elif sp and os.environ.get('SPT_DIAG'):
        # the vet_gd formal-solution moments (J_fs, H_fs) next to the M1 state
        vars_ += ['m1_fs']
    if os.environ.get('SPT_GDHR') and sp:
        vars_ += ['m1_gdhr']
    for v in vars_:
        if c['outs'] == 'tab':
            o = {'file_type': 'tab', 'variable': v, 'data_format': '%24.16e'}
            if sp:
                o.update({'slice_x2': '%.17g' % (c['th'][0] + 0.35*(c['th'][1] -
                                                                     c['th'][0])),
                          'slice_x3': '%.17g' % (c['ph'][0] + 0.1*(c['ph'][1] -
                                                                   c['ph'][0]))})
            else:
                o['slice_x2'] = 0.1
        else:
            o = {'file_type': 'bin', 'variable': v}
        if ot and not c['steady']:
            # equal spacing that hits every requested time
            g = ot[0]
            o['dt'] = '%.17g' % g
        else:
            o['dt'] = '%.17g' % dtout
        outs.append(o)
    blocks = [('job', {'basename': 'sp'}), ('mesh', mesh), ('meshblock', mb),
              ('time', time), ('hydro', hyd), ('rad_m1', rad), ('problem', prob)]
    for q, o in enumerate(outs):
        blocks.append(('output%d' % (q + 1), o))
    if c['geo'] == '3d' and 'beams' in c:
        (t0, t1), sg = c['beams'][0]
        rad.update({'vet_gd_beam': 'true', 'vet_gd_rotate_every': 0,
                    'vet_gd_gl_stagger': 'false',
                    'vet_gd_beam_th': '%.17g' % (0.5*math.pi - sg*0.15),
                    'vet_gd_beam_ph': '%.17g' % (math.pi/6.0),
                    'vet_gd_beam_th0': '%.17g' % t0, 'vet_gd_beam_th1': '%.17g' % t1,
                    'vet_gd_beam_ph0': '%.17g' % (math.pi/6.0 - 0.06),
                    'vet_gd_beam_ph1': '%.17g' % (math.pi/6.0 + 0.06),
                    'vet_gd_beam_e': 1.0})
        if len(c['beams']) > 1:
            (u0, u1), sg2 = c['beams'][1]
            rad.update({'vet_gd_beam2_th': '%.17g' % (0.5*math.pi - sg2*0.15),
                        'vet_gd_beam2_th0': '%.17g' % u0,
                        'vet_gd_beam2_th1': '%.17g' % u1})
    if c.get('two'):
        # two-stage steady state: A at implicit_cfl 1e6 (the deep diffusion relaxes), a
        # restart file, then B (launch) restarts at the case's production-like cfl
        rad['implicit_cfl'] = c.get('acfl', '1.0e6')
        time['nlim'] = c['two'][0]
        rotb = os.path.join(d, 'rotB.txt')
        if os.path.exists(rotb):
            os.remove(rotb)
        if sp and float(rad.get('vet_gd_rotate_dt', 0) or 0) > 0:
            # final_1012 FIN: stage A rotates at 20 x its own dt, stage B at 20 x its dt
            dxm = dxmin(c, n1)/float(rad['c_light'])
            rad['vet_gd_rotate_dt'] = '%.17g' % (20.0*float(rad['implicit_cfl'])*dxm)
            open(rotb, 'w').write('%.17g\n' % (20.0*float(c['two'][2])*dxm))
        blocks.append(('output%d' % (len(outs) + 1), {'file_type': 'rst', 'dt': 1.0e30}))
    inp = os.path.join(d, 'in.athinput')
    spl.mkinput(inp, blocks)
    return d, inp, c


def launch(case, arm, levs):
    procs = []
    for lev in levs:
        d, inp, c = mk(case, arm, lev)
        for f in ('tab', 'bin', 'log', 'rst', 'logA'):
            p = os.path.join(d, f)
            if os.path.isdir(p):
                os.system('rm -rf ' + p)
        procs.append((lev, spl.run(spl.BIN['new'], inp, d, c['np'], c['nthr'],
                                   wait=False)))
    for lev, p in procs:
        print(case, arm, lev, 'rc', p.wait(), flush=True)
    c = case_cfg(case)
    if c.get('two'):
        procs = []
        for lev in levs:
            d = os.path.join(spl.RUNS, case, '%s_%d' % (arm, lev))
            rs = sorted(f for f in os.listdir(os.path.join(d, 'rst')) if f.endswith('.rst'))
            os.rename(os.path.join(d, 'log'), os.path.join(d, 'logA'))
            xa = ('time/nlim=%d' % (c['two'][0] + c['two'][1]),
                  'rad_m1/implicit_cfl=%s' % c['two'][2])
            if os.path.exists(os.path.join(d, 'rotB.txt')):
                xa = xa + ('rad_m1/vet_gd_rotate_dt=%s'
                           % open(os.path.join(d, 'rotB.txt')).read().strip(),)
            procs.append((lev, spl.run(spl.BIN['new'], None, d, c['np'], c['nthr'],
                                       args=xa,
                                       wait=False, rst=os.path.join(d, 'rst', rs[-1]))))
        for lev, p in procs:
            print(case, arm, lev, 'B rc', p.wait(), flush=True)


def main():
    act = sys.argv[1]
    if act == 'run':
        case, arm = sys.argv[2:4]
        levs = [int(x) for x in sys.argv[4:]]
        launch(case, arm, levs)
    elif act == 'mk':
        case, arm = sys.argv[2:4]
        for lev in [int(x) for x in sys.argv[4:]]:
            print(mk(case, arm, lev)[1])


if __name__ == '__main__':
    main()
