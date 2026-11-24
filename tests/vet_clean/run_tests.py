#!/usr/bin/env python3
"""run_tests.py: the clean VET radiation test suite (CPU; Python standard library + numpy).

usage: run_tests.py --exe PATH [--quick] [--jobs N] [--threads T] [--only NAME,...]
                    [--workdir DIR] [--mpirun CMD] [--keep] [--record FILE]

Every test writes its runs under WORKDIR/<test>/ and prints one PASS / FAIL line with the
measured numbers and the criterion.  The reference values and tolerances are in
ref/references.json; see README.md for what each test checks and why.
"""
import argparse
import concurrent.futures as cf
import json
import os
import re
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vett  # noqa: E402

DECKS = os.path.join(HERE, 'decks')
REF = json.load(open(os.path.join(HERE, 'ref', 'references.json')))


class Ctx:
    """What every test needs: the binary, its work directory and the run options."""

    def __init__(self, a, name):
        self.exe = [os.path.abspath(a.exe)]
        self.dir = os.path.join(os.path.abspath(a.workdir), name)
        self.threads = a.threads
        self.quick = a.quick
        self.mpirun = a.mpirun
        self.ref = REF.get(name, {})
        shutil.rmtree(self.dir, ignore_errors=True)
        os.makedirs(self.dir)

    def run(self, deck, sub, over=None, args=(), restart=None, ranks=1):
        """Run deck (a file in decks/) with block/key overrides in DIR/sub; returns the
        run directory.  restart = a restart file to continue from instead."""
        d = os.path.join(self.dir, sub)
        os.makedirs(d, exist_ok=True)
        if restart:
            argv = ['-r', restart] + ['%s=%s' % kv for kv in (over or {}).items()]
        else:
            vett.make_deck(os.path.join(DECKS, deck), os.path.join(d, 'in.athinput'), over)
            argv = ['-i', 'in.athinput']
        rc, wall, log = vett.run(self.exe, d, argv + list(args), ranks=ranks,
                                threads=self.threads, mpirun=self.mpirun)
        if rc != 0:
            raise RuntimeError('%s: return code %s (see %s/run.log)' % (sub, rc, d))
        return d


def rel_drift(x):
    """max |x - x_0| / |x_0| over a history column."""
    return float(np.abs(x - x[0]).max()/abs(x[0]))


def hst_col(path, name):
    lab, a = vett.read_hst(path)
    return a[:, lab.index(name)]


# ---------------------------------------------------------------------------------------
# conservation: closed boxes, energy / angular momentum / mass to round-off

def am_budget(d, base):
    """energy, L_z and mass drifts of a closed sph_shell run (problem/am_hist history).
    dLz: the scheme's total L_z (gas + radiation + the L_z ledger of the clean solve,
    column Lz_tot); dLz_gr: gas + radiation only, for information (on the cubed sphere
    Lrz is the cell vector, not the face momenta the scheme conserves: use dLz)."""
    u = os.path.join(d, base + '.user.hst')
    e = hst_col(u, 'Eg') + hst_col(u, 'Er')
    lzs = np.abs(hst_col(u, 'Lgz')).max()
    lzt = hst_col(u, 'Lz_tot')[1:]   # (row 0: before the first solve, not yet summed)
    lgr = hst_col(u, 'Lgz') + hst_col(u, 'Lrz')
    m = hst_col(os.path.join(d, base + '.hydro.hst'), 'mass')
    return dict(dE=rel_drift(e), dM=rel_drift(m), dLz=float(np.abs(lzt - lzt[0]).max()/lzs),
                dLz_gr=float(np.abs(lgr - lgr[0]).max()/lzs))


def t_conservation_sp(c):
    n = 16 if c.quick else 40
    d = c.run('shell_sp.athinput', 'run', {'time/nlim': n})
    return am_budget(d, 'sp')


def t_conservation_cs(c):
    d = c.run('shell_cs.athinput', 'run', {'time/nlim': 12})
    return am_budget(d, 'cs')


# ---------------------------------------------------------------------------------------
# restart and dt_retry: bitwise

def restart_pair(c, deck, n, nr, over=None):
    """n cycles straight vs nr cycles + restart to n: largest difference of the dumps."""
    o = dict(over or {}, **{'time/nlim': n})
    a = c.run(deck, 'straight', o)
    o.update({'time/nlim': nr, 'output9/file_type': 'rst', 'output9/dcycle': nr})
    b = c.run(deck, 'restarted', o)
    rst = sorted(vett.glob_rst(b))[-1]
    for f in vett.bins(b) + vett.hsts(b):
        os.remove(f)
    c.run(deck, 'restarted', {'time/nlim': n}, restart=os.path.relpath(rst, b))
    for f in vett.glob_rst(b):
        os.remove(f)
    nf, mx, _ = vett.compare_dirs(a, b, last_only=True)
    return dict(files=nf, maxdiff=mx)


def t_restart_sp(c):
    return restart_pair(c, 'shell_sp.athinput', 16, 8)


def t_restart_cs(c):
    return restart_pair(c, 'shell_cs.athinput', 8, 4)


def t_dt_retry(c):
    """a forced dt_retry with factor 1 repeats the step exactly; factor 0.5 runs on."""
    o = {'time/nlim': 16, 'time/dt_retry': 'true'}
    a = c.run('shell_sp.athinput', 'straight', o)
    b = c.run('shell_sp.athinput', 'retry1', dict(o, **{'time/dt_retry_force_cycle': 4,
                                                        'time/dt_retry_force_fac': 1.0}))
    g = c.run('shell_sp.athinput', 'retry05', dict(o, **{'time/dt_retry_force_cycle': 4,
                                                         'time/dt_retry_force_fac': 0.5}))
    nf, mx, _ = vett.compare_dirs(a, b)
    nret = [open(os.path.join(x, 'run.log')).read().count('DT_RETRY cycle=') for x in (b, g)]
    return dict(files=nf, maxdiff=mx, retries_fac1=nret[0], retries_fac05=nret[1])


# ---------------------------------------------------------------------------------------
# layout invariance (tile mode): 1 MeshBlock vs 2 x 2 MeshBlocks, bitwise

def t_layout_tile(c):
    """tile mode with 8 x 8 tiles: one MeshBlock vs 1 x 2 MeshBlocks (along phi, mirror
    theta walls of the ray sweep) and vs 2 x 2 MeshBlocks (vacuum theta walls: the mirror
    walls need one MeshBlock along theta)."""
    o = {'time/nlim': 8, 'mesh/nghost': 3, 'mesh/nx2': 16, 'mesh/nx3': 16,
         'mesh/layout_invariant_mode': 'tile', 'mesh/layout_tile': 8,
         'meshblock/nx2': 16, 'meshblock/nx3': 16}
    out = {}
    for tag, wall, mb2 in (('1x2', 'mirror', 16), ('2x2', 'vacuum', 8)):
        o['rad_vet/vet_gd_theta_wall'] = wall
        a = c.run('shell_sp.athinput', 'mb1_' + wall, dict(o))
        b = c.run('shell_sp.athinput', 'mb' + tag, dict(o, **{'meshblock/nx2': mb2,
                                                              'meshblock/nx3': 8}))
        nf, mx, mr = vett.compare_dirs(a, b)
        out['maxdiff_' + tag], out['maxrel_' + tag] = mx, mr
    return out


# ---------------------------------------------------------------------------------------
# order of accuracy

def succ_orders(fields):
    """successive-difference relative L1 errors and orders of a list of arrays on the
    same grid (coarsest step first)."""
    e = [np.abs(fields[i] - fields[i + 1]).sum()/np.abs(fields[i + 1]).sum()
         for i in range(len(fields) - 1)]
    return e, [float(np.log2(e[i]/e[i + 1])) for i in range(len(e) - 1)]


def t_order_time(c):
    """deck D (spherical shell, kinematic gas, Marshak walls): time step halved 2 (quick)
    or 4 times at a fixed grid; E, F_r and T converge at second order as dt -> 0."""
    dts = ([] if c.quick else [0.0625, 0.03125]) + [0.015625, 0.0078125, 0.00390625]
    S = {'E': [], 'F': [], 'T': []}
    for dt in dts:
        d = c.run('deckD.athinput', 'dt%g' % dt, {'time/dt_fixed': dt})
        m, h = vett.last_bin(d, 'vet'), vett.last_bin(d, 'hydro_w')
        S['E'].append(m['vet_e'])
        S['F'].append(m['vet_f1'])
        S['T'].append(h['eint']/h['dens'])
    out = {}
    for k, v in S.items():
        e, p = succ_orders(v)
        out['err_' + k] = e
        out['order_' + k] = p
    return out


def coarsen(a):
    """average pairs of cells along x1 (the last axis)."""
    return 0.5*(a[..., 0::2] + a[..., 1::2])


def t_order_space(c):
    """linear radiation-hydro wave (Cartesian, periodic, kappa 0.1 per unit length): grid
    and time step refined together; the solution converges at second order."""
    ns = [16, 32, 64, 128]
    sol = {'E': [], 'v': []}
    for n in ns:
        d = c.run('wave_cart.athinput', 'N%d' % n, {'mesh/nx1': n, 'meshblock/nx1': n})
        m, h = vett.last_bin(d, 'vet'), vett.last_bin(d, 'hydro_w')
        sol['E'].append(m['vet_e'].mean(axis=(0, 1)))
        sol['v'].append(h['velx'].mean(axis=(0, 1)))
    out = {}
    for k, v in sol.items():
        # each level on the next coarser grid: differences of neighbouring levels
        e = []
        for i in range(len(v) - 1):
            f = coarsen(v[i + 1])
            e.append(float(np.abs(v[i] - f).sum()/np.abs(f - f.mean()).sum()))
        out['err_' + k] = e
        out['order_' + k] = [float(np.log2(e[i]/e[i + 1])) for i in range(len(e) - 1)]
    return out


# ---------------------------------------------------------------------------------------
# literature tests with reference solutions

def t_nlte(c, eps_tag, deck):
    """two-level-atom NLTE atmosphere (e.g. Jiang 2021 sec. 4.5): J against the exact
    transport solution of the same atmosphere, L1 over 1e-3 < tau < 1e5."""
    d = c.run(deck, 'run')
    m = vett.last_bin(d, 'vet')
    x1 = vett.cell_centres(d)[0]
    E = m['vet_e'].mean(axis=(0, 1))
    tau = 1.0e-3*(np.exp(10.0 - x1) - 1.0)
    t, J = np.loadtxt(os.path.join(HERE, 'ref', 'nlte_exact_eps%s.txt' % eps_tag)).T
    s = (tau > 1e-3) & (tau < 1e5)
    ex = np.interp(np.log(tau[s]), np.log(t), J)
    return dict(L1=float(np.abs(E[s] - ex).sum()/ex.sum()))


def t_nlte_eps1e_1(c):
    return t_nlte(c, '1e-1', 'nlte_eps1e-1.athinput')


def t_nlte_eps1e_4(c):
    return t_nlte(c, '1e-4', 'nlte_eps1e-4.athinput')


def t_thermal_eq(c):
    """uniform gas + radiation out of equilibrium relax at constant total energy
    (e.g. Jiang 2021 sec. 4.1), dt resolving the relaxation time (~3e-6): T(t) against
    the exact ODE at two time steps; the error falls at second order."""
    out = {}
    for tag, dt in (('h', 2.5e-7), ('h2', 1.25e-7)):
        d = c.run('thermal_eq.athinput', 'dt' + tag, {'time/dt_fixed': dt})
        p = vett.parse_deck_text(open(os.path.join(d, 'in.athinput')).read())
        cl, ka = float(p['rad_vet']['c_light']), float(p['rad_vet']['kappa_p'])
        arad = float(p['rad_vet']['arad'])
        rho, gm1 = float(p['problem']['gas_rho']), float(p['hydro']['gamma']) - 1.0
        T0, E0 = float(p['problem']['gas_eint'])*gm1/rho, float(p['problem']['e_rad'])
        tot = rho*T0/gm1 + E0
        errT = errE = 0.0
        for fm, fh in list(zip(vett.bins(d, 'vet'), vett.bins(d, 'hydro_w')))[1:]:
            m, h = vett.read_bin(fm), vett.read_bin(fh)
            Tex = vett.relax_ode(T0, tot, rho, gm1, cl*ka*rho, m['time'], arad=arad)
            T = (h['eint']*gm1/h['dens']).mean()
            errT = max(errT, abs(T/Tex - 1.0))
            errE = max(errE, abs(m['vet_e'].mean() + rho*T/gm1 - tot)/tot)
        out['errT_' + tag], out['dEtot_' + tag] = float(errT), float(errE)
    out['order'] = float(np.log2(out['errT_h']/out['errT_h2']))
    return out


# ---------------------------------------------------------------------------------------
# well-balancing: a radiative hydrostatic envelope stays at rest

def t_wb_sp(c):
    d = c.run('rhse_sp.athinput', 'run')
    v = hst_col(os.path.join(d, 'rh.user.hst'), 'vcs_max')
    return dict(vmax_over_cs=float(v.max()))


def t_wb_cs(c):
    d = c.run('rhse_cs.athinput', 'run')
    v = hst_col(os.path.join(d, 'rh.user.hst'), 'vcs_max')
    return dict(vmax_over_cs=float(v.max()))


# ---------------------------------------------------------------------------------------
# radiation-acoustic waves in every trapping regime

RADWAVE = (  # tag, P_rad/P_gas, x = tau_cell c_eq/c (grid-scale trapping), signal-speed taper
    ('thin', 10.0, 0.01, 'true'), ('thin_old', 10.0, 0.01, 'false'),
    ('mid', 100.0, 0.1, 'true'), ('one', 10.0, 1.0, 'true'), ('trapped', 1.0, 100.0, 'true'))


def radwave_fit(d, amp):
    """omega of the k = 1 mode of a radwave run: drho = A e^{-i w t} + B e^{-g t} (the
    acoustic and the non-propagating entropy mode), dv = V e^{-i w t}; joint fit."""
    k, ts, ar, av = 2.0*np.pi, [], [], []
    for f in vett.bins(d, 'hydro_w'):
        h = vett.read_bin(f)
        e = np.exp(-1j*k*vett.cell_centres(d)[0])
        ts.append(h['time'])
        ar.append(2.0*np.mean((h['dens'].mean(axis=(0, 1)) - 1.0)*e)/amp)
        av.append(2.0*np.mean(h['velx'].mean(axis=(0, 1))*e)/amp)
    ts, ar, av = np.array(ts), np.array(ar), np.array(av)
    w0 = complex(-np.polyfit(ts, np.unwrap(np.angle(av)), 1)[0],
                 np.polyfit(ts, np.log(np.abs(av)), 1)[0])

    def res(p):
        w, A, V, B = p[0] + 1j*p[1], p[2] + 1j*p[3], p[4] + 1j*p[5], p[6] + 1j*p[7]
        ea, sc = np.exp(-1j*w*ts), np.exp(w.imag*ts)   # relative to the decaying envelope
        return np.concatenate([((A*ea + B*np.exp(-p[8]*ts) - ar)/sc/abs(ar[0])).view(float),
                               ((V*ea - av)/sc/abs(av[0])).view(float)])
    p, _ = vett.lm_fit(res, [w0.real, w0.imag, ar[0].real, ar[0].imag, av[0].real,
                            av[0].imag, 0.0, 0.0, 0.1])
    return complex(p[0], p[1])


def t_radwave(c):
    """linear radiation-acoustic wave (64 x 4 cells, amplitude 1e-6, the exact eigenmode,
    2 periods, cfl 0.3): phase speed and damping rate against the analytic root, from the
    optically thin (x = 0.01) to the trapped (x = 100) regime; 'thin_old' runs without
    hydro/rad_signal_speed_diff, for the time-step gain."""
    cl, n, amp = 1.0e3, 64, 1.0e-6
    out = {'err_cph': [], 'err_damp': [], 'damp': [], 'warnings': 0}
    dt = {}
    for tag, prat, x, diff in RADWAVE:
        ceq = vett.radwave_root(prat, 1.0, cl)[3]
        kap = x*cl*n/ceq                     # tau_cell = kappa h = x c / c_eq
        w, ar, eig, _ = vett.radwave_root(prat, kap, cl)
        tl = 2.0*2.0*np.pi/w.real
        o = {'mesh/nx1': n, 'meshblock/nx1': n, 'time/cfl_number': 0.3, 'time/tlim': tl,
             'hydro/rad_signal_speed': 'true', 'hydro/rad_signal_speed_diff': diff,
             'rad_vet/c_light': cl, 'rad_vet/arad': ar, 'rad_vet/implicit_cfl': -1.0,
             'problem/radwave_amp': amp, 'output2/dt': tl/32}
        o.update({'rad_vet/kappa_' + s: kap for s in 'pef'})
        o.update({'problem/' + kk: v for kk, v in eig.items()})
        d = c.run('wave_cart.athinput', tag, o)
        out['warnings'] += vett.warnings(open(os.path.join(d, 'run.log')).read())
        h = vett.last_bin(d, 'hydro_w')
        dt[tag] = h['time']/h['cycle']
        if tag.endswith('_old'):
            continue
        wm = radwave_fit(d, amp)
        out['err_cph'].append(float(abs(wm.real/w.real - 1.0)))
        out['err_damp'].append(float(abs(wm.imag - w.imag)/wm.real))   # per radian
        out['damp'].append(float(-wm.imag/w.real))                     # > 0: decays
    out['dt_gain'] = float(dt['thin']/dt['thin_old'])
    return out


# ---------------------------------------------------------------------------------------
# phi mirror symmetry on narrow wedges (vet_gd_wall_sym)

def t_wedge_sym(c):
    """the closed shell at rest and spherically symmetric (no rotation, no lateral mode) on
    phi wedges with the production Gauss-Legendre 7 x 18 set (azimuthal period 20 deg):
    a 20-deg wedge (exact phi-wall maps) and narrow 3.7 and 7.3 deg wedges with
    vet_gd_wall_sym = true; the wedge-mean F_phi must vanish: fphi = |sum F3| / sum |F1|
    after 8 cycles.  The narrow wedge without the key must stop with the wedge FATAL."""
    o = {'time/nlim': 8, 'problem/gas_rot': 0.0, 'problem/lat_amp': 0.0,
         'rad_vet/vet_gd_gl_nmu': 7, 'rad_vet/vet_gd_gl_nphi': 18,
         'rad_vet/vet_gd_gl_stagger': 'true', 'rad_vet/vet_gd_wall_interp': 'true'}
    out = {'fphi': [], 'warnings': 0}
    for tag, w, sym in (('w20', np.pi/9, 'false'), ('w3.7', 0.064, 'true'),
                        ('w7.3', 0.128, 'true')):
        d = c.run('shell_sp.athinput', tag, dict(o, **{'mesh/x3max': w,
                                                        'rad_vet/vet_gd_wall_sym': sym}))
        m = vett.last_bin(d, 'vet')
        out['fphi'].append(float(abs(m['vet_f3'].sum())/np.abs(m['vet_f1']).sum()))
        out['warnings'] += vett.warnings(open(os.path.join(d, 'run.log')).read())
    try:
        c.run('shell_sp.athinput', 'w3.7_nosym', dict(o, **{'mesh/x3max': 0.064}))
        out['wedge_fatal'] = 0
    except RuntimeError:
        log = open(os.path.join(c.dir, 'w3.7_nosym', 'run.log')).read()
        out['wedge_fatal'] = int('vet_gd_wall_sym' in log and 'FATAL' in log)
    return out


# ---------------------------------------------------------------------------------------
# floors with a general EOS: the floor energy is paid from the radiation

def t_floor_pay(c):
    """the closed shell with hydro/eos = general (cgs units = code units) and a
    temperature floor just above the gas temperature, so the floor raises T in most cells
    every cycle; hydro/floor_pay = rad_column pays that energy from the radiation of the
    same column.  dE: total-energy drift from cycle 1 on (the payment of the initial
    conversion's floor energy comes at the top of cycle 1, after the t = 0 output: dE0 is
    that jump relative to E); paid / unpaid: floor energy over the run (floor_audit)."""
    d = c.run('shell_sp.athinput', 'run', {
        'hydro/eos': 'general', 'units/length_cgs': 1.0, 'units/mass_cgs': 1.0,
        'units/time_cgs': 1.0, 'hydro/tfloor': 1.00001e-4, 'hydro/floor_pay': 'rad_column',
        'hydro/floor_pay_skip': 'true', 'hydro/dfloor_flux_limit': 'true',
        'hydro/floor_audit': 1})
    u = os.path.join(d, 'sp.user.hst')
    e = hst_col(u, 'Eg') + hst_col(u, 'Er')
    log = open(os.path.join(d, 'run.log')).read()
    pay = [[float(x) for x in m] for m in
           re.findall(r'floor_audit cycle \d+ pay: paid (\S+) unpaid (\S+)', log)]
    return dict(dE=rel_drift(e[1:]), dE0=float(abs(e[1] - e[0])/e[0]),
                paid=sum(p[0] for p in pay), unpaid=sum(p[1] for p in pay),
                warnings=vett.warnings(log))


# ---------------------------------------------------------------------------------------
# linear momentum: a report, not a pass / fail criterion

def t_momentum_box(c):
    """periodic 3-D Cartesian box (8^3, the radwave along the xyz diagonal, amplitude 1e-3,
    8 cycles): total linear momentum = gas (hydro_u, cell sums) + radiation (the face
    momenta F_f / c^2 of the vet_face dump; each cell's left faces = every face once in a
    periodic box), every cycle.  dP_cyc: per cycle, the largest component of |P - P_0| /
    sum |p_gas| at t = 0.  The scheme does not conserve P in general (see README): the
    numbers are reported; the test passes when the run completes with finite numbers and
    no CLEAN WARNING."""
    n, cl = 8, 10.0
    o = {'mesh/nghost': 3, 'mesh/nx1': n, 'mesh/nx2': n, 'mesh/nx3': n,
         'meshblock/nx1': n, 'meshblock/nx2': n, 'meshblock/nx3': n, 'mesh/x3max': 1.0,
         'mesh/layout_invariant_mode': 'tile', 'mesh/layout_tile': 8,
         'time/tlim': 1e30, 'time/nlim': 8, 'problem/radwave_dir': 'xyz',
         'problem/radwave_amp': 1e-3, 'problem/radwave_eig': 'false',
         'output1/variable': 'vet_face', 'output1/dcycle': 1,
         'output2/variable': 'hydro_u', 'output2/dcycle': 1}
    d = c.run('wave_cart.athinput', 'run', o)
    dv = 1.0/n**3
    P, s0 = {}, None   # by cycle (the final dump repeats the last cycle)
    for fa, fh in zip(vett.bins(d, 'vet_face'), vett.bins(d, 'hydro_u')):
        A, H = vett.read_bin(fa), vett.read_bin(fh)
        pr = np.array([A['vet_fx%d' % k].sum() for k in (1, 2, 3)])*dv/cl**2
        pg = np.array([H['mom%d' % k].sum() for k in (1, 2, 3)])*dv
        if s0 is None:
            s0 = sum(np.abs(H['mom%d' % k]).sum() for k in (1, 2, 3))*dv
        P[A['cycle']] = pg + pr
    P = [P[k] for k in sorted(P)]
    dp = [float(np.abs(p - P[0]).max()/s0) for p in P[1:]]
    if len(dp) != o['time/nlim'] or not np.all(np.isfinite(dp)):
        raise RuntimeError('momentum_box: %d cycles read, dP %s' % (len(dp), dp))
    w = vett.warnings(open(os.path.join(d, 'run.log')).read())
    return dict(dP_max=max(dp), dP_cyc=dp, warnings=w)


# ---------------------------------------------------------------------------------------

TESTS = [  # name, function, in --quick
    ('conservation_sp', t_conservation_sp, True),
    ('conservation_cs', t_conservation_cs, True),
    ('restart_sp', t_restart_sp, True),
    ('restart_cs', t_restart_cs, False),
    ('dt_retry', t_dt_retry, True),
    ('layout_tile', t_layout_tile, True),
    ('order_time', t_order_time, True),
    ('order_space', t_order_space, False),
    ('nlte_eps1e-1', t_nlte_eps1e_1, True),
    ('nlte_eps1e-4', t_nlte_eps1e_4, False),
    ('thermal_eq', t_thermal_eq, False),
    ('wb_sp', t_wb_sp, True),
    ('wb_cs', t_wb_cs, False),
    ('momentum_box', t_momentum_box, True),
    ('radwave', t_radwave, False),
    ('wedge_sym', t_wedge_sym, True),
    ('floor_pay', t_floor_pay, True),
]


def check(name, res, ref):
    """PASS / FAIL of a test's measured values against its criteria in references.json:
    {"max": {k: bound}, "min": {k: bound}, "equal": {k: value},
     "near": {k: [reference, relative tolerance]}}; list-valued results are checked
    element by element, or by their last element only with "last": true ("near" always
    takes the last element: the finest level of a ladder)."""
    bad = []
    for k, b in ref.get('max', {}).items():
        v = np.atleast_1d(res[k])
        if ref.get('last'):
            v = v[-1:]
        if not np.all(v <= b):
            bad.append('%s %s > %g' % (k, fmt(res[k]), b))
    for k, b in ref.get('min', {}).items():
        v = np.atleast_1d(res[k])
        if ref.get('last'):
            v = v[-1:]
        if not np.all(v >= b):
            bad.append('%s %s < %g' % (k, fmt(res[k]), b))
    for k, b in ref.get('equal', {}).items():
        if res[k] != b:
            bad.append('%s %s != %s' % (k, fmt(res[k]), b))
    for k, (r, tol) in ref.get('near', {}).items():
        v = np.atleast_1d(res[k])[-1]
        if not abs(v - r) <= tol*abs(r):
            bad.append('%s %s not within %g of %g' % (k, fmt(res[k]), tol, r))
    return bad


def fmt(v):
    if isinstance(v, (list, tuple)):
        return '[' + ' '.join(fmt(x) for x in v) + ']'
    if isinstance(v, float):
        return '%.3g' % v
    return str(v)


def one(a, name, fn):
    t0 = time.time()
    c = Ctx(a, name)
    try:
        res = fn(c)
        bad = check(name, res, c.ref) if c.ref else ['no reference in references.json']
    except Exception as e:  # a failed run is a FAIL with its reason
        res, bad = {}, ['ERROR %s' % e]
    if not bad and not a.keep:
        shutil.rmtree(c.dir, ignore_errors=True)
    return name, res, bad, time.time() - t0


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--exe', required=True, help='the athena binary (CPU build)')
    p.add_argument('--quick', action='store_true', help='the short subset (11 tests)')
    p.add_argument('--jobs', type=int, default=max(1, (os.cpu_count() or 4)//4),
                   help='tests run at the same time')
    p.add_argument('--threads', type=int, default=4, help='OpenMP threads per run')
    p.add_argument('--only', default='', help='comma-separated test names')
    p.add_argument('--workdir', default='vet_clean_runs', help='where the runs go')
    p.add_argument('--mpirun', default='mpirun', help='MPI launcher (multi-rank tests)')
    p.add_argument('--keep', action='store_true', help='keep the run directories of passes')
    p.add_argument('--record', default='', help='write the measured values to this JSON')
    a = p.parse_args()
    sel = [t for t in TESTS if (not a.quick or t[2]) and
           (not a.only or t[0] in a.only.split(','))]
    print('clean VET test suite: %d tests, %s, %d at a time x %d threads, binary %s' % (
        len(sel), 'quick' if a.quick else 'full', a.jobs, a.threads, a.exe))
    t0, rows = time.time(), []
    with cf.ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(one, a, n, f) for n, f, _ in sel]
        for f in cf.as_completed(futs):
            name, res, bad, wall = f.result()
            rows.append((name, res, bad, wall))
            print('%-4s %-16s %6.0f s  %s%s' % (
                'PASS' if not bad else 'FAIL', name, wall,
                ' '.join('%s %s' % (k, fmt(v)) for k, v in res.items()),
                ('  <- ' + '; '.join(bad)) if bad else ''), flush=True)
    nf = sum(1 for r in rows if r[2])
    print('%d passed, %d failed, %.0f s' % (len(rows) - nf, nf, time.time() - t0))
    if a.record:
        json.dump({r[0]: dict(r[1], wall=r[3]) for r in rows}, open(a.record, 'w'), indent=1)
    sys.exit(1 if nf else 0)


if __name__ == '__main__':
    main()
