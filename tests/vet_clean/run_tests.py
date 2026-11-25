#!/usr/bin/env python3
"""run_tests.py: the clean VET radiation test suite (CPU; Python standard library + numpy).

usage: run_tests.py --exe PATH [--quick] [--jobs N] [--threads T] [--only NAME,...]
                    [--workdir DIR] [--mpirun CMD] [--keep] [--record FILE]
                    [--keys 'block/key=value ...']

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
        self.star = [os.path.abspath(a.star_exe)] if a.star_exe else None
        self.dir = os.path.join(os.path.abspath(a.workdir), name)
        self.threads = a.threads
        self.quick = a.quick
        self.mpirun = a.mpirun
        self.keys = dict(kv.split('=', 1) for kv in a.keys.split())
        self.ref = REF.get(name, {})
        shutil.rmtree(self.dir, ignore_errors=True)
        os.makedirs(self.dir)

    def run(self, deck, sub, over=None, args=(), restart=None, ranks=1):
        """Run deck (a file in decks/) with block/key overrides in DIR/sub; returns the
        run directory.  restart = a restart file to continue from instead.  The --keys
        overrides go into every fresh run (the test's own overrides win); a restart takes
        only the test's own overrides, so the --keys never override what the restart file
        stored from the run it continues."""
        own = dict(over or {})
        over = dict(self.keys, **own)
        d = os.path.join(self.dir, sub)
        os.makedirs(d, exist_ok=True)
        if restart:
            argv = ['-r', restart] + ['%s=%s' % kv for kv in own.items()]
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


# L_z-exact caps and cuts (read only when named)

LZ_KEYS = {'rad_vet/clean_stag_fcone': 'true', 'rad_vet/clean_stag_fcone_lz': 'true',
           'rad_vet/clean_stag_fcone_emax': 'true', 'rad_vet/clean_wb_mcap': 'true',
           'rad_vet/clean_lz_spin': 'true', 'hydro/eos': 'general',
           'hydro/floor_from_ekin': 'true', 'hydro/vceil': 5e-3, 'units/length_cgs': 1.0,
           'units/mass_cgs': 1.0, 'units/time_cgs': 1.0}
LZ_FC = {'rad_vet/floorcap_rho': 1e30, 'rad_vet/floorcap_v': 5e-3,
         'rad_vet/floorcap_lz': 'true'}
LZ_SP = dict(LZ_FC, **{'rad_vet/clean_stag_fwb_ring': 1})


def t_conservation_lz(c):
    """the closed rotating shells with every L_z-exact cap / cut on (clean_stag_fcone_lz,
    clean_lz_spin, clean_stag_fcone_emax, clean_wb_mcap, floor_from_ekin with vceil 5e-3
    below the shell's speeds; the floor-gas cap of every cell, floorcap_lz, on both
    meshes, on spherical polar with fwb_ring): energy, mass and L_z drifts as conservation_sp, plus
    dLz_gr (gas + radiation, no ledger); restart (8 + 8 = 16 cycles) and a forced dt_retry
    with factor 1 at cycle 4 with the same keys: bitwise."""
    sp = dict(LZ_KEYS, **LZ_SP)
    out = {'sp_' + k: v for k, v in
           am_budget(c.run('shell_sp.athinput', 'sp', dict(sp, **{'time/nlim': 40})),
                     'sp').items()}
    out.update({'cs_' + k: v for k, v in
                am_budget(c.run('shell_cs.athinput', 'cs',
                                dict(LZ_KEYS, **dict(LZ_FC, **{'time/nlim': 12}))),
                          'cs').items()})
    out['restart_maxdiff'] = restart_pair(c, 'shell_sp.athinput', 16, 8, sp)['maxdiff']
    o = dict(sp, **{'time/nlim': 16, 'time/dt_retry': 'true'})
    a = c.run('shell_sp.athinput', 'r_straight', o)
    b = c.run('shell_sp.athinput', 'r_retry1', dict(o, **{'time/dt_retry_force_cycle': 4,
                                                          'time/dt_retry_force_fac': 1.0}))
    out['retry_maxdiff'] = vett.compare_dirs(a, b)[1]
    out['warnings'] = sum(vett.warnings(open(os.path.join(c.dir, t, 'run.log')).read())
                          for t in ('sp', 'cs', 'r_straight', 'r_retry1'))
    return out


def t_floor_dt(c):
    """hydro/dt_floor_courant = 0.6 on the rotating shell with floor gas (floorcap_rho = 100,
    the cells that thin below the IC density, as floorcap_spin): floor gas takes dt <= C /
    sum_d (|v_d| + c_d)/h_d instead of cfl_number per direction.  dt_gain: mean dt with the
    key over dt without it (reported); energy as conservation_sp; restart (8 + 8) and a
    forced dt_retry with the key bitwise."""
    fc = dict(LZ_FC, **{'rad_vet/floorcap_rho': 100.0, 'hydro/dt_floor_courant': 0.6})
    on = c.run('shell_sp.athinput', 'on', dict(fc, **{'time/nlim': 16}))
    off = c.run('shell_sp.athinput', 'off', dict(fc, **{'time/nlim': 16, 'hydro/dt_floor_courant': 0.0}))
    dt = [hst_col(os.path.join(d, 'sp.user.hst'), 'dt')[1:] for d in (on, off)]
    out = dict(dE=am_budget(on, 'sp')['dE'], dt_gain=float(np.mean(dt[0])/np.mean(dt[1])))
    out['restart_maxdiff'] = restart_pair(c, 'shell_sp.athinput', 16, 8, fc)['maxdiff']
    o = dict(fc, **{'time/nlim': 16, 'time/dt_retry': 'true'})
    a = c.run('shell_sp.athinput', 'r_straight', o)
    b = c.run('shell_sp.athinput', 'r_retry1', dict(o, **{'time/dt_retry_force_cycle': 4,
                                                          'time/dt_retry_force_fac': 1.0}))
    out['retry_maxdiff'] = vett.compare_dirs(a, b)[1]
    out['warnings'] = sum(vett.warnings(open(os.path.join(c.dir, t, 'run.log')).read())
                          for t in ('on', 'r_straight', 'r_retry1'))
    return out


def t_floorcap_spin(c):
    """floorcap_lz_spin with the conservation_lz keys: the floor-gas cap scales the whole
    velocity and gives the L_z it removes to the gas at or above floorcap_rho as one rigid
    rotation.  floorcap_rho = gas_rho (sp 100, cs 1), so the cells that thin below the IC
    density are floor gas, capped below the shells' speeds (sp 5e-3, cs 0.05: a cap of the
    cs shell's 0.3 at 5e-3 puts so much kinetic energy into E that the first stage solve
    fails), and the rest take the rotation.  Drifts as conservation_lz on both meshes
    (dLz_gr without the ledger), n_v > 0 capped cells; restart (8 + 8) and a forced
    dt_retry bitwise; 1 vs 1 x 2 MeshBlocks in tile mode bitwise.  nospin_sp_dLz_gr: the
    same with floorcap_lz off, for information."""
    fc = dict(LZ_FC, **{'rad_vet/floorcap_rho': 100.0, 'rad_vet/floorcap_lz_spin': 'true',
                        'rad_vet/floorcap_report': 1})
    sp = dict(LZ_KEYS, **dict(fc, **{'rad_vet/clean_stag_fwb_ring': 1}))
    d = c.run('shell_sp.athinput', 'sp', dict(sp, **{'time/nlim': 40}))
    out = {'sp_' + k: v for k, v in am_budget(d, 'sp').items()}
    n = re.findall(r'floorcap cycle \d+ \(cumulative, all calls\): n_v (\S+)',
                   open(os.path.join(d, 'run.log')).read())
    out['n_v'] = float(n[-1]) if n else 0.0
    out.update({'cs_' + k: v for k, v in
                am_budget(c.run('shell_cs.athinput', 'cs',
                                dict(LZ_KEYS, **dict(fc, **{'time/nlim': 12,
                                                          'rad_vet/floorcap_rho': 1.0,
                                                          'rad_vet/floorcap_v': 0.05}))),
                          'cs').items()})
    o = dict(sp, **{'rad_vet/floorcap_lz': 'false', 'rad_vet/floorcap_lz_spin': 'false',
                    'time/nlim': 40})
    out['nospin_sp_dLz_gr'] = am_budget(c.run('shell_sp.athinput', 'nospin', o),
                                        'sp')['dLz_gr']
    out['restart_maxdiff'] = restart_pair(c, 'shell_sp.athinput', 16, 8, sp)['maxdiff']
    o = dict(sp, **{'time/nlim': 16, 'time/dt_retry': 'true'})
    a = c.run('shell_sp.athinput', 'r_straight', o)
    b = c.run('shell_sp.athinput', 'r_retry1', dict(o, **{'time/dt_retry_force_cycle': 4,
                                                          'time/dt_retry_force_fac': 1.0}))
    out['retry_maxdiff'] = vett.compare_dirs(a, b)[1]
    o = dict(sp, **{'time/nlim': 8, 'mesh/nghost': 3, 'mesh/nx2': 16, 'mesh/nx3': 16,
                    'mesh/layout_invariant_mode': 'tile', 'mesh/layout_tile': 8,
                    'meshblock/nx2': 16, 'meshblock/nx3': 16,
                    'rad_vet/vet_gd_theta_wall': 'mirror'})
    a = c.run('shell_sp.athinput', 'mb1', o)
    b = c.run('shell_sp.athinput', 'mb1x2', dict(o, **{'meshblock/nx3': 8}))
    out['layout_maxdiff'] = vett.compare_dirs(a, b)[1]
    out['warnings'] = sum(vett.warnings(open(os.path.join(c.dir, t, 'run.log')).read())
                          for t in ('sp', 'cs', 'r_straight', 'r_retry1', 'mb1', 'mb1x2'))
    return out


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

# star settings the layout test cannot take: the 2 x 2 difference (bound 1e-8) is 7e-8 on
# the 60 deg wedge with GL 7 x 18 and 8e-8 with the budget norms; dfloor_keep_velocity
# creates energy where the density floor acts (no floor_pay with this ideal gas)
LAYOUT_OFF = {'rad_vet/clean_clos1_budget': 'false', 'rad_vet/clean_bud_maxnorm': 'false',
              'rad_vet/vet_gd_wall_interp': 'false', 'rad_vet/vet_gd_wall_sym': 'false',
              'mesh/x3max': 1.5707963267948966, 'rad_vet/vet_gd_gl_nmu': 16,
              'rad_vet/vet_gd_gl_nphi': 32, 'hydro/dfloor_keep_velocity': 'false'}


def t_layout_tile(c):
    """tile mode with 8 x 8 tiles: one MeshBlock vs 1 x 2 MeshBlocks (along phi, mirror
    theta walls of the ray sweep) and vs 2 x 2 MeshBlocks (vacuum theta walls: the mirror
    walls need one MeshBlock along theta).  dfl: 1 vs 1 x 2 again with
    hydro/dfloor_flux_limit and dfloor just below the gas density (the outflow limiter
    acts in every cell; sdirk2's stage-2 donor across a MeshBlock face is a ghost), and
    that run's energy and L_z drifts (am_budget).  fofc: 1 vs 1 x 2 with hydro/fofc and
    that dfloor (fo_fofc first-order cell-stages > 0), with its energy and L_z drifts
    (dM_fofc is the density floor's own mass, reported)."""
    o = {'time/nlim': 8, 'mesh/nghost': 3, 'mesh/nx2': 16, 'mesh/nx3': 16,
         'mesh/layout_invariant_mode': 'tile', 'mesh/layout_tile': 8,
         'meshblock/nx2': 16, 'meshblock/nx3': 16}
    o.update(LAYOUT_OFF)
    out = {}
    for tag, wall, mb2 in (('1x2', 'mirror', 16), ('2x2', 'vacuum', 8)):
        o['rad_vet/vet_gd_theta_wall'] = wall
        a = c.run('shell_sp.athinput', 'mb1_' + wall, dict(o))
        b = c.run('shell_sp.athinput', 'mb' + tag, dict(o, **{'meshblock/nx2': mb2,
                                                              'meshblock/nx3': 8}))
        nf, mx, mr = vett.compare_dirs(a, b)
        out['maxdiff_' + tag], out['maxrel_' + tag] = mx, mr
    o.update({'rad_vet/vet_gd_theta_wall': 'mirror', 'hydro/dfloor_flux_limit': 'true',
              'hydro/dfloor': 99.9})
    a = c.run('shell_sp.athinput', 'mb1_dfl', dict(o))
    b = c.run('shell_sp.athinput', 'mb1x2_dfl', dict(o, **{'meshblock/nx3': 8}))
    out['maxdiff_dfl'] = vett.compare_dirs(a, b)[1]
    out.update({k + '_dfl': v for k, v in am_budget(b, 'sp').items() if k != 'dLz_gr'})
    # (fofc: the same with hydro/fofc and that dfloor without the outflow limiter, so the
    # trial state floors and FOFC replaces fluxes; its trial reads the u0 ghosts too; an
    # event log every 2 cycles, whose rows reset its interval counter: fo_fofc is the
    # run's total all the same)
    o.update({'hydro/dfloor_flux_limit': 'false', 'hydro/fofc': 'true',
              'output9/file_type': 'log', 'output9/dcycle': 2})
    a = c.run('shell_sp.athinput', 'mb1_fofc', dict(o))
    b = c.run('shell_sp.athinput', 'mb1x2_fofc', dict(o, **{'meshblock/nx3': 8}))
    out['maxdiff_fofc'] = vett.compare_dirs(a, b)[1]
    out.update({k + '_fofc': v for k, v in am_budget(b, 'sp').items() if k != 'dLz_gr'})
    lg = open(os.path.join(b, 'run.log')).read()
    out['fo_fofc'] = float(re.findall(r'FOFC first-order cell-stages=(\S+)', lg)[-1])
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
    or 4 times at a fixed grid; successive orders of E, F_r and T."""
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


def t_order_stiff(c):
    """linear radiation-acoustic wave in the trapped regime (P_rad/P_gas 1, x = 100; 32
    cells, one period, strict solves), time step T/64 ... T/512, with time2_gtime_slope:
    the stiff exchange of sdirk2 converges at second order."""
    cl, n = 1.0e3, 32
    ceq = vett.radwave_root(1.0, 1.0, cl)[3]
    kap = 100.0*cl*n/ceq
    w, ar, eig, _ = vett.radwave_root(1.0, kap, cl)
    tl = 2.0*np.pi/w.real
    o = {'mesh/nx1': n, 'meshblock/nx1': n, 'time/cfl_number': 10.0, 'time/tlim': tl,
         'hydro/rad_signal_speed': 'true', 'hydro/rad_signal_speed_diff': 'true',
         'rad_vet/c_light': cl, 'rad_vet/arad': ar, 'rad_vet/implicit_cfl': -1.0,
         'rad_vet/clean_conv': 'strict', 'rad_vet/clean_tol': 1e-14,
         'rad_vet/clean_tol_nl': 1e-13, 'rad_vet/time2_gtime_slope': 'true'}
    o.update({'rad_vet/kappa_' + s: kap for s in 'pef'})
    o.update({'problem/' + k: v for k, v in eig.items()})
    S = {'E': [], 'v': []}
    for m in (64, 128, 256, 512):
        d = c.run('wave_cart.athinput', 'T%d' % m, dict(o, **{'time/dt_fixed': tl/m}))
        S['E'].append(vett.last_bin(d, 'vet')['vet_e'])
        S['v'].append(vett.last_bin(d, 'hydro_w')['velx'])
    out = {}
    for k, v in S.items():
        out['err_' + k], out['order_' + k] = succ_orders(v)
    return out


def t_order_space(c):
    """linear radiation-hydro wave (Cartesian, periodic, kappa 0.1 per unit length): grid
    and time step refined together; orders of successive differences of E and v."""
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


def relax_errors(d):
    """(max |T/T_exact - 1|, max total-energy drift) over the dumps after t = 0 of a
    uniform relaxation run d (vet_test = equil) against the exact ODE."""
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
    return float(errT), float(errE)


def t_thermal_eq(c):
    """uniform gas + radiation out of equilibrium relax at constant total energy
    (e.g. Jiang 2021 sec. 4.1), dt resolving the relaxation time (~3e-6): T(t) against
    the exact ODE at two time steps and the order of its error."""
    out = {}
    for tag, dt in (('h', 2.5e-7), ('h2', 1.25e-7)):
        d = c.run('thermal_eq.athinput', 'dt' + tag, {'time/dt_fixed': dt})
        out['errT_' + tag], out['dEtot_' + tag] = relax_errors(d)
    out['order'] = float(np.log2(out['errT_h']/out['errT_h2']))
    return out


# ---------------------------------------------------------------------------------------
# well-balancing: a radiative hydrostatic envelope stays at rest

# the rotating-set arm: the set turns every dt_rot (sp: R136's GL 11 x 12 about z; cs:
# RSG's GL 10 x 12 about a 3-D axis), the fresh default EMA ring with the cross-fade, and
# vet_rotate_avg_ema_wb, so that a set change leaves a static state's closure unchanged
WB_ROT = {'rad_vet/vet_rotate_avg_ema_wb': 'true', 'rad_vet/vet_gd_gl_nphi': 12,
          'rad_vet/vet_gd_gl_stagger': 'true'}
WB_ROT_SP = dict(WB_ROT, **{'rad_vet/vet_gd_rotate_dt': 2.0, 'rad_vet/vet_gd_gl_nmu': 11,
                            'rad_vet/vet_gd_wall_interp': 'true',
                            'rad_vet/vet_gd_wall_sym': 'true'})
WB_ROT_CS = dict(WB_ROT, **{'rad_vet/vet_gd_rotate_dt': 0.5, 'rad_vet/vet_gd_gl_nmu': 10})


def t_wb(c, deck, rot):
    """max v / c_s of the envelope with the fixed set and with the rotating one"""
    out = {}
    for tag, o in (('', None), ('_rot', rot)):
        d = c.run(deck, 'run' + tag, o)
        v = hst_col(os.path.join(d, 'rh.user.hst'), 'vcs_max')
        out['vmax_over_cs' + tag] = float(v.max())
    return out


def t_wb_sp(c):
    return t_wb(c, 'rhse_sp.athinput', WB_ROT_SP)


def t_wb_cs(c):
    return t_wb(c, 'rhse_cs.athinput', WB_ROT_CS)


# ---------------------------------------------------------------------------------------
# radiation-acoustic waves in every trapping regime

RADWAVE = (  # tag, P_rad/P_gas, x = tau_cell c_eq/c (grid-scale trapping), signal-speed taper
    ('thin', 10.0, 0.01, 'true'), ('thin_old', 10.0, 0.01, 'false'),
    ('mid', 100.0, 0.1, 'true'), ('one', 10.0, 1.0, 'true'), ('trapped', 1.0, 100.0, 'true'),
    ('mid_phase', 100.0, 0.1, 'phase'), ('one_phase', 10.0, 1.0, 'phase'))


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
    """linear radiation-acoustic wave (80 x 4 cells, amplitude 1e-3, the exact eigenmode,
    2 periods, cfl 0.3; the thin wave decays by e^-9, so at 1e-6 it ends below the solver
    tolerance): phase speed and damping rate against the analytic root, from the optically
    thin (x = 0.01 on 64 cells) to the trapped (x = 100) regime; 'thin_old' runs without
    hydro/rad_signal_speed_diff, for the time-step gain; '*_phase' with
    hydro/rad_signal_speed_phase (in-phase weight 1 - exp(-x^2))."""
    cl, n, amp = 1.0e3, 80, 1.0e-3
    out = {'err_cph': [], 'err_damp': [], 'damp': [], 'warnings': 0}
    dt = {}
    for tag, prat, x, diff in RADWAVE:
        ceq = vett.radwave_root(prat, 1.0, cl)[3]
        kap = x*cl*64/ceq                    # tau_cell = kappa h = x c / c_eq on 64 cells
        w, ar, eig, _ = vett.radwave_root(prat, kap, cl)
        tl = 2.0*2.0*np.pi/w.real
        o = {'mesh/nx1': n, 'meshblock/nx1': n, 'time/cfl_number': 0.3, 'time/tlim': tl,
             'hydro/rad_signal_speed': 'true',
             'hydro/rad_signal_speed_diff': 'false' if diff == 'false' else 'true',
             'rad_vet/clean_clos1_budget': 'false', 'rad_vet/clean_bud_maxnorm': 'false',
             'rad_vet/c_light': cl, 'rad_vet/arad': ar, 'rad_vet/implicit_cfl': -1.0,
             'problem/radwave_amp': amp, 'output2/dt': tl/32}
        if diff == 'phase':
            o['hydro/rad_signal_speed_phase'] = 'true'
        else:
            o['hydro/rad_signal_speed_phase'] = 'false'
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
    out['dt_gain_phase'] = [float(dt['mid_phase']/dt['mid']), float(dt['one_phase']/dt['one'])]
    return out


# ---------------------------------------------------------------------------------------
# phi mirror symmetry on narrow wedges (vet_gd_wall_sym)

def t_wedge_sym(c):
    """the closed shell at rest and spherically symmetric (no rotation, no lateral mode) on
    phi wedges with the production Gauss-Legendre 7 x 18 set (azimuthal period 20 deg):
    a 20-deg wedge (exact phi-wall maps) and narrow 3.7 and 7.3 deg wedges with
    vet_gd_wall_sym = true; the wedge-mean F_phi must vanish: fphi = |sum F3| / sum |F1|
    after 8 cycles.  The narrow wedge without the key must stop with the wedge FATAL.
    The same shell made theta-periodic (20 deg wedge) with vet_gd_theta_seam_sym = true: the
    seam face from both sides (the equatorial mirror), fth = |sum F_theta| / sum |F_r|;
    without the key it must stop with the theta-image FATAL.  fth_rot: the theta-periodic
    30 deg wedge with GL 11 x 12 (odd nmu: every ring pairs with its equatorial mirror)
    turning about z with the production rotation stack of vet_gd_rotate_dt (EMA ring,
    cross-fade), both seam keys, after 8 cycles; GL 10 x 12 staggered (even nmu, no
    mirror) must stop with the mirror FATAL (rot_fatal)."""
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
        c.run('shell_sp.athinput', 'w3.7_nosym', dict(o, **{
            'mesh/x3max': 0.064, 'rad_vet/vet_gd_wall_sym': 'false',
            'rad_vet/vet_gd_wall_interp': 'false'}))
        out['wedge_fatal'] = 0
    except RuntimeError:
        log = open(os.path.join(c.dir, 'w3.7_nosym', 'run.log')).read()
        out['wedge_fatal'] = int('vet_gd_wall_sym' in log and 'FATAL' in log)
    th = dict(o, **{'mesh/x3max': np.pi/9, 'mesh/ix2_bc': 'periodic',
                    'mesh/ox2_bc': 'periodic', 'rad_vet/vet_gd_wall_sym': 'true'})
    d = c.run('shell_sp.athinput', 'th20', dict(th, **{'rad_vet/vet_gd_theta_seam_sym': 'true'}))
    m = vett.last_bin(d, 'vet')
    out['fth'] = float(abs(m['vet_f2'].sum())/np.abs(m['vet_f1']).sum())
    out['warnings'] += vett.warnings(open(os.path.join(d, 'run.log')).read())
    try:
        c.run('shell_sp.athinput', 'th20_nosym', th)
        out['theta_fatal'] = 0
    except RuntimeError:
        log = open(os.path.join(c.dir, 'th20_nosym', 'run.log')).read()
        out['theta_fatal'] = int('vet_gd_theta_seam_sym' in log and 'FATAL' in log)
    r = dict(th, **{'mesh/x3max': np.pi/6, 'rad_vet/vet_gd_gl_nphi': 12,
                    'rad_vet/vet_gd_rotate_dt': 2.0, 'rad_vet/vet_gd_theta_seam_sym': 'true',
                    # (strict solves: with the budget keys the iterate is mirror-symmetric only
                    # to the time-error tolerance, and the turning set carries that seed on)
                    'rad_vet/clean_clos1_budget': 'false', 'rad_vet/clean_bud_maxnorm': 'false'})
    d = c.run('shell_sp.athinput', 'th30rot', dict(r, **{'rad_vet/vet_gd_gl_nmu': 11}))
    m = vett.last_bin(d, 'vet')
    out['fth_rot'] = float(abs(m['vet_f2'].sum())/np.abs(m['vet_f1']).sum())
    out['warnings'] += vett.warnings(open(os.path.join(d, 'run.log')).read())
    try:
        c.run('shell_sp.athinput', 'th30rot_even', dict(r, **{'rad_vet/vet_gd_gl_nmu': 10}))
        out['rot_fatal'] = 0
    except RuntimeError:
        log = open(os.path.join(c.dir, 'th30rot_even', 'run.log')).read()
        out['rot_fatal'] = int('not mirror-symmetric' in log and 'FATAL' in log)
    return out


# ---------------------------------------------------------------------------------------
# floors with a general EOS: the floor energy is paid from the radiation

def t_floor_pay(c):
    """the closed shell with hydro/eos = general (cgs units = code units) and a
    temperature floor just above the gas temperature, so the floor raises T in most cells
    every cycle; hydro/floor_pay = rad_column pays that energy from the radiation of the
    same column.  dE: total-energy drift from t = 0 (the initial conversion's floor energy
    is paid before the t = 0 output; dE0 = the t = 0 -> 1 step alone); paid / unpaid:
    floor energy over the run (floor_audit).  The ekin run: hydro/floor_from_ekin = true
    with a velocity ceiling below the shell's speeds (vceil 1.2e-2): the ceiling's kinetic
    energy becomes heat, phi momentum kept (L_z): ekin_dE, ekin_dM, ekin_dLz as
    conservation_sp.  worst_lines / worst_bitwise: the same run with
    hydro/floor_audit_worst prints its cells and leaves the history bitwise."""
    g = {'hydro/eos': 'general', 'units/length_cgs': 1.0, 'units/mass_cgs': 1.0,
         'units/time_cgs': 1.0}
    d = c.run('shell_sp.athinput', 'run', dict(g, **{
        'hydro/tfloor': 1.00001e-4, 'hydro/floor_pay': 'rad_column',
        'hydro/floor_pay_skip': 'true', 'hydro/dfloor_flux_limit': 'true',
        'hydro/floor_audit': 1}))
    u = os.path.join(d, 'sp.user.hst')
    e = hst_col(u, 'Eg') + hst_col(u, 'Er')
    log = open(os.path.join(d, 'run.log')).read()
    pay = [[float(x) for x in m] for m in
           re.findall(r'floor_audit cycle \d+ pay: paid (\S+) unpaid (\S+)', log)]
    out = dict(dE=rel_drift(e), dE0=float(abs(e[1] - e[0])/e[0]),
               paid=sum(p[0] for p in pay), unpaid=sum(p[1] for p in pay),
               warnings=vett.warnings(log))
    o = dict(g, **{'hydro/tfloor': 1.00001e-4, 'hydro/floor_pay': 'rad_column',
                   'hydro/floor_pay_skip': 'true', 'hydro/dfloor_flux_limit': 'true',
                   'hydro/floor_audit': 1, 'hydro/floor_audit_worst': 'true'})
    dw = c.run('shell_sp.athinput', 'worst', o)
    out['worst_lines'] = open(os.path.join(dw, 'run.log')).read().count('floor_audit_worst')
    out['worst_bitwise'] = all(open(os.path.join(d, f), 'rb').read() ==
                               open(os.path.join(dw, f), 'rb').read()
                               for f in ('sp.user.hst', 'sp.hydro.hst'))
    d = c.run('shell_sp.athinput', 'ekin', dict(g, **{
        'time/nlim': 16, 'hydro/floor_from_ekin': 'true', 'hydro/vceil': 1.2e-2}))
    out.update({'ekin_' + k: v for k, v in am_budget(d, 'sp').items() if k != 'dLz_gr'})
    out['warnings'] += vett.warnings(open(os.path.join(d, 'run.log')).read())
    return out


# the velocity ceiling in the metric: floor_from_ekin keeps u_phi (L_z), cuts the rest

CS_FRAME = {0: ((0, 1, 0), (0, 0, 1), (1, 0, 0)), 1: ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
            2: ((0, -1, 0), (0, 0, 1), (-1, 0, 0)), 3: ((0, 1, 0), (-1, 0, 0), (0, 0, 1)),
            4: ((1, 0, 0), (0, 0, 1), (0, -1, 0)), 5: ((0, 1, 0), (1, 0, 0), (0, 0, -1))}


def phys_velocity(d, cs):
    """the t = 0 hydro_w velocity of run dir d per cell as Cartesian-frame (V.r_hat,
    V.theta_hat, V.phi_hat): sp the orthonormal components; the cubed sphere the
    contravariant v^a on the panel's unit tangents (cubed_sphere::PanelTangents,
    cubed_sphere::PanelFrame), rebuilt as V = v^1 r_hat + v^2 e_xi + v^3 e_eta"""
    out = []
    for p, (ext, a) in enumerate(vett.read_blocks(vett.bins(d, 'hydro_w')[0])):
        v = np.stack([a['velx'], a['vely'], a['velz']])
        if not cs:
            out.append(v.reshape(3, -1))
            continue
        n3, n2, _ = a['velx'].shape
        c = [ext[2*m] + (np.arange(n) + 0.5)*(ext[2*m + 1] - ext[2*m])/n
             for m, n in ((1, n2), (2, n3))]
        x, y = np.meshgrid(np.tan(0.25*np.pi*c[0]), np.tan(0.25*np.pi*c[1]))
        fa, fb, fn = [np.array(f, float)[:, None, None] for f in CS_FRAME[p]]
        dl = np.sqrt(1 + x*x + y*y)
        q = (fa*x + fb*y + fn)/dl
        e1 = (fa*dl - q*x)/np.sqrt(1 + y*y)
        e2 = (fb*dl - q*y)/np.sqrt(1 + x*x)
        e1, e2 = e1/np.linalg.norm(e1, axis=0), e2/np.linalg.norm(e2, axis=0)
        ph = np.stack([-q[1], q[0], 0*q[0]])/np.hypot(q[0], q[1])
        th = np.cross(ph, q, axis=0)
        V = v[0]*q[..., None] + v[1]*e1[..., None] + v[2]*e2[..., None]
        out.append(np.stack([(V*b[..., None]).sum(0) for b in (q, th, ph)]).reshape(3, -1))
    return np.concatenate(out, 1)


def t_ceiling_metric(c):
    """the rotating shells (gas_rot 0.3 about z, gas_vr 0.2) with hydro/eos = general,
    floor_from_ekin and vceil = 0.5, below |v| in most cells, against the same run without
    the ceiling, cell by cell at t = 0 (after the initial inversion): u_phi kept (L_z),
    the rest scaled along itself (theta part stays 0), |V| (cs: the metric norm of the
    contravariant v^a) = max(vceil, |v_phi|) where the ceiling acts and unchanged
    elsewhere.  A ceiling on Sum v_a^2 of mixed components, or a split along the wrong
    (co- / contravariant) phi component, fails it on the cubed sphere"""
    vc = 0.5
    o = {'hydro/eos': 'general', 'units/length_cgs': 1.0, 'units/mass_cgs': 1.0,
         'units/time_cgs': 1.0, 'hydro/floor_from_ekin': 'true', 'output2/variable':
         'hydro_w', 'time/nlim': 1, 'problem/gas_rho': 1.0, 'problem/gas_temp': 0.5,
         'problem/e_out': 1.0, 'problem/r0': 2.0, 'problem/width': 0.3,
         'problem/gas_rot': 0.3, 'problem/gas_vr': 0.2, 'problem/rho_lat_amp': 0.0,
         'problem/lat_amp': 0.0, 'mesh/x1max': 3.0}
    out = {}
    for m, deck in (('sp', 'shell_sp.athinput'), ('cs', 'shell_cs.athinput')):
        a = phys_velocity(c.run(deck, m + '_free', dict(o, **{'hydro/vceil': 0.0})),
                          m == 'cs')
        b = phys_velocity(c.run(deck, m + '_ceil', dict(o, **{'hydro/vceil': vc})),
                          m == 'cs')
        va, vb = np.linalg.norm(a, axis=0), np.linalg.norm(b, axis=0)
        cap = np.where(va > vc, np.maximum(vc, np.abs(a[2])), va)
        rad = np.abs(a[0]) > 1e-8*va.max()   # (the radial part's scale, where it has one)
        sr = b[0][rad]/a[0][rad]
        out.update({m + '_dvphi': float(np.abs(b[2] - a[2]).max()/va.max()),
                    m + '_dcap': float(np.abs(vb - cap).max()/va.max()),
                    m + '_vtheta': float(np.abs(b[1]).max()/va.max()),
                    m + '_scale_ok': int(np.all((sr >= -1e-13) & (sr <= 1 + 1e-13))),
                    m + '_capped': int((va > vc).sum())})
    return out


# ---------------------------------------------------------------------------------------
# a non-finite iterate is rejected, never floored

def t_inf_reject(c):
    """the closed shell with rad_vet/clean_dbg_inf = 4: the first stage-1 write-back of
    cycle 4 gives one cell a non-finite E iterate.  implicit_pos_floor must not raise it
    (no booking of a non-finite raise: nonfinite = implicit_pos_floor creation warnings
    that print inf or nan); the write-back's check rejects the stage, which is redone
    (rejected >= 1), and the run stays energy-exact (dE as conservation_sp)."""
    d = c.run('shell_sp.athinput', 'run', {
        'time/nlim': 8, 'rad_vet/implicit_pos_floor': 'true',
        'rad_vet/implicit_pos_floor_cons': 'true', 'rad_vet/clean_dbg_inf': 4})
    log = open(os.path.join(d, 'run.log')).read()
    cw = [ln for ln in log.splitlines() if 'CLEAN WARNING' in ln]
    out = am_budget(d, 'sp')
    fl = [ln for ln in cw if 'implicit_pos_floor created' in ln]
    return dict(dE=out['dE'], nonfinite=sum(bool(re.search(r'\b(inf|nan)\b', ln,
                                                           re.I)) for ln in fl),
                rejected=sum(('stage fallback' in ln) or ('dt_retry' in ln) for ln in cw))


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
# IDORT literature tests (Ma et al. 2025; Jiang 2021): beams, shadow, emitting cylinder and
# radiative shocks against exact or semi-analytic solutions

def beam_list(d):
    """[(E_b, angle deg, x2 lo, x2 hi)] of the implicit_beam keys of run d, with the in-plane
    angle the code snapped each beam to (its log line), which the exact solution uses."""
    p = vett.parse_deck_text(open(os.path.join(d, 'in.athinput')).read())['rad_vet']
    ang = dict(re.findall(r'implicit_beam(\d+): .*in-plane ang=(\S+) deg',
                          open(os.path.join(d, 'run.log')).read()))
    return [(float(p['implicit_beam%d_e' % b]), float(ang[str(b)]),
             float(p['implicit_beam%d_y0' % b]), float(p['implicit_beam%d_y1' % b]))
            for b in range(int(p['implicit_beam_n']))]


def mesh_of(d):
    return vett.parse_deck_text(open(os.path.join(d, 'in.athinput')).read())['mesh']


def beams_exact(d, ns=16):
    """cell averages (ns x ns sub-samples) of straight beams entering at x1min, x2 periodic"""
    p, (x1, x2, _) = mesh_of(d), vett.cell_centres(d)
    x1a, x2a = float(p['x1min']), float(p['x2min'])
    L = float(p['x2max']) - x2a
    o = (np.arange(ns) + 0.5)/ns - 0.5
    xs, ys = (x1[:, None] + o*(x1[1] - x1[0])).ravel(), (x2[:, None] + o*(x2[1] - x2[0])).ravel()
    E = np.zeros((len(x2), len(x1)))
    for e, ang, y0, y1 in beam_list(d):
        y = (ys[:, None] - (xs[None, :] - x1a)*np.tan(np.radians(ang)) - y0) % L
        E += e*(y < y1 - y0).reshape(len(x2), ns, len(x1), ns).mean(axis=(1, 3))
    return E


def l1(a, b):
    return float(np.abs(a - b).sum()/np.abs(b).sum())


def t_idort_beams(c):
    """two crossing beams in vacuum: E against the exact straight beams (L1 over the box;
    lum_top: the energy of the third row below the x1max face relative to the exact);
    flux_err: the largest |x1 energy flux through a row (sum of vet_phi1) / injected - 1|"""
    d = c.run('idort_beams.athinput', 'run')
    m, Ex = vett.last_bin(d, 'vet'), beams_exact(d)
    E, x2 = m['vet_e'][0], vett.cell_centres(d)[1]
    cl = float(vett.parse_deck_text(open(os.path.join(d, 'in.athinput')).read())
               ['rad_vet']['c_light'])
    fin = sum(e*np.cos(np.radians(a))*(y1 - y0) for e, a, y0, y1 in beam_list(d))
    fx = m['vet_phi1'][0].sum(axis=0)*(x2[1] - x2[0])/cl
    return dict(L1=l1(E, Ex), lum_top=float(E[:, -3].sum()/Ex[:, -3].sum()),
                flux_err=float(np.abs(fx/fin - 1.0).max()))


def t_idort_beams_mb(c):
    """the beams on 4 MeshBlocks against 1 MeshBlock: along x2 in tile mode (identical),
    along x1 in block mode (tile mode needs one MeshBlock along x1; the result depends on
    the MeshBlock size at the solver tolerance)."""
    t = {'mesh/layout_invariant_mode': 'tile', 'mesh/layout_tile': 16}
    bl = {'mesh/layout_invariant_mode': 'block'}
    a, b = c.run('idort_beams.athinput', 'tile1', t), c.run('idort_beams.athinput', 'block1', bl)
    x2 = c.run('idort_beams.athinput', 'tile4x2', dict(t, **{'meshblock/nx2': 16}))
    x1 = c.run('idort_beams.athinput', 'block4x1', dict(bl, **{'meshblock/nx1': 64}))
    return dict(maxdiff_tile_x2=vett.compare_dirs(a, x2)[1],
                maxdiff_block_x1=vett.compare_dirs(b, x1)[1])


def shadow_exact(d, nsub=4096):
    """E of the beams attenuated along straight rays back to x1min through the ellipse
    rho = rho0 (1 + (ratio - 1)/(1 + exp(delta (r^2 - 1)))), r^2 = (x/ax)^2 + (y/ay)^2,
    rho kappa = rho^2 T^-3.5 = rho^5.5 (T = 1/rho); point values at the cell centres."""
    p, (x1, x2, _) = mesh_of(d), vett.cell_centres(d)
    q = vett.parse_deck_text(open(os.path.join(d, 'in.athinput')).read())['problem']
    ax, ay, dl, r0, rt = (float(q['shadow_' + k]) for k in ('ax', 'ay', 'delta', 'rho0', 'ratio'))
    x1a, x2a = float(p['x1min']), float(p['x2min'])
    L, lx = float(p['x2max']) - x2a, float(p['x1max']) - x1a
    E = np.zeros((len(x2), len(x1)))
    for e, ang, _, _ in beam_list(d):
        t = np.tan(np.radians(ang))
        for i, xc in enumerate(x1):
            n = int(np.ceil(nsub*(xc - x1a)/lx)) + 8
            s = x1a + (np.arange(n) + 0.5)*(xc - x1a)/n
            y = (x2[:, None] - (xc - s)*t - x2a) % L + x2a
            z = np.clip(dl*((s/ax)**2 + (y/ay)**2 - 1.0), -700.0, 700.0)
            rho = r0*(1.0 + (rt - 1.0)/(1.0 + np.exp(z)))
            E[:, i] += e*np.exp(-(rho**5.5).sum(axis=1)*(xc - x1a)/n/np.cos(np.radians(ang)))
    return E


def t_idort_shadow(c, deck='idort_shadow.athinput'):
    """two beams past an opaque ellipse: E against the exact attenuated beams (L1)"""
    d = c.run(deck, 'run')
    return dict(L1=l1(vett.last_bin(d, 'vet')['vet_e'][0], shadow_exact(d)))


def t_idort_shadow_gl(c):
    return t_idort_shadow(c, 'idort_shadow_gl.athinput')


def t_idort_cylinder(c):
    """homogeneous emitting cylinder: E inside and outside, the energy flux Phi_x (vet_phi1)
    and the face momentum F_x (vet_f1, FH_out) outside, along y = 0 against the exact
    transport solution (ref/idort_cylinder_exact.npz, x > 0)."""
    d = c.run('idort_cylinder.athinput', 'run')
    m, (x1, x2, _) = vett.last_bin(d, 'vet'), vett.cell_centres(d)
    cl = float(vett.parse_deck_text(open(os.path.join(d, 'in.athinput')).read())
               ['rad_vet']['c_light'])
    ex = np.load(os.path.join(HERE, 'ref', 'idort_cylinder_exact.npz'))
    j, s = np.argsort(np.abs(x2))[:2], x1 > 0          # the two rows next to y = 0
    E, F, H = (m[v][0][j].mean(axis=0)[s] for v in ('vet_e', 'vet_phi1', 'vet_f1'))
    i, o = ex['x'] < 1.0, ex['x'] > 1.05
    return dict(E_in=l1(E[i], ex['E'][i]), E_out=l1(E[o], ex['E'][o]),
                F_out=l1(F[o]/cl, ex['F'][o]), FH_out=l1(H[o]/cl, ex['F'][o]))


def t_idort_radshock(c, m, deck='idort_radshock_m%s.athinput'):
    """Lowrie & Edwards radiative shock started from its semi-analytic profile: rho, v,
    T_gas, T_rad at the end against the profile (no shift), L1 / mean |reference|."""
    ref = os.path.join(HERE, 'ref', 'idort_radshock_m%s.txt' % m)
    d = c.run(deck % m, 'run', {'problem/vet_shock_ref': ref})
    p = vett.parse_deck_text(open(os.path.join(d, 'in.athinput')).read())
    tu, ar = float(p['problem']['vet_shock_tunit']), float(p['rad_vet']['arad'])
    h, r = vett.last_bin(d, 'hydro_w'), vett.last_bin(d, 'vet')
    rho = h['dens'].mean(axis=(0, 1))
    q = dict(rho=rho, v=h['velx'].mean(axis=(0, 1)),
             Tg=(float(p['hydro']['gamma']) - 1.0)*h['eint'].mean(axis=(0, 1))/rho/tu,
             Tr=(np.maximum(r['vet_e'].mean(axis=(0, 1)), 0.0)/ar)**0.25/tu)
    a = np.loadtxt(ref)
    x = vett.cell_centres(d)[0] - float(p['problem']['vet_shock_xs'])
    out = {}
    for n, (k, v) in enumerate(q.items()):
        e = np.interp(x, a[:, 0], a[:, n + 1])
        out['L1_' + k] = float(np.abs(v - e).mean()/np.abs(e).mean())
    return out


def t_idort_radshock_m2(c):
    return t_idort_radshock(c, '2')


def t_idort_radshock_m5(c):
    return t_idort_radshock(c, '5')


def t_idort_thermal_eq2(c):
    """hot gas cooling into radiation at large steps (implicit_cfl 1.6): T(t) and the
    total energy against the exact ODE"""
    errT, errE = relax_errors(c.run('idort_thermal_eq2.athinput', 'run'))
    return dict(errT=errT, dEtot=errE)


# ---------------------------------------------------------------------------------------
# classic radiation tests

def warn(d):
    return vett.warnings(open(os.path.join(d, 'run.log')).read())


def x1prof(a, k):
    """x2-x3 mean profile of dump variable k."""
    return a[k].mean(axis=(0, 1))


def t_diffusion(c):
    """Gaussian pulse in a static pure scatterer: the variance grows at 2 D t,
    D = c/(3 rho kappa_s), from tau_cell = 1 (transport, for information) to 1e5
    (asymptotic preserving); errD = growth/(2 D t) - 1, L1 against the exact Gaussian
    (diffusion equation), the same pulse isotropic in 2-D (64^2) and a top-hat at
    tau_cell = 1e4: E must stay within [background, background + amplitude]."""
    out = {'errD': [], 'L1': [], 'warnings': 0}
    cases = [('tau1', 1.0, {}), ('tau10', 10.0, {}), ('tau1e3', 1e3, {}),
             ('tau1e5', 1e5, {}), ('2d', 1e3, {'mesh/nx2': 64, 'meshblock/nx2': 64,
                                                'mesh/nx1': 64, 'meshblock/nx1': 64,
                                                'mesh/layout_tile': 16,
                                                'problem/pulse_1d': 'false'}),
             ('tophat', 1e4, {'problem/vet_test': 'tophat'})]
    for tag, tc, o in cases:
        n = o.get('mesh/nx1', 128)
        D = 1.0/(3.0*tc*n)
        t = 5e-4/(2.0*D)                    # 2 D t = 5e-4: the variance grows by 40 %
        d = c.run('rad_diffusion.athinput', tag, dict({
            'rad_vet/kappa_s': tc*n, 'time/tlim': t, 'time/dt_fixed': t/100}, **o))
        out['warnings'] += warn(d)
        x, y, _ = vett.cell_centres(d)
        E = [m['vet_e'][0] if tag == '2d' else x1prof(m, 'vet_e')
             for m in (vett.read_bin(vett.bins(d, 'vet')[0]), vett.last_bin(d, 'vet'))]
        E = [e - 1e-3 for e in E]           # minus the background pulse_bg
        if tag == 'tophat':
            out['th_min'], out['th_max'] = float(E[1].min()), float(E[1].max() - 1.0)
            continue
        r2, nd = ((x[None, :] - 0.5)**2 + y[:, None]**2, 2) if tag == '2d' else \
            ((x - 0.5)**2, 1)
        var = [(e*r2).sum()/e.sum()/nd for e in E]   # per dimension
        s2 = 0.05**2 + 4.0*D*t
        ex = (0.05**2/s2)**(0.5*nd)*np.exp(-r2/s2)
        err = [float((var[1] - var[0])/(2.0*D*t) - 1.0),
               float(np.abs(E[1] - ex).sum()/ex.sum())]
        if tag in ('tau1', '2d'):
            out['errD_' + tag], out['L1_' + tag] = err
        else:
            out['errD'].append(err[0])
            out['L1'].append(err[1])
    return out


def t_radshock(c, m=2):
    """the Lowrie & Edwards shocks of idort_radshock_m* (same reference) with the production
    closure (clean_closure default, not eddington): L1 of rho, v, T_gas, T_rad; E > 0"""
    out = t_idort_radshock(c, m, 'rad_radshock_m%s.athinput')
    d = os.path.join(c.dir, 'run')
    out.update(Emin=float(vett.last_bin(d, 'vet')['vet_e'].min()), warnings=warn(d))
    return out


def t_radshock_m5(c):
    return t_radshock(c, 5)


def t_marshak(c):
    """non-equilibrium Marshak wave (constant c_v): no analytic solution, so robustness +
    a regression value: E and T_gas positive and bounded by the bath (maximum principle;
    Emax, Tmax relative to E_b, T_b), E non-increasing into the slab up to round-off
    (Erise: largest increase / E_b), and the front x_f (T_gas = T_b / 2) at t = 200."""
    d = c.run('rad_marshak.athinput', 'run')
    m, h = vett.last_bin(d, 'vet'), vett.last_bin(d, 'hydro_w')
    p = vett.parse_deck_text(open(os.path.join(d, 'in.athinput')).read())
    eb = float(p['rad_vet']['implicit_ebath_x1min'])
    tb = (eb/float(p['rad_vet']['arad']))**0.25
    E, T = x1prof(m, 'vet_e'), (2.0/3.0)*x1prof(h, 'eint')/x1prof(h, 'dens')
    x = vett.cell_centres(d)[0]
    i = int(np.argmax(T < 0.5*tb))
    xf = x[i - 1] + (x[i] - x[i - 1])*(T[i - 1] - 0.5*tb)/(T[i - 1] - T[i])
    return dict(Emin=float(m['vet_e'].min()), Tmin=float(T.min()), Emax=float(E.max()/eb),
                Tmax=float(T.max()/tb), Erise=float(max(np.diff(E).max(), 0.0)/eb),
                x_front=float(xf), warnings=warn(d))


def t_advect_uniform(c):
    """uniform thick medium moving at 0.01 c, at the coupling's fixed point: the largest
    relative change of T_gas and E over the run (dumps every 2e-8 s)."""
    d = c.run('rad_advect_uniform.athinput', 'run')
    T0, dT, dE = None, 0.0, 0.0
    for fm, fh in zip(vett.bins(d, 'vet'), vett.bins(d, 'hydro_w')):
        m, h = vett.read_bin(fm), vett.read_bin(fh)
        T = h['eint']/h['dens']
        if T0 is None:
            T0, E0 = T.mean(), m['vet_e'].mean()
        dT = max(dT, float(np.abs(T/T0 - 1.0).max()))
        dE = max(dE, float(np.abs(m['vet_e']/E0 - 1.0).max()))
    return dict(dT=dT, dE=dE, warnings=warn(d))


def t_advect_pulse(c):
    """advecting radiation pulse (dynamic diffusion): the pulse moving at v against the same
    pulse at rest, shifted by v t (periodic): L1 = sum |T_v(x + v t) - T_0(x)| /
    sum |T_0 - T_bg| of the gas temperature at each dump (0.6e-6 s apart)."""
    a = c.run('rad_advect_pulse.athinput', 'moving')
    b = c.run('rad_advect_pulse.athinput', 'rest', {'problem/pulse_v': 0.0})
    p = vett.parse_deck_text(open(os.path.join(a, 'in.athinput')).read())
    v, Tbg = float(p['problem']['pulse_v']), float(p['problem']['pulse_t0'])
    x = vett.cell_centres(a)[0]
    L = float(p['mesh']['x1max']) - float(p['mesh']['x1min'])
    out = {'L1': [], 'warnings': warn(a) + warn(b)}
    for fa, fb in list(zip(vett.bins(a, 'hydro_w'), vett.bins(b, 'hydro_w')))[1:]:
        ha, hb = vett.read_bin(fa), vett.read_bin(fb)
        Ta, Tb = [x1prof(h, 'eint')/x1prof(h, 'dens')*2.0/3.0 for h in (ha, hb)]
        Ts = np.interp(x + v*ha['time'], x, Ta, period=L)
        out['L1'].append(float(np.abs(Ts - Tb).sum()/np.abs(Tb - Tbg).sum()))
    return out


def shadow_cells(x, y, x0=0.5, ax=0.1, ay=0.06, ymax=0.12, n=241):
    """(umbra, lit) masks: cells outside the clump (D >= 1.5) every straight sight line of
    which (in the x-y plane) to the illuminated face x = 0, |y0| <= ymax, crosses the clump
    core D < 1 (umbra); cells no sight line of which comes within D < 2 (lit)."""
    X, Y = np.meshgrid(x, y)
    y0 = np.linspace(-ymax, ymax, n)
    s = np.linspace(0.0, 1.0, 400)
    dmin = np.full(X.shape + (n,), np.inf)
    for k in range(n):   # D along the segment (0, y0) -> (x, y)
        px = X[..., None]*s
        py = y0[k] + (Y[..., None] - y0[k])*s
        dmin[..., k] = np.sqrt(((px - x0)/ax)**2 + (py/ay)**2).min(axis=-1)
    own = np.sqrt(((X - x0)/ax)**2 + (Y/ay)**2)
    return (dmin < 1.0).all(axis=-1) & (own >= 1.5), (dmin >= 2.0).all(axis=-1)


def t_shadow(c):
    """shadow of an opaque clump against the same run without the clump (shadow_ratio =
    1) after 3 light-crossing times (steady): umbra = sum E / sum E_noclump over the
    geometric umbra behind the clump (the radiation left in its shadow); lit_err = max |E / E_noclump - 1|
    where no sight line to the illuminated face comes near the clump (D >= 2)."""
    a = c.run('rad_shadow.athinput', 'clump')
    b = c.run('rad_shadow.athinput', 'noclump', {'problem/shadow_ratio': 1.0})
    x, y, _ = vett.cell_centres(a)
    um, lit = shadow_cells(x, y)
    ea, eb = (vett.last_bin(z, 'vet')['vet_e'][0] for z in (a, b))
    return dict(umbra=float(ea[um].sum()/eb[um].sum()),
                lit_err=float(np.abs(ea[lit]/eb[lit] - 1.0).max()),
                n_umbra=int(um.sum()), n_lit=int(lit.sum()), warnings=warn(a) + warn(b))


# ---------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------
# the star tests: the production star deck's settings on a small spherical-polar wedge
# (a binary built with -D PROBLEM=star_vet, --star-exe; the opacity and profile tables of
# inputs/stars/star3d_data are linked into the run directory)

STAR_DATA = os.path.abspath(os.path.join(HERE, '..', '..', 'inputs', 'stars', 'star3d_data'))
# columns of the histories compared with ref/star_smoke.*.hst (Lz and the momenta of the
# unrotating star are sums of cancelling noise and are left out)
STAR_COLS = {'hydro': ('mass', 'tot-E', '1-KE', '2-KE', '3-KE'),
             'user': ('time', 'dt', 'L_bot', 'L_mid', 'L_int', 'L_top', 'L_in', 'E_rad',
                      'e_gas', 'M_tot', 'Etot', 'KE_int', 'L_midE', 'L_intE')}


def star_run(c, sub, over):
    d = os.path.join(c.dir, sub)
    os.makedirs(d, exist_ok=True)
    os.symlink(STAR_DATA, os.path.join(d, 'star3d_data'))
    exe, c.exe = c.exe, c.star
    try:
        return c.run('star_small.athinput', sub, over)
    finally:
        c.exe = exe


def t_star_smoke(c):
    """the reduced star3d_sp deck (256 x 8 x 8, the decided production scheme keys and
    tables) for 20 cycles against the stored history: dev = largest |x - x_ref| / max|x_ref|
    over the columns STAR_COLS and all rows."""
    d = star_run(c, 'run', {})
    out = {'dev': 0.0, 'warnings': vett.warnings(open(os.path.join(d, 'run.log')).read())}
    for kind, cols in STAR_COLS.items():
        lab, a = vett.read_hst(os.path.join(d, 'bsg3d.%s.hst' % kind))
        rl, r = vett.read_hst(os.path.join(HERE, 'ref', 'star_smoke.%s.hst' % kind))
        if a.shape != r.shape:
            raise RuntimeError('%s history: %s rows x columns, reference %s' % (
                kind, a.shape, r.shape))
        for q in cols:
            x, y = a[:, lab.index(q)], r[:, rl.index(q)]
            v = float(np.abs(x - y).max()/np.abs(y).max())
            if v > out['dev']:
                out['dev'], out['worst'] = v, q
    return out


def t_wb_star(c):
    """the well-balance mode of the stars (wellbalance_dynamic, wb_x1, force_reference =
    wb_arad, fixed GL 7 x 18 rays) on the unseeded BSG column, 20 cycles.  A report, not a
    criterion: the BSG column is not in discrete equilibrium with its own opacities, so
    v/c_s does not stay at 1e-8 (see README); it passes when the run completes with no
    CLEAN WARNING.  v_in: largest |v|/c_s where r < 2.9e12 (the interior); v_all: over the
    mesh; lat_mid: largest lateral |v|/c_s at 2.44e12 < r < 2.95e12 (outside the
    merged-cell radius), ~1e-11 in a symmetric star."""
    d = star_run(c, 'run', {'problem/star_ic_profile_amp': 0.0})
    r = vett.cell_centres(d)[0]
    h = vett.last_bin(d, 'hydro_w')
    cs = np.sqrt((5.0/3.0)*(2.0/3.0)*h['eint']/h['dens'])
    vt = np.sqrt(h['vely']**2 + h['velz']**2)
    v = np.sqrt(h['velx']**2 + vt**2)/cs
    m = (r > 2.44e12) & (r < 2.95e12)
    return dict(v_in=float(v[:, :, r < 2.9e12].max()), v_all=float(v.max()),
                lat_mid=float((vt/cs)[:, :, m].max()),
                warnings=vett.warnings(open(os.path.join(d, 'run.log')).read()))


TESTS = [  # name, function, in --quick
    ('conservation_sp', t_conservation_sp, True),
    ('conservation_cs', t_conservation_cs, True),
    ('conservation_lz', t_conservation_lz, False),
    ('floorcap_spin', t_floorcap_spin, False),
    ('floor_dt', t_floor_dt, False),
    ('restart_sp', t_restart_sp, True),
    ('restart_cs', t_restart_cs, False),
    ('dt_retry', t_dt_retry, True),
    ('layout_tile', t_layout_tile, True),
    ('order_time', t_order_time, True),
    ('order_stiff', t_order_stiff, False),
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
    ('ceiling_metric', t_ceiling_metric, True),
    ('inf_reject', t_inf_reject, False),
    ('star_smoke', t_star_smoke, False),
    ('wb_star', t_wb_star, False),
    ('idort_beams', t_idort_beams, False),
    ('idort_beams_mb', t_idort_beams_mb, False),
    ('idort_shadow', t_idort_shadow, True),
    ('idort_shadow_gl', t_idort_shadow_gl, True),
    ('idort_cylinder', t_idort_cylinder, False),
    ('idort_radshock_m2', t_idort_radshock_m2, False),
    ('idort_radshock_m5', t_idort_radshock_m5, False),
    ('idort_thermal_eq2', t_idort_thermal_eq2, False),
    # classic radiation tests
    ('diffusion', t_diffusion, True),
    ('radshock', t_radshock, False),
    ('radshock_m5', t_radshock_m5, False),
    ('marshak', t_marshak, False),
    ('advect_uniform', t_advect_uniform, True),
    ('advect_pulse', t_advect_pulse, False),
    ('shadow', t_shadow, False),
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
    p.add_argument('--star-exe', default='', help='a binary built with -D PROBLEM=star_vet '
                   '(star_smoke, wb_star; without it they are skipped)')
    p.add_argument('--quick', action='store_true', help='the short subset (16 tests)')
    p.add_argument('--jobs', type=int, default=max(1, (os.cpu_count() or 4)//4),
                   help='tests run at the same time')
    p.add_argument('--threads', type=int, default=4, help='OpenMP threads per run')
    p.add_argument('--only', default='', help='comma-separated test names')
    p.add_argument('--workdir', default='vet_clean_runs', help='where the runs go')
    p.add_argument('--mpirun', default='mpirun', help='MPI launcher (multi-rank tests)')
    p.add_argument('--keep', action='store_true', help='keep the run directories of passes')
    p.add_argument('--record', default='', help='write the measured values to this JSON')
    p.add_argument('--keys', default='', help="extra overrides of every run, e.g. "
                   "'rad_vet/clean_tol=1e-9 hydro/rad_signal_speed_diff=true'")
    a = p.parse_args()
    sel = [t for t in TESTS if (not a.quick or t[2]) and
           (not a.only or t[0] in a.only.split(','))]
    if not a.star_exe:
        for n in ('star_smoke', 'wb_star'):
            if any(t[0] == n for t in sel):
                print('SKIP %-16s        needs --star-exe (a binary built with '
                      '-D PROBLEM=star_vet)' % n)
        sel = [t for t in sel if t[0] not in ('star_smoke', 'wb_star')]
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
