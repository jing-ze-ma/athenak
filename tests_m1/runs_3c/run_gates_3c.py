#!/usr/bin/env python3
"""run_gates_3c.py -- the gate driver of milestone 3c (implicit_flux = blend).

Runs every case of the 3c gate table, analyses it with the tests_m1 gate scripts and
DELETES its dumps immediately (the viper inode quota).  Numbers land in
tests_m1/runs_3c/results_<gate>.json and the profiles for the results page in
tests_m1/plots/impl3c_*.json.

    python3 run_gates_3c.py G1 G2 G3 G4 G5 G6 G7 G9 REG
    python3 run_gates_3c.py G8              # needs build_cpu_box

Every run goes into a scratch directory that is removed as soon as it is measured, so
the repository never holds a dump.
"""

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
ROOT = os.path.dirname(TESTS)
X = os.path.join(ROOT, 'build_cpu_m1', 'src', 'athena')
XB = os.path.join(ROOT, 'build_cpu_box', 'src', 'athena')
INP = os.path.join(ROOT, 'inputs', 'tests')
SCRATCH = os.environ.get('M1_SCRATCH', '/tmp/m1_3c')
PLOTS = os.path.join(TESTS, 'plots')

PICARD = re.compile(r'solves=(\S+) Picard iterations mean=(\S+) max=(\S+) '
                    r'NON-CONVERGED=(\S+)')


def run(binary, athinput, outdir, extra, timeout=3600):
    """run one case; return (stdout, Picard statistics dict)."""
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    cmd = [binary, '-i', athinput, '-d', outdir] + list(extra)
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if p.returncode != 0:
        sys.stderr.write(p.stdout[-3000:] + p.stderr[-3000:])
        raise RuntimeError('run failed: ' + ' '.join(cmd))
    m = PICARD.search(p.stdout)
    st = {}
    if m:
        st = {'solves': float(m.group(1)), 'mean': float(m.group(2)),
              'max': float(m.group(3)), 'nfail': float(m.group(4))}
    return p.stdout, st


def gate(script, args):
    """run a tests_m1 gate script and return its stdout."""
    cmd = [sys.executable, os.path.join(TESTS, script)] + list(args)
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.stdout + p.stderr


def tab(path_glob):
    """last column block of a .tab file as a list of rows of floats."""
    f = sorted(glob.glob(path_glob))[-1]
    rows = []
    for ln in open(f):
        if ln.startswith('#'):
            continue
        rows.append([float(v) for v in ln.split()])
    return rows


def fluxopts(name):
    """the <rad_m1> overrides that select one entry of the 3c table."""
    if name == 'central':
        return ['rad_m1/implicit_flux=central']
    if name == 'berthon':
        return ['rad_m1/implicit_flux=berthon']
    parts = name.split('-')
    kind = {'W1': 'tau', 'W2': 'f', 'W3': 'tau_f'}[parts[0]]
    o = ['rad_m1/implicit_flux=blend', 'rad_m1/implicit_blend=' + kind]
    for p in parts[1:]:
        if p == 'mean':
            o.append('rad_m1/implicit_blend_fmode=mean')
        elif p == 'dis':
            o.append('rad_m1/implicit_blend_mode=dissipation')
        else:
            k, v = p.split('=')
            o.append('rad_m1/implicit_blend_%s=%s' % (k, v))
    return o


IMPL = ['rad_m1/transport=implicit_x1', 'time/cfl_number=1e6']
EXTRA = []


# ---------------------------------------------------------------------------- G1
def g1(forms, cfls, out):
    res = {}
    for f in forms:
        for c in cfls:
            d = os.path.join(SCRATCH, 'g1')
            _, st = run(X, os.path.join(INP, 'rad_m1_atmosphere.athinput'), d,
                        IMPL + EXTRA + fluxopts(f)
                        + ['rad_m1/implicit_cfl=%g' % c, 'time/tlim=2000',
                           'output1/dt=2000', 'output2/dt=2000'])
            js = os.path.join(PLOTS, 'impl3c_atm_%s.json' % f) if c == cfls[-1] else None
            a = ['--quiet', '--label', f]
            if js:
                a += ['--json', js]
            o = gate('t9_atmosphere.py',
                     sorted(glob.glob(d + '/tab/*.m1.*.tab'))[-1:] + a)
            res['%s@%g' % (f, c)] = {'out': o.strip(), 'picard': st}
            print(f, c, o.strip().replace('\n', ' | '))
            shutil.rmtree(d)
    json.dump(res, open(out, 'w'), indent=1)


# ---------------------------------------------------------------------------- G2
def g2(forms, cfls, recons, out):
    res = {}
    for f in forms:
        for rc in recons:
            for c in cfls:
                d = os.path.join(SCRATCH, 'g2')
                _, st = run(X, os.path.join(INP, 'rad_m1_pulse1d.athinput'), d,
                            IMPL + EXTRA + fluxopts(f)
                            + ['rad_m1/implicit_cfl=%g' % c,
                               'rad_m1/implicit_recon=' + rc])
                fs = sorted(glob.glob(d + '/tab/*.m1.*.tab'))
                e0 = [r[3] for r in tab(fs[0])]
                e1 = [r[3] for r in tab(fs[-1])]
                amp = (max(e1) - min(e1))/(max(e0) - min(e0))
                key = '%s/%s@%g' % (f, rc, c)
                res[key] = {'amp': amp, 'picard': st}
                if c == 0.4:
                    res[key]['profile'] = [[r[0], r[3]] for r in tab(fs[-1])]
                print(key, amp, st)
                shutil.rmtree(d)
    json.dump(res, open(out, 'w'), indent=1)


# ---------------------------------------------------------------------------- G3
MARSH = ['rad_m1/implicit_bc_x1min=marshak', 'rad_m1/implicit_ebath_x1min=1.0e-10',
         'rad_m1/implicit_bc_x1max=marshak']
T6 = ['--table', os.path.join(TESTS, 'runs_3a', 't6_ref_sn.txt'), '--marshak',
      '--a-rad', '1e30', '--cv', '1.5', '--rho', '1', '--centre', '0.0', '--quiet']


def g3(forms, cfls, recons, nxs, out):
    res = {}
    for f in forms:
        for rc in recons:
            for nx in nxs:
                for c in cfls:
                    d = os.path.join(SCRATCH, 'g3')
                    _, st = run(X, os.path.join(INP, 'rad_m1_marshak.athinput'), d,
                                IMPL + EXTRA + MARSH + fluxopts(f)
                                + ['rad_m1/implicit_cfl=%g' % c,
                                   'rad_m1/implicit_recon=' + rc,
                                   'mesh/nx1=%d' % nx, 'meshblock/nx1=%d' % nx])
                    hy = sorted(glob.glob(d + '/tab/*.hydro_w.*.tab'))[-1]
                    o = gate('t6_marshak.py',
                             [sorted(glob.glob(d + '/tab/*.m1.*.tab'))[-1],
                              '--hydro', hy] + T6)
                    key = '%s/%s/n%d@%g' % (f, rc, nx, c)
                    res[key] = {'out': o.strip(), 'picard': st}
                    print(key, o.strip().replace('\n', ' | '), st)
                    shutil.rmtree(d)
    json.dump(res, open(out, 'w'), indent=1)


# ---------------------------------------------------------------------------- G4
JUMP = ['rad_m1/implicit_bc_x1min=efix', 'rad_m1/implicit_bc_x1max=efix']
T3B = ['--c', '1', '--kappa', '0.64', '--ratio', '1000', '--flux', '1.0e-3',
       '--e-left', '1.0', '--quiet']


def g4(forms, out):
    res = {}
    for f in forms:
        for c in (1.0, 100.0):
            d = os.path.join(SCRATCH, 'g4')
            _, st = run(X, os.path.join(INP, 'rad_m1_jump.athinput'), d,
                        IMPL + EXTRA + JUMP + fluxopts(f)
                        + ['rad_m1/implicit_cfl=%g' % c])
            o = gate('t3b_jump.py',
                     [sorted(glob.glob(d + '/tab/*.m1.*.tab'))[-1]] + T3B)
            res['%s@%g' % (f, c)] = {'out': o.strip(), 'picard': st}
            print(f, c, o.strip().replace('\n', ' | '))
            shutil.rmtree(d)
    json.dump(res, open(out, 'w'), indent=1)


# ---------------------------------------------------------------------------- G5
def g5(forms, cfls, kappas, out):
    res = {}
    for f in forms:
        for kap in kappas:
            for c in cfls:
                d = os.path.join(SCRATCH, 'g5')
                _, st = run(X, os.path.join(INP, 'rad_m1_thick_pulse.athinput'), d,
                            IMPL + EXTRA + fluxopts(f)
                            + ['rad_m1/implicit_cfl=%g' % c,
                               'rad_m1/kappa_s=%g' % kap])
                o = gate('t3_pulse.py',
                         sorted(glob.glob(d + '/bin/*.bin'))
                         + ['--c', '1', '--rho', '1', '--kappa', '%g' % kap,
                            '--subtract-min', '--quiet'])
                res['%s/k%g@%g' % (f, kap, c)] = {'out': o.strip(), 'picard': st}
                print(f, kap, c, o.strip().replace('\n', ' | '))
                shutil.rmtree(d)
    json.dump(res, open(out, 'w'), indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gates', nargs='*')
    ap.add_argument('--forms', default='central,berthon,W1,W2,W3')
    ap.add_argument('--cfls', default='')
    ap.add_argument('--recons', default='dc')
    ap.add_argument('--tag', default='')
    ap.add_argument('--extra', default='')
    a = ap.parse_args()
    if a.extra:
        EXTRA.extend(a.extra.split(','))
    forms = a.forms.split(',')
    recons = a.recons.split(',')
    os.makedirs(SCRATCH, exist_ok=True)
    os.makedirs(PLOTS, exist_ok=True)
    for g in a.gates:
        o = os.path.join(HERE, 'results_%s%s.json' % (g, a.tag))
        if g == 'G1':
            g1(forms, [float(x) for x in (a.cfls or '1,100,1e4').split(',')], o)
        elif g == 'G2':
            g2(forms, [float(x) for x in (a.cfls or '0.4,1,10').split(',')], recons, o)
        elif g == 'G3':
            g3(forms, [float(x) for x in (a.cfls or '1,10,100').split(',')], recons,
               [128, 64], o)
        elif g == 'G4':
            g4(forms, o)
        elif g == 'G5':
            g5(forms, [float(x) for x in (a.cfls or '1,100,1e4').split(',')],
               [1280.0, 128000.0, 128000000.0], o)
        else:
            raise SystemExit('unknown gate ' + g)


if __name__ == '__main__':
    main()
