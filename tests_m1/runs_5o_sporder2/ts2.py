#!/usr/bin/env python3
"""T-S2 free-streaming point source on the wedge (runs_5b: m1 closure, kappa = 0, F_in = 1
at r_in, Marshak q = 1 at r_out), steady.  Exact: E = F_in r_in^2/(c r^2) (f = 1).

usage: ts2.py run|eval <tag> <arm> [n ...]
arm: old (the input as is) | lin (implicit_marshak_face = linear); suffix _s = the he4
stretched radial grid.  Error: L1_V and Linf of E/E_exact - 1 at the centroids
(r^2 E const there).
"""
import os
import subprocess
import sys

import numpy as np

import solib as so

SRC = '/viper/ptmp2/jinma/sph2_0924/inp/sp_sph_fs.athinput'
ADD = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts', 'addkeys.py')


def mk(tag, arm, n):
    d = os.path.join(so.W, 'cpu', 'ts2', tag, '%s_%d' % (arm, n))
    os.makedirs(d, exist_ok=True)
    inp = os.path.join(d, 'in.athinput')
    rad = []
    if arm.startswith('lin'):
        rad.append('implicit_marshak_face=linear')
    subprocess.check_call(['python3', ADD, SRC, inp + '.1', 'rad_m1'] + rad)
    mesh = []
    if arm.endswith('_s'):
        mesh = ['use_grid_stretch_r_poly=true'] + ['%s=%s' % kv for kv in so.POLY.items()]
    subprocess.check_call(['python3', ADD, inp + '.1', inp, 'mesh'] + mesh)
    return d, inp


def main():
    act, tag, arm = sys.argv[1:4]
    ns = [int(x) for x in sys.argv[4:]] or [32, 64, 128, 256]
    if act == 'run':
        b = os.path.join(so.W, 'bin', 'athena_%s_none_cpu' % tag)
        for n in ns:
            d, inp = mk(tag, arm, n)
            args = ['mesh/nx1=%d' % n, 'meshblock/nx1=%d' % n]
            print(arm, n, 'rc', so.run(b, inp, d, 1, args), flush=True)
        return
    e1, ei = [], []
    for n in ns:
        d = mk(tag, arm, n)[0]
        fs = sorted(f for f in os.listdir(os.path.join(d, 'tab')) if f.endswith('.tab'))
        t = so.tab(os.path.join(d, 'tab', fs[-1]))
        rf = so.rfaces(1.0, 3.0, n, arm.endswith('_s'))
        rl, rr = rf[:-1], rf[1:]
        q = rl/rr
        rc = 0.25*(q*q + 1.0)/((q*q + q + 1.0)/3.0)*(rr + rl)
        # the discrete steady state keeps r^2 E constant at the centroids
        ea = 1.0/(rc*rc)
        err = np.abs(t['m1_e']/ea - 1.0)
        w = rr**3 - rl**3
        e1.append(np.sum(w*err)/np.sum(w))
        ei.append(err.max())
    print('%s %s: L1 %s  Linf %s' % (tag, arm, ' '.join('%.2e' % x for x in e1),
                                     ' '.join('%.2e' % x for x in ei)))


if __name__ == '__main__':
    main()
