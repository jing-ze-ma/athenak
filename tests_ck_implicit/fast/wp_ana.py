"""ck-fast well-posed gates: tables for A1, A2, C, D at dt = 20 s.

usage: wp_ana.py <root> [<root2> ...]   (arms are the sub-directories of <root>/<test>;
a later root's arm of the same name overrides an earlier one)
A1: max/mean |T - T**|/T** over the ck cells after 200 calls, T** = the t4x arm (the
    state every stable arm settles on after the restart step; see wellposed/README.md),
    and against the t4 arm; drift per call over the last 10 calls.
A2: max/mean relative T error after 100 calls against the wellposed t4x dt = 4 s
    reference (/viper/ptmp2/jinma/wellposed_0923/A/a2_ref), and against the t4 arm.
D:  eos_h2 = true, mu0 = 1: max relative T error vs the t4x arm after 10 and 100 calls.
C:  energy budget per call (the hook's column budget, sum V (e - e*) - dt sum V S, over
    |dt sum V S|; and over the column energy) and the final residual, from A1 and A2.
Columns 'passes' / 'nonconv' are the per-call pass counts and non-converged calls.
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rec  # noqa: E402

ROOTS = sys.argv[1:]
A2REF = '/viper/ptmp2/jinma/wellposed_0923/A/a2_ref/rec.txt'


def arms(test):
    d = {}
    for r in ROOTS:
        for p in sorted(glob.glob(os.path.join(r, test, '*', 'rec.txt'))):
            d[os.path.basename(os.path.dirname(p))] = p
    return d


def rel(a, b, ic):
    e = np.abs(a[ic:] - b[ic:])/b[ic:]
    return e.max(), e.mean()


def stats(r):
    s = r['s']
    return s['passes'].mean(), int(s['passes'].max()), int(s['nonconv'].sum())


def A1():
    a = arms('A1')
    if 't4x' not in a:
        print('A1: no t4x arm')
        return
    ref = rec.read(a['t4x'])
    ic = ref['icut']
    tr = rec.last(ref)
    t4 = rec.last(rec.read(a['t4'])) if 't4' in a else None
    print('A1 (200 calls at dt 20 s from the A0b state), ck cells i >= %d' % ic)
    print('%-5s %10s %10s %10s %10s %10s %7s %4s %7s' % (
        'arm', 'max vs T**', 'mean', 'max vs t4', 'mean', 'drift/cl', 'passes', 'maxp',
        'nonconv'))
    for k, f in a.items():
        r = rec.read(f)
        t = rec.last(r)
        ks = sorted(r['T'])
        dr = rel(r['T'][ks[-1]], r['T'][ks[-2]], ic)[0]/(ks[-1] - ks[-2])
        m1 = rel(t, tr, ic)
        m2 = rel(t, t4, ic) if t4 is not None else (np.nan, np.nan)
        p = stats(r)
        print('%-5s %10.3e %10.3e %10.3e %10.3e %10.3e %7.2f %4d %7d' % (
            k, m1[0], m1[1], m2[0], m2[1], dr, p[0], p[1], p[2]))


def A2():
    a = arms('A2')
    ref = rec.read(A2REF)
    ic = ref['icut']
    tr = rec.last(ref)
    t4 = rec.last(rec.read(a['t4'])) if 't4' in a else None
    print('A2 (100 calls at dt 20 s from the IC) vs t4x dt 4 s (500 calls), i >= %d' % ic)
    print('%-5s %10s %10s %10s %10s %7s %4s %7s' % (
        'arm', 'max vs ref', 'mean', 'max vs t4', 'mean', 'passes', 'maxp', 'nonconv'))
    for k, f in a.items():
        r = rec.read(f)
        t = rec.last(r)
        m1 = rel(t, tr, ic)
        m2 = rel(t, t4, ic) if t4 is not None else (np.nan, np.nan)
        p = stats(r)
        print('%-5s %10.3e %10.3e %10.3e %10.3e %7.2f %4d %7d' % (
            k, m1[0], m1[1], m2[0], m2[1], p[0], p[1], p[2]))


def D():
    a = arms('D')
    if 't4x' not in a:
        print('D: no t4x arm')
        return
    ref = rec.read(a['t4x'])
    ic = ref['icut']
    print('D (eos_h2 = true, mu0 = 1, dt 20 s) vs the t4x arm, i >= %d' % ic)
    print('%-5s %11s %11s %11s %7s %4s %7s %10s' % (
        'arm', 'maxrel@10', 'maxrel@100', 'mean@100', 'passes', 'maxp', 'nonconv',
        'Trange'))
    for k, f in a.items():
        r = rec.read(f)
        e10 = rel(r['T'][10], ref['T'][10], ic)[0]
        e100 = rel(r['T'][100], ref['T'][100], ic)
        p = stats(r)
        tl = r['T'][100][ic:]
        print('%-5s %11.3e %11.3e %11.3e %7.2f %4d %7d %4.0f-%4.0f' % (
            k, e10, e100[0], e100[1], p[0], p[1], p[2], tl.min(), tl.max()))


def C():
    print('C: energy budget per call, max over calls')
    print('%-9s %9s %9s %11s %11s %11s' % ('run', 'res(max)', 'resfin', '|bud|/|dtS|',
                                           'budmax', '|bud|/E'))
    for t in ('A1', 'A2'):
        for k, f in arms(t).items():
            s = rec.read(f)['s']
            be = np.abs(s['dE0'] - s['dtS0'])/s['sumVe0']
            print('%-9s %9.2e %9.2e %11.3e %11.3e %11.3e' % (
                t + '_' + k, s['res'].max(), s['resfinal'].max(),
                np.abs(s['bud0']).max(), s['budmax'].max(), be.max()))


for f in (A1, A2, D, C):
    f()
    print()
