#!/usr/bin/env python3
"""tables_ext 1-D (480x8x8, seed off) A/B: new-table arm vs old-table arm, same binary/keys.
usage: ana_run1d.py X  (X = A or B) -> prints per dump: interior max |d rho|, |dT|, |dE| (new vs old,
r < 0.9 R_ph); per band (0.9-1, 1-1.1, 1.1-1.5, 1.5-3 R_ph) median T, E ratios new/old and the
T ranges; L(r) = 4 pi r^2 F_r at 0.9, 1, 1.1, 1.5, 2, 2.9 R_ph (/L*); R_ph(tau_R 2/3, each arm's own
table); cost: cycles, mean dt, wall, NEWTON-FALLBACK and FATAL counts; writes checks/run1d_X.png."""
import glob
import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.dont_write_bytecode = True
sys.path.insert(0, '/viper/ptmp2/jinma/lbv_1008/agcar/geos/scripts')
sys.path.insert(0, '/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/scripts')
import bin_convert_agcar as bc  # noqa: E402
import make_ic_mlt_star_ge as mk  # noqa: E402
import column_compare as cc  # noqa: E402

RS = 6.957e10
R = '/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/run1d'
X = sys.argv[1]
RPH, LST = {'A': (388.3007, 3.411709e39), 'B': (101.2962, 5.662027e39)}[X]
EOS = mk.TableEOS(R + '/A_new/eos_dump.txt' if os.path.exists(R + '/A_new/eos_dump.txt')
                  else R + '/smoke2A_new/eos_dump.txt')   # log T 3.4..6.6, log rho -21..-2


def prof(D, k):
    hb = sorted(glob.glob(D + '/bin/*.hydro_w.*.bin'))
    mb = sorted(glob.glob(D + '/bin/*.m1.*.bin'))
    gf = glob.glob(D + '/*.x1grid.txt')[0]
    with open(gf) as fh:
        is_ = int(fh.readline().split()[3])
    grid = np.array([[float(v) for v in ln.split()] for ln in open(gf)
                     if not ln.startswith('#') and len(ln.split()) == 4])
    f, g = bc.read_binary(hb[k]), bc.read_binary(mb[k])
    mp = lambda h, v: np.mean([h['mb_data'][v][m].mean(axis=(0, 1))  # noqa: E731
                               for m in range(h['n_mbs'])], axis=0)
    rho, eint = mp(f, 'dens'), mp(f, 'eint')
    E, F1 = mp(g, 'm1_e'), mp(g, 'm1_f1')
    r = grid[is_:is_ + len(rho), 2]
    dr = grid[is_:is_ + len(rho), 3]
    lr = np.log10(np.maximum(rho, 1e-21))
    # bisection in log T on the EOS dump (e(T) monotonic at fixed rho); Newton from a fixed
    # start jumped across the ionisation plateaus at the floor density
    lo, hi = np.full_like(lr, 3.4), np.full_like(lr, 6.6)
    le = np.log10(eint/rho)
    for _ in range(50):
        mid = 0.5*(lo + hi)
        up = EOS.le(mid, lr) < le
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    lT = 0.5*(lo + hi)
    res = np.abs(EOS.le(lT, lr) - le)
    if res.max() > 1e-3:
        print('   (T inversion: %d cells with |dlog e| > 1e-3, max %.2e: clamped at the EOS edge)'
              % ((res > 1e-3).sum(), res.max()))
    return dict(t=f['time'], r=r, dr=dr, x=r/RS/RPH, rho=rho, T=10**lT, E=E, L=4*np.pi*r**2*F1/LST)


def rph(p, tab):
    k = cc.look(tab, p['T'], p['rho'])
    tau = np.cumsum((k*p['rho']*p['dr'])[::-1])[::-1]
    j = np.where(tau >= 2/3)[0]
    return p['x'][j[-1]] if j.size else np.nan


def cost(D):
    log = open(D + '/run.log').read()
    cyc = re.findall(r'cycle=(\d+) time=([0-9.e+-]+)', log)
    wall = re.search(r'Elapsed \(wall clock\) time.*: (\S+)', open(D + '/time.txt').read())
    h = np.loadtxt(glob.glob(D + '/*.hydro.hst')[0])
    return (int(cyc[-1][0]) if cyc else -1, float(cyc[-1][1]) if cyc else np.nan,
            h[:, 1].mean(), h[:, 1].min(), wall.group(1) if wall else '?',
            log.count('NEWTON-FALLBACK'), log.count('FATAL'))


def main():
    NEW = sys.argv[2] if len(sys.argv) > 2 else 'new'
    Do, Dn = R + '/%s_old' % X, R + '/%s_%s' % (X, NEW)
    n = min(len(glob.glob(Do + '/bin/*.hydro_w.*.bin')), len(glob.glob(Dn + '/bin/*.hydro_w.*.bin')))
    bands = ((0.9, 1.0), (1.0, 1.1), (1.1, 1.5), (1.5, 3.0))
    xs = (0.5, 0.9, 1.0, 1.1, 1.5, 2.0, 2.9)
    for nm, D in (('old', Do), ('new', Dn)):
        c = cost(D)
        print('%s %s: last cycle %d t %.4e s | dt mean %.1f min %.1f s | wall %s | '
              'NEWTON-FALLBACK lines %d | FATAL %d' % (X, nm, *c))
    fig, ax = plt.subplots(2, 2, figsize=(12, 7))
    for k in range(n):
        po, pn = prof(Do, k), prof(Dn, k)
        assert np.allclose(po['r'], pn['r'])
        m = po['x'] < 0.9
        print('t %.3e s | interior r<0.9 R_ph max |drho/rho| %.2e |dT/T| %.2e |dE/E| %.2e | '
              'R_ph old %.4f new %.4f (own tables)' % (
                  po['t'], np.abs(pn['rho'][m]/po['rho'][m] - 1).max(),
                  np.abs(pn['T'][m]/po['T'][m] - 1).max(),
                  np.abs(pn['E'][m]/po['E'][m] - 1).max(),
                  rph(po, cc.T_['old'][0]), rph(pn, cc.T_['new'][0])))
        for a, b in bands:
            mm = (po['x'] >= a) & (po['x'] < b)
            print('   %.1f-%.1f R_ph: T old %6.0f-%6.0f new %6.0f-%6.0f K, median T new/old %.3f '
                  '(min %.3f max %.3f), E new/old median %.3f (min %.3f max %.3f), '
                  'rho new/old median %.3f' % (
                      a, b, po['T'][mm].min(), po['T'][mm].max(), pn['T'][mm].min(),
                      pn['T'][mm].max(), np.median(pn['T'][mm]/po['T'][mm]),
                      (pn['T'][mm]/po['T'][mm]).min(), (pn['T'][mm]/po['T'][mm]).max(),
                      np.median(pn['E'][mm]/po['E'][mm]), (pn['E'][mm]/po['E'][mm]).min(),
                      (pn['E'][mm]/po['E'][mm]).max(), np.median(pn['rho'][mm]/po['rho'][mm])))
        idx = [np.argmin(abs(po['x'] - v)) for v in xs]
        print('   L/L* at r/R_ph %s: old %s | new %s' % (
            ' '.join('%.2f' % v for v in xs), ' '.join('%.3f' % po['L'][i] for i in idx),
            ' '.join('%.3f' % pn['L'][i] for i in idx)))
        if k == n - 1:
            for j, (key, lab) in enumerate((('T', 'T [K]'), ('E', 'E [erg/cm3]'),
                                            ('L', 'L_rad/L*'), ('rho', 'rho'))):
                a = ax.flat[j]
                for p, nm, ls in ((po, 'old', '-'), (pn, 'new', '--')):
                    (a.semilogx if key == 'L' else a.loglog)(p['x'], p[key], ls,
                                                             label='%s t=%.2e' % (nm, p['t']))
                a.axvspan(0.9, 3.0, color='0.92')
                a.set_xlabel('r/R_ph')
                a.set_ylabel(lab)
                a.legend(fontsize=7)
    fig.suptitle('AG Car %s 1-D (480x8x8, seed off): old vs extended tables' % X)
    fig.tight_layout()
    out = '/viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/checks/run1d_%s%s.png' % (X, '' if len(sys.argv) < 3 else '_' + sys.argv[2])
    fig.savefig(out, dpi=110)
    print('wrote', out)


if __name__ == '__main__':
    main()
