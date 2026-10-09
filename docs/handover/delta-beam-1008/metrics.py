"""rad-beam-1008 AG Car A metrics (TASK-2026-10-08-delta-beam-batchN).
usage: python3 metrics.py <rundir> [<refdir>]
From the LAST bin set (m1, m1_fs, hydro_w of the same output number) of <rundir>:
  shell means at r/R_ph = 0.9 ... 2.9: L_M1 = 4 pi r^2 <F_r>/L0, L_fs = 4 pi r^2 c <H_r>/L0,
  f_M1 = <F_r>/(c<E>), f_fs = <H>/<J>, <E>/<J>, inward fraction of cells (F_r < 0), S/E;
  lateral rms/mean of E; and with <refdir> the max |E/E_ref - 1| and |rho/rho_ref - 1| over
  r < 0.9 R_ph at the same output number (if present).
Cost: the end-of-run <rad_m1> lines (Picard, inner iterations, fallbacks) and the cycle /
cpu-time summary of run.log."""
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..',
                                'vis', 'python'))
import bin_convert as bc  # noqa: E402

C, RPH, L0 = 2.99792458e10, 2.70140777e13, 3.41170859e39
RS = (0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.1, 1.2, 1.5, 2.0, 2.5, 2.9)


def last(d, tag):
    fs = sorted(glob.glob(d + '/bin/*.%s.*.bin' % tag))
    return fs[-1] if fs else None


def main():
    d = sys.argv[1]
    ref = sys.argv[2] if len(sys.argv) > 2 else None
    fm = last(d, 'm1')
    if fm is None:
        print('no m1 bin in', d)
        return
    num = fm.split('.')[-2]
    m = bc.read_binary_as_athdf(fm)
    ffs = fm.replace('.m1.', '.m1_fs.')
    s = bc.read_binary_as_athdf(ffs) if os.path.exists(ffs) else None
    r = np.asarray(m['x1v'], dtype=float)
    E, F = m['m1_e'], m['m1_f1']
    print('# %s  bin %s  t %.4e  cycle %d' % (d, num, m['Time'], m['NumCycles']))
    print(' r/Rph   L_M1   L_fs  f_M1   f_fs   E/J   S/E  inward  Erms/E')
    for x in RS:
        i = int(np.argmin(np.abs(r - x*RPH)))
        e = E[:, :, i].mean()
        fr = F[:, :, i].mean()
        line = '%5.3f %6.3f' % (r[i]/RPH, 4*np.pi*r[i]**2*fr/L0)
        if s is not None and s['m1_jfs'][:, :, i].mean() > 0:
            j = s['m1_jfs'][:, :, i].mean()
            h = s['m1_hfs'][:, :, i].mean()
            line += ' %6.3f %5.3f %6.3f %5.3f %5.3f' % (4*np.pi*r[i]**2*C*h/L0, fr/(C*e), h/j,
                                                        e/j, s['m1_sfs'][:, :, i].mean()/e)
        else:
            line += '    -   %5.3f    -      -     -  ' % (fr/(C*e))
        line += ' %6.3f %6.3f' % (np.mean(F[:, :, i] < 0), E[:, :, i].std()/e)
        print(line)
    if ref is not None:
        fr_ = glob.glob(ref + '/bin/*.m1.%s.bin' % num)
        fw = glob.glob(ref + '/bin/*.hydro_w.%s.bin' % num)
        if fr_ and fw:
            mr = bc.read_binary_as_athdf(fr_[0])
            w0 = bc.read_binary_as_athdf(fw[0])
            w1 = bc.read_binary_as_athdf(fm.replace('.m1.', '.hydro_w.'))
            k = r < 0.9*RPH
            de = np.abs(E[:, :, k]/mr['m1_e'][:, :, k] - 1).max()
            dr = np.abs(w1['dens'][:, :, k]/w0['dens'][:, :, k] - 1).max()
            print('# vs %s (bin %s, t %.4e): r < 0.9 R_ph max|dE/E| %.3e max|drho/rho| %.3e'
                  % (ref, num, mr['Time'], de, dr))
        else:
            print('# vs %s: no bin %s' % (ref, num))
    lg = os.path.join(d, 'run.log')
    if os.path.exists(lg):
        for ln in open(lg):
            if (ln.startswith('<rad_m1> implicit transport') or 'inner iterations' in ln
                    or ln.startswith('<rad_m1> offdiag') or re.search(r'cpu time used|cycles =|'
                                                                     r'zone-cycles/cpu', ln)):
                print('# ' + ln.rstrip()[:200])


main()
