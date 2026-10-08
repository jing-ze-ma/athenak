"""SPBLEND2: final-dump deviation of arm runs against a reference arm (same restart, same cycle).
Per band of r/R_ph: max and rms |E/E_ref - 1|, |rho/rho_ref - 1|, |eint/eint_ref - 1|, and max |dv1|/c_ref-scale
(km/s); plus L_cell(r)/L0 (shell mean of m1_f1) at XS for both.
usage: [CASE=B] python3 cmp.py <ref_rundir> <rundir> [<rundir> ...]"""
import glob
import os
import sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../vis/python'))
import bin_convert as bc  # noqa: E402
CASE = os.environ.get('CASE', 'A')   # CASE=B for AG Car B
RPH = {'A': 2.70140777e13, 'B': 101.2962*6.957e10}[CASE]
L0 = {'A': 3.41170859e39, 'B': 10**6.17*3.828e33}[CASE]
BANDS = ((0.0, 0.9), (0.9, 1.0), (1.0, 1.2), (1.2, 3.1))
XS = (1.0, 1.02, 1.05, 1.5, 2.0, 2.9)


def last(d, v):
    f = sorted(glob.glob(d + '/bin/*.%s.*.bin' % v))
    return bc.read_binary_as_athdf(f[-1]) if f else None


ref = sys.argv[1]
mr, wr = last(ref, 'm1'), last(ref, 'hydro_w')
r = mr['x1v'].astype(float)
print('# ref %s t=%.6e cyc %d' % (ref, mr['Time'], mr['NumCycles']))
print('# bands r/R_ph', BANDS, ': max|dE/E| rms | max|drho/rho| | max|de/e| | max|dv1| km/s ; then L_cell at', XS)


def lrow(m):
    f = m['m1_f1'].astype(float).mean(axis=(0, 1))
    return ' '.join('%.3f' % (4*np.pi*r[i]**2*f[i]/L0)
                    for i in [int(np.argmin(np.abs(r - x*RPH))) for x in XS])


print('%-10s %s' % ('ref', lrow(mr)))
for d in sys.argv[2:]:
    m, w = last(d, 'm1'), last(d, 'hydro_w')
    if m is None or m['NumCycles'] != mr['NumCycles']:
        print('%-10s no matching dump' % d.split('/')[-1])
        continue
    s = '%-10s' % d.split('/')[-1]
    for a, b in BANDS:
        k = (r >= a*RPH) & (r < b*RPH)
        de = np.abs(m['m1_e'][:, :, k].astype(float)/mr['m1_e'][:, :, k] - 1)
        dr = np.abs(w['dens'][:, :, k].astype(float)/wr['dens'][:, :, k] - 1)
        dei = np.abs(w['eint'][:, :, k].astype(float)/wr['eint'][:, :, k] - 1)
        dv = np.abs(w['velx'][:, :, k].astype(float) - wr['velx'][:, :, k]).max()/1e5
        s += ' | %.1e %.1e %.1e %.1e %.1e' % (de.max(), np.sqrt((de**2).mean()), dr.max(), dei.max(), dv)
    print(s + ' | ' + lrow(m))
