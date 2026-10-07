#!/usr/bin/env python3
"""Stage-2 envelope analysis (ry_per_accretor, inner = envelope).
usage: ana_env.py <rundir> [spin]
 - hst: Menv, Jenv drift, wall / R_acc / stream budgets (ryper.user.hst)
 - bin: max |v - v_init| (rotating frame) in r < R_acc, in r < R_acc - 0.2 (resolved
   part) and over gas with rho > 1e-3 rho_ph; shell-mean inertial Omega/Omega_orb.
Code units: Rsun, km/s; Omega_orb = 7.375171809554155, R_acc 4.06."""
import glob
import sys
import numpy as np
sys.path.insert(0, '/viper/ptmp2/jinma/wt_accretor/vis/python')
import bin_convert as bc  # noqa: E402

OM = 7.375171809554155
RACC = 4.06
d = sys.argv[1]
spin = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0


def hst(fn):
    rows, lab = [], None
    for line in open(fn):
        if line.startswith('# Athena'):
            rows = []
        elif line.startswith('#'):
            lab = [w.split('=')[1] for w in line[1:].split() if '=' in w]
        else:
            rows.append([float(x) for x in line.split()])
    return lab, np.array(rows)


def load(fn):
    fd = bc.read_binary(fn)
    bl = []
    for b in range(fd['n_mbs']):
        md = fd['mb_data']
        bl.append((fd['mb_x3v'][b][0], fd['mb_x1v'][b], fd['mb_x3v'][b], md['dens'][b],
                   md['velx'][b], md['velz'][b]))
    bl.sort(key=lambda t: t[0])
    r = bl[0][1]
    ph = np.concatenate([b[2] for b in bl])
    D, VR, VP = [np.concatenate([b[q] for b in bl], axis=0) for q in (3, 4, 5)]
    return fd.get('time', np.nan), r, ph, D, VR, VP


lab, a = hst(glob.glob(d + '/*.user.hst')[0])
c = {n: i for i, n in enumerate(lab)}
t = a[:, c['time']]
print(f"{d}: hst t {t[0]:.4f}..{t[-1]:.5f} (= {t[-1]/0.85194:.3f} orbits)")
for n in ('Menv', 'Jenv', 'Mdom', 'Jdom'):
    print(f"  {n}: {a[0, c[n]]:.6e} -> {a[-1, c[n]]:.6e}  rel {a[-1, c[n]]/a[0, c[n]]-1:+.3e}")
for n in ('Min', 'Mout', 'Mwal', 'Jwal', 'Jstr', 'Jout', 'MR', 'JR', 'Jin'):
    print(f"  cum {n}: {a[-1, c[n]]:.5e}", end='')
print()
bins = sorted(glob.glob(d + '/bin/*.bin'))
if bins:
    t0, r, ph, D0, VR0, VP0 = load(bins[0])
    t1, r, ph, D1, VR1, VP1 = load(bins[-1])
    dv = np.sqrt((VR1 - VR0)**2 + (VP1 - VP0)**2)
    R = r[None, None, :]*np.ones_like(D1)
    rhoph = D0[:, :, np.argmin(abs(r - RACC + 0.006))].mean()
    for name, msk in (('r<R_acc', R < RACC), ('r<R_acc-0.2', R < RACC - 0.2),
                      ('rho>1e-3 rho(R_acc)', D1 > 1e-3*rhoph)):
        i = np.unravel_index(np.argmax(np.where(msk, dv, -1)), dv.shape)
        print(f"  bin t {t0:.4f}->{t1:.4f} max|v-v0| {name}: {dv[i]:.4g} km/s at r "
              f"{r[i[2]]:.4f} phi {ph[i[0]]:.3f}; max|drho/rho| "
              f"{np.max(np.where(msk, abs(D1/D0-1), 0)):.3e}")
    # shell-mean inertial angular velocity (density weighted)
    for rr in (RACC - 1.0, RACC - 0.5, RACC - 0.2, RACC - 0.1, RACC - 0.05, RACC - 0.02,
               RACC - 0.006):
        ii = np.argmin(abs(r - rr))
        w = D1[:, :, ii]
        om = ((VP1[:, :, ii]/r[ii] + OM)*w).sum()/w.sum()/OM
        print(f"  r {r[ii]:.4f}: <Omega>/Omega_orb {om:.4f} (spin {spin})")
