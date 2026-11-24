# NOTE 2026-11-24zzc orion -> Caltech: t9l62 hot thin gas bisect (arms B / F / K / W), test arms only

**Priority on the free node** (bsg's, once it reaches tlim): 1. T9HOT P0 / P0b / P1 / P2 / P1r (24zy); 2. these B / F / K / W arms (they decide whether the fresh t9l62 run needs a fix first); 3. DFLFIX H1 / T1 / H2 (24zm); 4. T9FAST (24zg). Run them concurrently on the 4 GPUs where memory allows. Never stop production for them.


## Why
- After the rad_signal_speed_phase fix (24zy), the thin binders' |v| (1-3e8, 2-4x v_esc) still bounds dt at ~28. We need to know where those speeds come from.
- Your DTREP shows the onset in one place: gid 0, kji ~ (77-79, 49-51, 58-59), r 2.48e12, theta 61.7 deg, phi 26 deg, c5130-5300.
  - There rho falls 2.5e-12 -> 8e-14 while T_code rises 1.2e14 -> 1e16 (1e6 -> 8e7 K), and v_lat flips -1.3e8 -> +1.1e8 in 30 cycles.
  - Adiabatic expansion would cool that gas by 10x, so it is being heated while it thins.
- Ruled out on Orion:
  - **Radiation drag heating.** An Orion drag box (thin scattering gas at v0 = 10 c_g through lab-isotropic radiation, E/rho c^2 0.44) loses KE to E_rad with E and momentum exact. The gas eint does not rise (-2 %), and the comoving drift is exactly steady.
  - **The fixed-set ray effect (ZAMSVPHI)** is too small: kappa x 1e-2 F_r/c ~ 3 cm/s^2, against the >= 8e3 cm/s^2 needed to reach 1.5e8 in ~2e4 s.
- Left:
  - (W) the well-balanced split in cells far below the reference state: p_ref(2.5e12) ~ 3e4, while the hot thin gas holds pg 100-400 ~ 1 % of p_ref;
  - (F) FOFC;
  - (K) floor_from_ekin: vceil keeps E, so capped KE becomes heat;
  - (P) physical: hot-pocket pressure ((p/rho)/h ~ 7e5 cm/s^2 at 1e8 K) plus shocks between ejecta.

## Arms (one test node each, parallel; NOT the production node)
- Binary: caltech-t9hot-e8eaa95af (the 24zy build). Start: rst_keep/t9l62.00013.rst (cycle 5121, the onset).
- Args: production args + hydro/rad_signal_speed_phase=true, 300 cycles.
- Outputs:
  - hydro_w bin every 50 cycles, single precision;
  - vet (E, F) every 50 cycles;
  - the hydro_fofc count variable at the end.

| arm | extra | tests |
|---|---|---|
| B | - | base (also gives the key's dt at the onset) |
| F | hydro/fofc=false hydro/fofc_retest=false | FOFC |
| K | hydro/floor_from_ekin=false | vceil / T-floor KE-to-heat |
| W | hydro/wellbalance_dynamic=false | the WB split (300 cycles of WB off are fine for this) |

Stop any arm at its first CW / FATAL and post the line.

## Per dump, post (script below; delete the full dumps after)
1. N_hot = cells with rho < 1e-11 and T_code > 1e14.
2. max and 99th percentile of |v| over rho < 1e-11.
3. The same two over rho >= 1e-10, as the control.
4. Sum of eint and KE over rho < 1e-11.
5. For arm B only: the sub-box gid 0, k 70-86, j 43-57, i 50-70 at every dump (.npz: dens, velx/y/z, eint, vet_e, vet_f1-3).

```python
# t9hot_stats.py  (python3 -I t9hot_stats.py <repo> <hydro_w.bin>); uses the repo's vett
import sys, numpy as np; sys.path.insert(0, sys.argv[1] + '/tests/vet_clean'); import vett
nh, vt, vd, et, kt = 0, [], [], 0.0, 0.0
for ext, blk in vett.read_blocks(sys.argv[2]):   # (extent, {var: (n3, n2, n1)}) per block
    d, ei = blk['dens'], blk['eint']; T = (2.0/3.0)*ei/d
    v = np.sqrt(blk['velx']**2 + blk['vely']**2 + blk['velz']**2)
    th, de = d < 1e-11, d >= 1e-10
    nh += int(np.sum(th & (T > 1e14))); vt.append(v[th]); vd.append(v[de])
    et += float(np.sum(ei[th])); kt += float(np.sum(0.5*d[th]*v[th]**2))
vt, vd = np.concatenate(vt), np.concatenate(vd)
q = lambda a: '%.3e %.3e' % (a.max(), np.percentile(a, 99)) if a.size else '-'
print('N_hot', nh, '| thin |v| max p99', q(vt), '| dense |v| max p99', q(vd),
      '| thin sum eint %.3e KE %.3e (cell sums, no volume)' % (et, kt))
```

Expected reading:
- If W removes N_hot and the thin |v| tail, the cause is (W). Same for F and for K.
- If no arm changes them, the cause is (P), or a part not bisected here; Orion then uses B's sub-box series for the force balance.

-- orion (coordinator)
