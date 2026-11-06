# TASK 2026-11-06 Caltech -> viper (USER): early kinetic energy in your real Athena++ 1-D BSG runs

Numbers only, from the existing Raven runs (`/raven/ptmp/jinma/athpp_audit_1015`, `athpp_fix_1016`); no new runs needed unless nothing suitable exists.

## Why
Linhao's 3-D Fig. 2 starts at E_k ~ 1e40 erg and is ~7e40 at 0.5 d, quiet until ~3 d. Our AthenaK arms with the Athena++ velocity mimic
(`rad_m1/doppler_vel = vstar_n`) get a coherent radial kick right after t = 0, which he does not show (KE inside r_int, erg, whole wedge pi/3 sr):

| t | `none` (no v*) | vstar_n CFL 0.3 | vstar_n CFL 0.35 |
|---|---|---|---|
| ~1000 s | 1.25e40 | 6.7e42 | 1.9e43 |
| ~2000 s | 2.4e40 | 9.7e41 | 7.9e42 |
| ~3000 s | 3.2e40 | 2.9e41 | 2.0e42 |

We need to know whether the REAL v* ordering in Athena++ produces such a kick at all.

## Please report, per arm (at least his unmodified code/BC arm, plus any arm closest to his 3-D setup)
1. Arm name, nx1, nmu, dt, and whether hydro + radiation are his unmodified code.
2. Column-integrated gas KE vs time at t = 0, then every ~500 s to 3 h, then at 6 h, 12 h, 1 d, 2 d, 3 d (from .hst or dumps), with the integration range and any volume scaling so it can be compared with a pi/3 sr wedge.
3. Max |v_r| and where (r) at those times, and the wall luminosity (inner-face flux) at the same times.
4. Does KE have a sharp peak in the first ~1000 s, or does it grow smoothly from the seed?

The user has confirmed you do not have his `density_unit`; no need to look for it.

Reply as a NOTE on this branch.
