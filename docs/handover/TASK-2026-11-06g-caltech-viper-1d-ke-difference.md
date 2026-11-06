# TASK 2026-11-06 Caltech -> viper (USER): what makes af1_cflrad 13x higher in KE than his1d?

## Why
From `cflrad_1106/cflrad_summary_interim.txt`, at cfl_rad 0.025, KE in the pi/3 wedge at 0.5 d:
- `his1d` 1.17e40
- `af1_cflrad` 1.58e41 (13.5x higher; also 5.43e40 vs 5.15e39 at 1 d)

Our AthenaK 1-D at the same step (dt 4.4 s, no FOFC, his floors, 256 grid, unseeded) follows `af1_cflrad`: 7.8e40 at 0.1 d, 1.56e41 at 0.3 d. KE barely depends on dt in our code (dt 171 s gives 7.4e40 and 2.0e41). So whatever separates your two arms may be what separates us from his setup.

**It is probably not the tables.** We compared his original `BSG/log_kappa_ross.txt` and `log_kappa_planck.txt` with our uniform-log-T padded tables along the IC at 20-80 Rsun, bilinear in logs: maximum difference 6e-5 (Rosseland), 2e-4 (Planck). Please check that your padded table agrees with his the same way, since yours may have been built differently from ours.

## Please do
Starting from `his1d` (cfl_rad 0.025, 384 uniform), change **one item at a time** toward `af1_cflrad`, and run each to 0.5-1 d:
1. units: old density_unit and T_unit vs his (`density_unit 6.15245e-5`, `T_unit 1.07887e6`); also GM and molecular_weight if they differ
2. nzeta 4 vs 5
3. floors: dfloor, pfloor, tfloor
4. tables: padded vs original
5. the IC conversion: profile, values at the bottom, and the seed if it differs between the arms (threshold and random stream)
6. anything else that differs in the two decks or pgen; please diff them and list everything

If one item doesn't explain it, say which pair does, and go on only as far as is cheap.

## Please report
- Table per arm: KE at 0.1, 0.25, 0.5, 1 d (KE of the same integration range), T/IC of cells 0-1, wall flux, dE.
- Which item carries the 13x, and through what: for example where KE sits (r) and whether it is the coherent infall front. At 0.5 d in your two arms, give max|v_r| and its r.
- Whether the item is something physical (such as a different value of g or L, or the IC) or numerical (such as a tolerance or floor that depends on units).
- If it's the units: which quantity scales with them, for example error_limit, nlimit, the Jacobi stop test, or a floor.

Short NOTE on this branch, please; don't delay the 3-D `prod_seeded` for it.
