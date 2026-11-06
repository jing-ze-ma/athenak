# NOTE 2026-11-06 Orion -> Caltech: first 3-D Athena++ data, t = 0-0.17 d (reply to TASK-2026-11-06m)

**This is an early push so you can check the format and have Q1-Q3; the requested times beyond t = 0 don't exist yet.** The next push is at ~1 d, about 4 h away. Files are in `docs/handover/athenapp-bsg/orion3d/`: A-E are done with all columns, F is skipped. Radiation fields are in the dumps, so nothing is derived indirectly. Normalization: `README_normalization.txt`. Checks: `verify.txt`, `seed_check.txt`. `extra/` is a test at dump 13 (0.174 d), not a requested time. `summary.txt` has dumps 0 and 13 only for now; `bsg.hst` is the current history.

## Questions
1. **Seed: yes.** rms(drho/rho) at t = 0 is 1.446, 1.437, 1.441, 1.440 % in the first four shells, 1.443 % mean over the seeded shells (i = 0-36, r 20.16-31.41; theory 1.4434 %). It is **0 at r = 41 and 70**, because rho is below the threshold there. r = 70 is at the density floor.
2. **Non-spherical KE at r < 26 exceeds KE_sym at the first dump after t = 0 (0.0134 d), against ~2.5 d in yours.** It can't be bracketed earlier. KE_ns(r<26) / KE_sym(whole domain): 18 at 0.0134 d, 7.8 at 0.174 d. Top v_r modes at r = 21 are k ~ 54-63 (Nyquist 64), so this is grid-scale motion from the seed, and it may not be the same thing as your event.
3. **Job 227350:** 8 nodes, 896 ranks, 01:12:48 elapsed, cycle=3707 time=3.0957944873060201e+00. Disk 56G. Rate so far ~4.3 star-d per wall-d, and faster lately (Jacobi mean 1, recent 0.98 s per cycle); expect it to fall when iterations rise.

## Things you should know
- **Job plan (set before your task, USER):** job 227360 is queued at 17 nodes, 1904 ranks. When it starts, it stops the 8-node run and restarts from its final restart. This changes rank count only; the deck is unchanged. The stop and restart path passed a small test, but the takeover itself hasn't happened yet.
- **Disk limit:** a 2 TB guard on `output/` stops the chain near **10.5 d** at 184 GB per star-day, before his ~18 d event, unless the USER frees space.
- **Definitions:** dE is wedge here (x12 for viper's). L_in, L_top are cell-centred, because dumps hold no face fluxes. N^2 uses `wall.py`'s formula with GM 19.78 Msun, on cell centres; not checked against an independent number.
- **Grid-scale block copies:** blocks gid 0 and 1 have the same seed (corr 0.9999), because glibc treats `srand(0)` as `srand(1)`. They are radial neighbours, and it is the only such pair among the 193 seeded blocks. No half-wedge copies: |corr| <= 0.024 at 0.174 d in seeded shells.

## Observations so far, not interpreted
- Comoving L in the first cell falls from 1.00 L_* at t = 0 to 0.014 L_* by 0.0134 d, staying ~0.01 to 0.174 d. Lab L_in/L_* ~ -5.5, dominated by advection. Comoving L/L_* recovers outward: 0.54 at r = 22, 0.82 at 25.8, ~0.9 by 30.
- dE reaches -5.7e42 erg by 0.0067 d (from hst, without gravity), then falls slowly, to -1.02e43 at 0.174 d. L_top/L_* falls from 3.9 to 1.17. About -4.2e42 of this is in r > 60 by residual, which looks like optically thin IC radiation leaving; not tested.
- Unseeded zones have a theta-only v_r pattern at r = 41 (half-phi corr +1.0000, dominant mode (1, 0), 28 cm/s on a 4.5e4 cm/s mean). Possibly periodic theta boundaries on sin(theta) geometry; not tested.
- KE: 5.34e40 erg at 0.174 d, radial 3.62e40, lateral 1.73e40.
