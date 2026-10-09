# NOTE 2026-10-09 Caltech -> viper: BSG true repro stays exactly symmetric -> proposal: restart fresh with a tiny seed (user OK to ask)

Run: chain 4286643-45 (bsg_1009/prod), t = 6.24 d of 57.4, clean (0 FATAL/NaN/NON-CONVERGED, dt 85 s, ~1.24 s/cycle).
Page: https://claude.ai/artifact/BQ8LZRU8uV1Zw83GtZQQxK (viper page pipeline, t 6.0 d).

## Evidence: no lateral motion at all
- user.hst: lateral KE = KE_int - KEr_int <= 1.9e25 erg, i.e. <= 2.6e-16 of KE_int over all 540 rows (non-zero in only 55
  rows, pure round-off of the difference); the latest rows are exactly 0.
- Dumps (t 4.5 and 5.5 d): Fig 4 shell r = 40.59 Rsun has v_r = -2.448 km/s, rho = 4.606e-8, T = 1.804e5 K identical in
  every cell of the float32 dump; the intensity map at each column's tau_R = 1 is flat (all at 59.22 Rsun, rms 0.000).
- KE_int: 1.16e41 (1 d) -> 4.5e40 (2 d) -> 2.0e40 (4-5 d) -> 1.1e40 (6 d): the start-up radial adjustment decaying, no growth
  (paper ~1e43 erg by ~20 d).
- Interpretation: laterally uniform IC + periodic theta/phi + deterministic per-column arithmetic -> the symmetry is kept
  EXACTLY, so there is no round-off seed for the instability to grow from; it will stay 1-D. Same behaviour as repro_4n
  (symmetric to ~1e-6, KE 4e39 by day 7), here even more exact. The paper's convection without a stated seed most likely
  came from its numerics (120 discrete ordinates not invariant under lateral shifts, SMR/MPI block structure).

## Proposal (user 10-09: "push the note to viper"; decision is viper's / the user's)
- Stop the unseeded chain and restart FRESH (t = 0) with a tiny seed that mimics numerical noise rather than the old 1 %:
  e.g. `he_seed = 1.0e-6` (same 16-mode T seed at fixed rho, he_seed_nk 16), everything else unchanged; or tell us if
  he_seed can be applied on a restart (then seed the 6 d state instead of restarting).
- Caltech keeps the current chain running until viper answers. Binary/input unchanged otherwise (6c5d8fb2_nofma, md5 e9873137);
  a seeded fresh start would need only the one key in the input file (command-line keys absent from the input are FATAL).
