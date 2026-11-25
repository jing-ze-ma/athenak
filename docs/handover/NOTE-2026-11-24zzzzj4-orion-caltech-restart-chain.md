# NOTE 2026-11-24zzzzj4 orion -> Caltech: production-restart test of the release candidate (b09a7b3c1), 02-09 now

USER asks when the updated scheme can be tested as a production restart: start on 02-09 as soon as the b09a7b3c1 star_vet (CUDA) is built and its tiny gate passes (the full suite stays on 02-03; if the build is not there yet, build first). Diagnostics only: production runs are not touched; use hard-linked copies of their rst (newest complete one), 4 ranks.
Order: m35zams, m35mms, bsg, t9l62, t9mlt (cheapest first; stop and post on the first FATAL).
Per star, 4 H200:
1. +1 step with `time/nlim` from the rst and the run's own command-line keys, else no changes; list every NOTE line of the restart rule ("the fix applies", "is the old behaviour", "added from the command line"), rc, FATAL / CW, clean_init passes, dt.
2. 50 cycles straight; 25 cycles + rst + 25 cycles from it: hst and final dump bitwise (state and hst).
3. against the running binary (c4speed or sdorder3, from the same rst, same 50 cycles): L_top, KE_int, dt, s/cyc, passes; the first 5 cycles must show the new defaults' expected change only (radial merge for bsg, FOFC flag exchange, wall flux, gtime_pos): say which of these NOTE lines apply.
4. Mass / energy per-cycle residual with your END-row formulas, and L_z (bsg).
Any failure: post the log and stop that star; go on with the next.

-- orion (coordinator)
