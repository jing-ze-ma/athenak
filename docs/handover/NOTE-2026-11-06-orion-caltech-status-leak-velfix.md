# NOTE 2026-11-06 Orion -> Caltech: status of TASK-2026-11-06t; a USER velocity-fix test is also running

No new results here, this is so you know what is coming and don't duplicate it.

## TASK-2026-11-06t, all running
1. **Column, cfl_rad 0.011, 0.05, 0.1, 0.25, at 1e-6 and 1e-9 (job 227416):** the three larger-step 1e-6 runs have reached 1 d; the 1e-9 runs are at ~0.05-0.09 d. Results when all are at 1 d.
2. **Short 3-D, 1/16-lateral wedge, 112 ranks (227414 at 1e-6, 227415 at 1e-9):** 1e-6 needs ~2 h to 0.25 d. 1e-9 runs at ~10 s per cycle with ~59 iterations, so ~12.5 h to 0.25 d. **It may be stopped near 0.1 d** if that holds, to free a job slot; I'd say so if it is. Early check: the wedge's KE x 16 at 1e-6 is 2.2e40 at 0.027 d, against 2.1e40 in `prod_seeded`.
3. **`prod_seeded` 1 d data:** about 1.4 h away, in `orion3d/`, with dE, dU, dEgrav. It's at 0.72 d with KE 7.06e40, past its 7.29e40 peak at 0.56 d.

## Also running: USER's proposed Athena++ fix (job 227417)
The USER asked to try fixing Athena++ itself: **correct the velocity in the implicit solve, use a tight error_limit, and use the normal CFL step.** A separate source copy adds an opt-in `radiation/vel_mode` that takes the same-step pre-radiation kick out of the frame velocity, in the 4 x 4 column only. No results yet, and the design isn't reviewed yet.

| Arm | cfl_rad | error_limit | vel_mode |
|---|---|---|---|
| A | 1 | 1e-9 | 0 |
| B (the proposal) | 1 | 1e-9 | 1 |
| C | 1 | 1e-6 | 1 |
| D | 0.025 | 1e-9 | 1 |

The test is whether B matches the converged small-step reference (`el1e-9`) at dt ~20x larger, with no wall film to 10 d. It doesn't change `prod_seeded`, which is still Linhao's deck at 1e-6 and cfl_rad 0.025.

If you have already tried this in AthenaK or in a toy, or know of a reason it should fail, please say so here.
