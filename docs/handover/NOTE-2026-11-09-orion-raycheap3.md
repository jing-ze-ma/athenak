# NOTE 2026-11-09 Orion -> viper: raycheap3 base-scheme time-order isolation (B0 C1 G1 C2e)

**Summary:**
- B0 and C1 are 1st order (p 1.01 / 1.07). G1 is below 1st order (0.43 / 0.56).
- **C2e (closure = eddington + implicit_hr_model = true) is 2nd order: p 2.24 / 1.80, with d ~ 3000x smaller than B0.**
- **Correction to the interim NOTE 9bc6c6545:** the hesdirk2 transport stepping and the inflow forcing are NOT the 1st-order source, since C2e uses the same stepping and the same forcing.
- The 1st-order term comes with what the vet_col arms bring in: the VET tensor and/or the vet_gd half-range face data. It is independent of how often the VET is updated, because C1 freezes it after step 1 and is still p = 1.
- C2e changes two things against B0 (the closure and the hr_model beam instead of the vet_gd ray data). A next arm that separates them would be vet_col + implicit_hr_model = true (no vet_gd half-range data), if the code allows it.
- Possibly relevant to raycheap4 (P3 = C2e): the C2e d values (7e-8) are tiny compared with B0 (2e-4). Check that the C2e solution actually responds to the A = 0.5 forcing at the same amplitude as B0, i.e. that the hr_model test star is not much less sensitive to the inner flux.

## Binary
- /orion/ptmp/jinma/raycheap3_1109/athena_9fbd2d3c6: CPU MPI, PROBLEM none, gcc 13.1, OpenMPI 4.1.8, Kokkos SERIAL, -ffp-contract=off.
- Its src/, CMakeLists and kokkos are identical to 9b9951606 (`git diff --stat` is empty); decks and conv.py are from 9b9951606.

## Deviations from the TASK
1. **No staging (USER: "make tests fast").** All 4 levels of every arm were submitted at once, so c0.01875 ran regardless of the 3-level p. conv.py is reported on both 3 and 4 levels.
2. **Layout.** A p.exclusive node is 112 physical cores, so 4 runs of 24 ranks fit per node. The account has a 6-job running limit. Jobs: B0 227947 and C1 227948 (one node each), G1 inside the 3-node job 227957, C2e 227958.
3. **C2 / C3 (closure = m1) cannot run on the cubed sphere at this commit.** Each workaround hit the next fatal error:
   - (a) `implicit_precond = rbgs needs implicit_krylov_fuse >= 1`;
   - (b) a command-line `rad_m1/implicit_krylov_fuse=1` on a fresh start is refused: `Parameter 'implicit_krylov_fuse' in block 'rad_m1' on command line not found`;
   - (c) the key in a deck copy: `implicit_flux_beam = halfrange on the cubed sphere needs vet_gd ... or implicit_hr_model`;
   - (d) `implicit_flux_beam = closure`: `implicit_blend = idort_f ... needs implicit_flux_beam = halfrange`;
   - (e) `implicit_hr_model = true`: `rad_m1_sph.cpp:198 ... cubed sphere supports only ... closure = eddington | vet_col`.
   - **Substitute C2e** = the C2 deck with `closure = eddington` + `implicit_hr_model = true`; viper accepted it in 4949ab359.
   - C3 has no separate meaning with a fixed Eddington tensor, since there is nothing to lag.

## Results (t_end = 400.0 at every level; rc = 0 and 0 FATAL in every run)
| arm | d (c.15/.075, .075/.0375, .0375/.01875) | p 3-level | p 4-level | walls c.15/.075/.0375/.01875 (s) | NOT ADMISSIBLE |
|---|---|---|---|---|---|
| B0 | 1.805e-4, 8.932e-5, 4.258e-5 | 1.01 | 1.01, 1.07 | 565 / 1015 / 1845 / 3662 | 0 |
| C1 (vet_col frozen after step 1) | 1.952e-4, 9.690e-5, 4.628e-5 | 1.01 | 1.01, 1.07 | 380 / 682 / 1229 / 2299 | 0 |
| G1 (coupling = false) | 2.529e-4, 1.878e-4, 1.271e-4 | 0.43 | 0.43, 0.56 | 584 / 827 / 1566 / 3048 | 6 / 0 / 0 / 0 |
| C2e (eddington + hr_model) | 7.165e-8, 1.518e-8, 4.372e-9 | 2.24 | 2.24, 1.80 | 1271 / 1970 / 2837 / 4360 | 1 / 0 / 0 / 1 |

time_scheme=hesdirk2 counters:
- **B0, C1:** stage steps 298 / 597 / 1195 / 2391; BE steps 1; fallbacks 0. C1_c0.15 has one time2_posfix re-solve of 40 cells.
- **G1:** c0.15 has stage steps 299, BE 13, fallbacks 5. The other levels have BE 1 and 0 fallbacks.
- **C2e:**
  - c0.15: stage steps 296, BE 3, fallbacks 1, posfix 1 re-solve (51472 cells).
  - c0.075 / c0.0375: BE 1, 0 fallbacks.
  - c0.01875: stage steps 2389, BE 3, fallbacks 1, posfix 1 re-solve (20 cells).
  - The finest-pair p of 1.80 may be affected by the BE steps / fallback at c0.01875, or by d approaching the 1e-9 level.

Bin outputs and test runs were deleted after conv.py. Logs are kept in /orion/ptmp/jinma/raycheap3_1109.
