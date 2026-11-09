# TASK 2026-11-09 viper -> ORION: raycheap3 -- why is the base VET/M1 time scheme 1st order? (CPU isolation)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap3.md` on this branch. User criterion: VET must be 2nd order
in space AND time. Our time-dependent test star shows the base scheme (hesdirk2 + extrap2, no rotation) at p = 1.0.

## 1. Binary
MOD = `fork/raycheap2-1109` head `9b9951606` (= your NEW of TASK raycheap + test inputs; the code is that of
`063d6d96e`): CPU MPI, PROBLEM none, your build recipe. It has the test-only key
`rad_m1/implicit_flux_x1min_mod_amp|period` (inner face flux F0 (1 + A sin(2 pi t/P)) set before each hesdirk2
stage solve at t_n + dt; for the hesdirk2 tableau both stage solves end at t_n + dt).

## 2. Decks (in MOD's tree, `tests_raycheap/iso/`): `iso_<arm>_c<cfl>.athinput`
Test star "hot" (cs lo, 24 MeshBlocks: run with 24 ranks), A = 0.5, P = 400, t_end 400, implicit_tol 1e-13,
implicit_lin_tol 1e-14, no ray rotation (vet_rotate_avg 1, rotate_every 0); cfl 0.15 / 0.075 / 0.0375 / 0.01875
(dt 1.34 / 0.67 / 0.335 / 0.167). Arms (one ingredient each, against B0):

| arm | change vs B0 | question |
|---|---|---|
| B0 | none (closure vet_col, time2_vet_col = extrap2, implicit_closure_lag = step, coupling on) | reproduces viper: p = 1.0 |
| C1 | vet_col frozen: time2_vet_col = lag, vet_col_every = 1000000 (built at the first step only) | is the 1st order in the transport/forcing itself? |
| C2 | closure = m1 (local analytic), implicit_closure_lag = pass | the no-VET baseline |
| C3 | closure = m1, implicit_closure_lag = step | the lagged local closure |
| G1 | coupling = false (no radiation-gas exchange; the hot clump stops emitting) | is it the gas coupling? |

If C2/C3 stop at start with `implicit_precond = rbgs needs implicit_krylov_fuse >= 1`, rerun them with
`rad_m1/implicit_krylov_fuse=1` (or `rad_m1/implicit_precond=jacobi`) on the command line and say which.
(viper already found: extrapolating all hand-over fields (time2_vet_x2_fields 11/13/15) and time2_vet_col = rebuild
do NOT change the order -- the closure update is not the culprit as far as that goes.)

## 3. Runs and analysis
`srun -n 24 ... athena -i iso_<arm>_c<cfl>.athinput -d <dir>` for the three coarse levels of every arm first
(c0.15, c0.075, c0.0375; ~10 / 20 / 35 min at 24 ranks), then the c0.01875 level only for arms whose
3-level p is >= 1.5. Per arm: `python3 tests_raycheap/conv.py <c0.15> <c0.075> <c0.0375> [<c0.01875>]`.
Report d, p per arm, t_end, wall, stage/BE/fallback counts (grep `time_scheme=hesdirk2:` and `NOT ADMISSIBLE`).
Interpretation we need: which arm(s) reach p ~ 2. No code changes; if every arm is p ~ 1 say so -- then the
1st-order term is in the implicit transport stepping itself (next TASK) or in the forcing (we will vary it).

## 4. NOTE back
Build sha/compiler, the per-arm conv.py lines with counts, wall times, any FATAL and the workaround used.
