# NOTE: time order of the implicit VET scheme (rad_m1, hesdirk2), torder-1010, 10-10

Branch `torder-1010` (local, not pushed) in worktree /raven/u/jinma/ATHENAK/wt-torder, based on
fork/xthinfix-1009 e50ec7d6. Fix commit: **8ee1e217**. CPU binary: /raven/ptmp/jinma/torder_1010/bin/athena_cpu_8ee1e217
(md5 79b9a167e9c7122707ef2f8b5e03ddd3; gcc/13 + openmpi/5.0, MPI, Release, no PROBLEM).

## Verdict

**There were two root causes. Both are now fixed by options, and with them every T ladder reaches 2.**

1. **The vet_sc tensor is frozen at U^n for both hesdirk2 stages.**
   - Both implicit stages sit at t^{n+1} (c2 = c3 = 1).
   - Stage 1 builds the formal solution at U^n: rad_m1_implicit.cpp:9432 calls `Time2VetStart` (rad_m1_time2.cpp:444, which reads the hydro u1 = U^n at :459). Stage 2 keeps that tensor.
   - The extrapolation that would correct it is off by default: rad_m1_time2.cpp:159 sets `t2_vext = false`, so r = 0 at :593.
   - D is therefore lagged by O(dt). That gives an O(dt) error in K2 and K3, an O(dt^2) local error, and a **first-order** step wherever D changes. Pulses and the thin and intermediate regions are affected; the thick limit, where D = I/3, is not.
   - **Fix (8ee1e217):** a new key `<rad_m1>/time2_vet_sc = lag | predict | rebuild`.
     - `predict` builds the step's single formal solution at the predicted state P = U^n + dt K1 (the FSAL slope). The gas there is the Heun predictor + dt K1 of the coupling. P is at t^{n+1} to O(dt^2).
     - `rebuild` uses predict for stage 1 and adds a second build at Y1 for stage 2.
     - The code reuses the vet_col machinery (`Time2VetColAt`, which now calls VetShortChar for vet_sc).
     - The key is read only when named. The default is `lag`, which is bitwise unchanged: cen and sq reproduce the e50ec7d6 numbers digit for digit, and G1 / G3 / G5 are identical.
2. **rw: the battery runs with `implicit_vimp = false` (od.py BASE), the explicit-velocity coupling. It is first order in time at P_rad/P_gas = 100.**
   - Even the Eddington closure gives 1.00 on rw_t10 and rw_t1000.
   - The production default for hesdirk2 is `implicit_vimp = true` (rad_m1_implicit.cpp:1788). With it, rw_t1000 is about 2.0.
   - rw_t10 also needs (1) for dens and E, and `time2_vstage = true` for F1. time2_vstage is an existing key; its default is true only on the sph wedge. Without it F1 is 1.0, as already documented in rad_m1_time2.cpp.
3. **The rw_t10 E and F1 T ladders at amplitude 1e-5 are not measurable below about 1e-4. This is a ladder artifact, not the scheme.**
   - At amplitude 1e-5, the E perturbation is about 7e-8 of the background E.
   - The Picard test (|dE|/E) and the Krylov stop (lin_tol x max|b|) are both relative to the background. So once the error drops below about 1e-4, the solver noise dominates. At small dt, Picard is about 1.0 and BiCGStab about 0.85 iterations per solve.
   - The symptoms:
     - every second-order arm's rw_t10 E / F1 at fine dt is erratic;
     - this is the source of sq's "hd2 T -0.12 -0.16" on rw_t10 in the Raven battery NOTE;
     - tightening the tolerances (implicit_tol 1e-14, lin_tol 1e-15, EW off) does NOT cure it, because it hits round-off.
   - At amplitude 1e-3 (new cases rw_t10a / rw_t1000a) the same arms give clean orders of 2.0.

**Recommended settings for 2nd order in time:** `time2_vet_sc = predict` + `implicit_vimp = true` (already the default) + `time2_vstage = true` (not yet the default off the wedge). On top of these, use a ladder amplitude of 1e-3 for rw.

The default is not flipped on this branch. Flipping `time2_vet_sc = predict` (and time2_vstage) to the default for fresh hesdirk2 + vet_sc runs is recommended once a GPU production smoke run (AG Car) is done. A restart without the key would keep `lag`, as time2_vet_col does.

### Evidence: T orders (hesdirk2, nx1 = 256, or 128 for rw; dt/1..16; last three orders)

The `be` row is for reference. Old = e50ec7d6 binary; new = 8ee1e217. Arms are on top of cen = `closure=vet_sc vet_tensor=full implicit_flux=central`.

| experiment | pulse_k0.128 E | pulse_k12.8 E | pulse_k1280 E | rw_t10 (a = amp 1e-3) | rw_t1000 |
| --- | --- | --- | --- | --- | --- |
| cen (baseline, lag, vimp off) | 1.79 1.71 1.57 | 1.05 1.03 1.02 | 1.96 1.94 1.90 | E 1.01 1.00 1.00, dens 1.00 1.00 1.00 | E/dens 1.01 1.00 |
| be (cen) | 0.91 0.95 0.97 | 0.91 0.95 0.97 | 0.99 0.99 1.00 | E 0.87 0.93 0.96 | 0.97 0.99 |
| (a) closure = eddington | **2.01 2.01 2.00** | **2.00 2.00 2.00** | 1.98 1.99 2.00 | E 1.02 0.99 1.02, dens 1.00 (vimp off) | 1.01 1.01 (vimp off) |
| (a) eddington + vimp | n/a (static gas) | n/a | n/a | dens 2.23 2.02 2.14 | E 1.98 2.21 |
| (b) time2_vet_extrap = true | 1.94 1.98 2.00 | 2.07 2.04 2.02 | 1.98 1.99 2.00 | E 1.03 1.01 1.00 | n/a |
| (c) radiation only: pulses (gas_feedback off) | = cen | = cen | = cen | rw gas_feedback off: E diverges (no wave) | same |
| (c) vimp on (lag) | n/a | n/a | n/a | 10a: E 2.93 0.96 0.98, dens 1.85 1.90 1.73 | 10a/1000a E 2.04 2.07, dens 1.99 2.00 |
| (d) trbdf2 tableau | 1.81 1.74 1.61 | 1.07 1.04 1.02 | 1.96 1.95 1.92 | E 1.01 1.01 1.00 | n/a |
| (e) tight tol (1e-14, lin 1e-15, EW off) | = cen, bitwise | = cen | = cen | 10: still 1.0; 10 vimp: erratic (round-off) | = |
| **FIX cenP (time2_vet_sc = predict)** | **1.94 1.97 (2.00 at 16, tight tol)** | **2.06 2.03 2.02** | **1.98 1.99 2.00** | vimp off: still 1.0 (cause 2) | 1.01 (cause 2) |
| cenB (rebuild) | 1.94 1.97 | 2.04 2.02 2.01 | 1.98 1.99 2.00 | 1.0 (vimp off) | 1.0 |
| **cenPv (predict + vimp)** | n/a | n/a | n/a | 10a: **E 2.67 2.01 2.00, dens 1.95 2.08 2.08**, velx 1.87 2.02 1.97, F1 1.07 1.04 1.02 | 1000a: **E 2.05 2.09, dens 1.99 1.99** |
| **cenPvs (predict + vimp + vstage)** | n/a | n/a | n/a | 10a: **E 2.67 2.01 2.00, F1 1.58 2.32 1.97, dens 2.08 2.08, velx 2.02 1.97** | 1000a: **E 2.05 2.09, F1 2.00 1.99, dens 1.99** |
| cenPvs, original amp 1e-5 | | | | 10: dens 2.19 1.99 2.12, E 0.26 3.85 1.69, F1 noisy (cause 3) | 1000: E 1.97 2.21, dens 1.88 2.65 |

sq and hr checks (new binary):

| arm | pulse k0.128 T | k12.8 T | k1280 T | rw_t10a T (with vimp + vstage) | rw_t1000a T |
| --- | --- | --- | --- | --- | --- |
| sq (lag) | 1.79 1.71 | 1.05 1.03 1.02 | 1.96 1.94 1.90 | sqv0: E 2.94 0.96 0.98, dens 1.78 1.74 | 2.04 2.07 |
| **sqP** | **1.94 1.97** | **2.06 2.03 2.02** | **1.98 1.99 2.00** | sqPvs: E 2.00 2.01, dens 1.96 2.06, F1 1.98 1.83 1.64 | 2.05 2.09 |
| hr | 1.24 1.13 | 1.03 1.02 1.01 | 2.08 2.55 -0.97 | hrv0: E -0.24 0.08 0.40 (no convergence) | 2.09 1.88 |
| hrP | 1.25 1.13 | 1.25 1.11 1.06 | 2.10 2.70 -1.17 | hrPvs: same as hrv0 | 2.10 1.88 |
| hrx0P (hr, xthin 0) | 1.94 1.97 2.00 | 2.06 2.03 2.02 | | | |

**No regression of X / C, gates or cost:**
- **X ladders.** cenP and sqP equal cen and sq to 3 digits on every pulse. rw_t10a / rw_t1000a X with vimp is 1.96-2.07.
- **C ladders improve:**
  - k12.8 C: cen 1.68 1.95 1.67 -> cenP 1.90 2.24 2.24.
  - rw_t1000 C: 1.05 -> 1.97 (cenPv).
  - rw_t10a C, cenPvs: E 2.15 1.74 2.03, dens 2.05 2.02 2.03.
- **Gates** (G1 cfl 1e2 / 1e4, G3 c1 / c10, G5, with cen / cenP / sq / sqP / hr / hrP): all rc 0.
  - Hopf L1, G3 L1 and the G5 ratios are identical to the old binary to 3-4 digits.
  - NON-CONVERGED is unchanged, i.e. the known G1 cfl 1e4 at maxit 30.
- **Cost.**
  - Picard per solve: lag vs predict is identical to 3 digits in every ladder (pulses 1.995-1.999, rw_t10 1.33 / 2.0, rw_t1000 1.92-2.0). BiCGStab inner iterations are equal or slightly lower. NON-CONVERGED 0 everywhere.
  - predict does ONE formal solution per step (the same as lag), plus one fused kernel for T and the opacities at P.
  - rebuild does two formal solutions per step and gives no better orders than predict.

## Remaining issues

- **hr stays first order in T (1.0-1.25) with predict.** hr differs from cen on the pulses only by the `implicit_blend_xthin` override. That override's face weight depends on dt (X = (c dt/dx)/(1 + c dt chi_f), M1HrFace in rad_m1_implicit.cpp, about line 296). Under a T ladder the operator therefore changes with dt. With xthin = 0, hr equals cen and goes to 2.00.
  - hr1 (pure half-range, berthon) is 1.04-1.21 at k0.128 but 2.0 at k12.8. Its absolute T error is 70x smaller than cen's. Not investigated further.
  - hr on rw_t10 does not converge in T (dens 0.4-0.5), with or without the fix. This is the known hr rw_t10 failure.
- **Opacities and the EOS inside the stage solve are still evaluated at the stage START state.**
  - Stage 1 uses the gas Heun predictor, which lacks the coupling increment dt K1.
  - With T-dependent opacities, this is the same O(dt) lag as the tensor. The tested cases have constant kappa, so they cannot show it.
  - The fix would be to evaluate the solve's opacities at P as well, or to repair `implicit_opac_update` (HANDOVER-09-26: "v_h 0.88 from the frozen vet_sc tensor + opacities").
- **rw at amplitude 1e-5 cannot measure E / F1 time errors below about 1e-4** (cause 3). The battery should use amplitude 1e-3 (inputs cpu/inp/rw_t10a, rw_t1000a) and `implicit_vimp = true` (or not force it off in od.py BASE).
- **Not tested:**
  - GPU;
  - 2-D / 3-D production (AG Car);
  - restart (time2_vet_sc keeps no extra restart state; K1 is already restart state);
  - `implicit_mr_every` + predict (the code path is that of vet_col and handles mr_dt);
  - vet_sc_every > 1 with predict (refused, FATAL).
- **Multi-block / MPI check: PASS.** See the MPI section below. The pulse and rw cases cannot be split, because vet_x1_periodic needs a single MeshBlock.

## Details

### Setup
- Worktree: `git worktree add -b torder-1010 /raven/u/jinma/ATHENAK/wt-torder origin/xthinfix-1009`, with kokkos copied from xthinfix_1010/src/kokkos.
- Run tree /raven/ptmp/jinma/torder_1010:
  - cpu/ holds copies of od.py (new arms and cases), one.sh, fixkeys.py, pybin and inp, plus inp/rw_t10a and rw_t1000a (`radwave_amp = 1e-3`).
  - run/ holds the old-binary diagnosis (e1, e2).
  - run2/ holds the new-binary runs (e3-e6, eval_e3/e4/e6.txt).
  - gates_run/ holds the gates (battery_torder.sh).
  - logs/.
- Slurm jobs (general, 72 cores): 31050708 (e1), 31050889 (e2), 31051135 (build), 31051238/43/47/53 (e3), 31051265 (gates), 31051909/10 (e4), 31051919 (e5), 31052354/55 (e6).
- **Ladder hygiene (e).**
  - T mode: dt = T/(nT l), which is constant and hits tlim exactly (pulse_k12.8 l = 4: 160 cycles, dt 1.25e-3).
  - The space error is common to all levels (same grid), so it cancels in successive differences.
  - The reference is the next finer level, not the finest one.
  - The pulse ladders are not tolerance-limited: the tight-tolerance arms are bitwise equal, or equal to 3 digits.
  - rw at amplitude 1e-5 is tolerance- and round-off-limited (cause 3).
  - rw_t1000 level 1 is hydro-cfl-limited, so its orders use levels 2-16 (as in the battery).

### Arms (cpu/ana/od.py, all on top of cen)
| arm | keys |
| --- | --- |
| cenE | closure = eddington |
| cenX | time2_vet_extrap = true |
| cenG | gas_feedback = false |
| cenS | implicit_tol 1e-14, implicit_lin_ew_max 0 |
| cenR | time2_tableau = trbdf2 |
| cenv / cenvs | implicit_vimp = true (+ time2_vstage = true) |
| cenEv / cenEvs | eddington + vimp (+ vstage) |
| cenP / cenB | time2_vet_sc = predict / rebuild |
| cenPS / cenPv / cenPvS / cenPvs / cenPvsS | predict + tight tol / vimp / both / vimp + vstage / all |
| sqP, hrP, sqPv, hrPv, sqPvs, hrPvs, sqv0, hrv0, hrx0P, hr1P | sq / hr analogues |

### Why predict, not the extrapolation
`time2_vet_extrap` also gives 2.0 on the pulses (row (b)). However, it was measured to fail G1 at P = 100, tau = 10 (tests_m1/runs_3x_hesdirk2: 1.48). Its history (vet_prev) also depends on the previous dt. It extrapolates a realizability-clipped tensor, and it is first order after every backward-Euler fallback.

predict evaluates the true formal solution at a state consistent to O(dt^2), with no history. It is the default mode on the sph wedge for vet_col (time2_vet_col = predict, m1-sp-order2b).

### MPI / multi-block (mpi/mpi2.sbatch, job 31054010, interactive partition)
- Case: Marshak G3 c10 (hesdirk2, vet_sc, cen), widened to mesh nx2 = 32 on x2 in [-4, 4] (the same dy as the gate).
- Runs: 1 rank with one MeshBlock, against 4 ranks with 4 MeshBlocks along x2 (nx2 = 8 each).
- Result: the final E and F1 are **bitwise identical** between np1 and np4, both for lag (cen) and for predict (cenP).
- cenP vs cen differs by an L1 of 1.3e-5 in E and 1.6e-5 in F1, i.e. the change in the time error.
- Picard mean is 2.5347 (lag) vs 2.5423 (predict); NON-CONVERGED 0; time2_vet_sc fallbacks 0.
