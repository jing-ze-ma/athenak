# NOTE rotslope: removing the first-order FSAL term of the vet_gd rotation (vet_gd_rotate_slope = correct)

Worker R2, 2026-10-12, MPCDF Raven (one GPU cost run on DeltaAI). Work dir `/raven/ptmp/jinma/rotslope_1012`.

Branch **rotslope-1012**, worktree `/raven/u/jinma/ATHENAK/wt-rotslope` (kokkos is a symlink to wt-prodfinal; not committed).

Commits:
- 1bfe213a (on method-1012 e39ce690, pushed by main): `vet_gd_rotate_slope = correct`. Branch `rotslope-1012-v1` still points to it.
- Rebased onto prodfinal-1011 ed32f0fa:
  - 5f0a8c5b: the same change.
  - a088a43c: a Picard-pass counter for the correction solves, and the align counter relabelled.
- Main merged a088a43c into bsgcand-1012 d3bf42d7.

Binaries in `bin/`:
- c1 = 1bfe213a
- c2 = rebased onto f9d56e2d
- c3 = a088a43c

Each tag has CPU, he-CPU and A100 builds. On reduced AG Car A, keep and correct with c3 match c1 bitwise in the data (only the headers differ, by the new keys). Keep with c1 matches F's 78b2c961 FIN run bitwise in all 17 m1 dumps.

## 0. Test or method? (answered before any code was written)
The rlatt_om5_T remainder is a METHOD term, not a ladder artefact.
- Same rotation instants on every level: vet_gd_rotate_dt is fixed at 0.15625 (10 x the level-1 dt), t0 = 0, aligned. Every level made 25 rotations, and dt_l divides dt_rot exactly.
- Not start-up: rotate 0 is 2.0 on the same ladder.
- Per-jump scaling (job 31104094, Z's r4 binary, fixed physical dt_rot). Orders are for the level pairs 2-4, 4-8, 8-16; finest-pair error in brackets:

  | jumps | E | F1 | velx |
  |---|---|---|---|
  | 0 (rotate 0) | 2.03 [1.06e-8] | 1.98 [3.53e-6] | [4.8e-6] |
  | 1 (t = 2) | 1.93 [1.14e-8] | 1.87 [3.92e-6] | [5.1e-6] |
  | 3 (t = 1, 2, 3) | 1.85 [1.23e-8] | 1.20 [1.47e-5] | [1.7e-5] |
  | 25 (rt) | 1.08 [1.3e-7] | 1.04 [7.0e-5] | [2.2e-5] |

  The error grows with the number of jumps (and is larger for later jumps). This is the stale FSAL slope K1 = f_old(U^n) at each jump.

## 1. Method (`rad_m1_time2.cpp` Time2RotCorrect; `rad_m1_vetgd.cpp` VetGdRotStep; `driver.cpp`, 2 lines)
`vet_gd_rotate_slope = correct` (opt-in; `keep` stays the default):
1. On the first step of a new physical-time block, before the Heun stages, take one stage-type implicit solve.
   - Size dt_c = `vet_gd_rotate_slope_frac` x dt (default g = 1 - 1/sqrt 2, the stage size).
   - It starts from U^n (radiation, faces, gas; no hydro) and uses the NEW set.
   - Its old vector is U^n - dt_c K1 (the cell F of conservative2 likewise: -dt_c csk1).
2. Its slope (Y - old)/dt_c is stored as K1.
3. To first order, the stored slope is K1 + (I - dt_c A_new)^{-1} (f_new - f_old)(U^n):
   - where dt_c A is small (thin), it is the exact closure jump of the slope;
   - where the cell is stiff, it is keep's K1. keep is 2nd order there, and the L-stable stages forget a stage-1 slope error.
4. Why Z's `probe` broke the stiff ladders: it replaced the whole K1 by a BE slope (I - dt_p A)^{-1} f, which is wrong where A dt_p >> 1. Here only the jump is filtered.
5. Only the jump steps change, so an O(dt) error of the correction stays O(dt^2) globally.

Implementation details:
- The solve runs as a STAGE1 solve: its slope lands in t2k2 / csk2 and is copied to t2k1 / csk1.
- implicit_predictor is off for that solve, so ipred is kept. U^n is put back afterwards.
- `vet_gd_rotate_slope_reuse` (default true): the step's stage-1 solve keeps the correction solve's Time2VetColAt(1) build, so there is no extra vet_gd sweep.
  - That build is the same set at P = U^n + theta dt K1_old, with the gas at U^n instead of after Heun stage 1. The difference is O(dt), on jump steps only.
  - Reuse off gives the same errors (rlatt_om5_T E [1.37e-8] vs [1.41e-8]).
- A failed correction solve keeps K1 (counted). No correction runs after a stage fallback (t2_ok false).
- Report line: `vet_gd_rotate_slope=correct solves= failed= Picard passes=`.

Design (2), blending the closure over k steps, was not implemented, and (3), design (d), was not run: (1) meets every target.

## 2. Orders (FIN keys, F's harness `sp/anaF`, arm `HFfZ+vet_gd_rotate_slope~correct`; T ladders from the actual mean dt (ordh.py); finest-pair L1/mean in [])

### 2a. Fresh T ladders, levels 1-16, rotate_dt aligned at 10 x dt0

rlatt_om5_T (thin), orders for pairs 2-4, 4-8, 8-16:

| arm | E | F1 | velx | dens |
|---|---|---|---|---|
| keep | 1.17 1.06 [1.77e-7] | 1.11 1.04 [9.80e-5] | 1.01 [3.0e-5] | 1.46 [2.7e-8] |
| **correct** | **2.11 2.25 [1.41e-8]** | **2.09 2.16 [3.78e-6]** | 1.20 [2.25e-6] | 1.74 [1.56e-8] |
| rotate 0 | 2.13 2.41 [1.32e-8] | 2.14 2.27 [3.15e-6] | 1.31 [2.67e-6] | 1.84 [1.86e-8] |

correct - keep halves per level (1.48e-3 ... 9.6e-5 in F1): the first-order term is removed.

Other fresh T ladders:
- rlat_k100_T (stiff), F1:
  - keep 1.93 1.92 1.87 [2.14e-6]
  - correct 1.95 1.97 1.95 [1.95e-6]
  - every variable is 1.9-2.2
- rlat_T: correct = keep to 4 digits, dens 1.84 1.82 1.68 [7.65e-7]; rotate 0 gives 1.84 1.80 1.70 [7.71e-7] (start-up layer, F's known item).
- rlat_om5_T: correct = keep, E 1.89 1.82 1.62 [4.67e-7]; rotate 0 gives [4.71e-7].
- rsw_T, rswv_T (auto f0 and f0 zero), dif_r1_T, relax1_T: identical to keep.
- marsh_r1_T: flat (known open), errors 7.0e-6 vs 6.6e-6.
- rswt_T, rswtv_T:
  - correct - keep is O(dt^2) (dens 2.0e-7, 4.6e-8, 1.1e-8, 2.8e-9 at levels 1, 2, 4, 8), from the closure build at a slightly different predicted state.
  - Finest errors are equal or smaller. rswt_T: E [1.42e-8] vs [2.35e-8], dens [6.2e-9] vs [1.27e-8].
  - Orders on levels 1-8 differ within noise:
    - rswt_T E: correct 1.37 / 2.20 vs keep 1.89 / 2.03
    - rswtv_T F0 E: correct 1.24 / 1.59 vs keep 1.44 / 1.73 (errors 6.2e-8 vs 5.7e-8)

### 2b. Restart T ladders (F's pre-run rst files, levels 1-16, rord.py)

| ladder | keep | correct | rotate 0 |
|---|---|---|---|
| rlat_Tr dens | 2.07 2.03 1.77 | 2.07 2.03 1.79 | 2.09 1.97 1.99 |
| rlat_Tr E | 2.01 1.92 1.69 | 2.01 1.92 1.67 | 2.07 1.89 1.81 |
| rlat_Tr F1 | 1.95 1.94 1.66 | 1.95 1.94 1.66 | 2.00 1.93 1.66 |
| rsw_Tr | 2.0 every variable | identical | identical |
| rswv_Tr F0 | dens 2.03 2.01 1.68, E 1.93 1.71 1.56 | identical (to 3 digits) | identical |

### 2c. X ladders and exact-solution tests
Every X ladder is identical to keep to 3-4 digits: rlat, rlat_om5, rlatt, rlatt_om5, rsw, rswt, rswv, rswtv, pulse_r1, rswt_q. rlatt X F1 0.04 is F's open X-ladder item and is unchanged.

The exact-solution tests are identical too:
- patm / patmT orders
- thin_u / thin_L16_u: E L1 3.93e-6
- core_u / core_L16_u
- dif_r1: 2.0
- T7 shd: umbra 0.144
- P2 pinc: L/L_ex 0.614 / 0.624 / 0.634
- pinc: 1 NON-CONVERGED with keep against 2 with correct, both within its known run of non-admissible stages; the result is unchanged.

## 3. Production path
**Reduced CPU AG Car A** (480x32x32, 300 cycles, 14 rotations; job 31104737):

| quantity | keep | correct |
|---|---|---|
| rc | 0 | 0 |
| Picard mean / max | 3.149 / 6 | 3.234 / 6 |
| NON-CONVERGED | 0 | 0 |
| correction solves / failed | | 14 / 0 |
| L_top/L_in, end (min) | 1.00022 (0.99965) | 1.00022 (0.99965) |
| E/J-1 at R_ph | -1.61e-4 | -1.61e-4 |
| wall | 733 s | 726 s |

correct - keep:
- E, shell mean at R_ph: 8e-8.
- E, pointwise max: 6.3e-6 at 1.004 R_ph, growing by about 4e-7 per rotation. This is the stale-slope term at the thin/thick transition; for comparison, rotate 0 vs rotate 10 differs by 5e-5 there (Z).
- dens: only in the floor/sponge gas above 1.03 R_ph (at most 6e-4); below 0.98 R_ph it is 1.2e-7.

**GPU cost, full AG Car A** (DeltaAI 4 GH200, job 3366761, 120 cycles, 1bfe213a; Raven gpudev copy cancelled):

| arm | avg s/cycle | steady | rotation cycles | Picard |
|---|---|---|---|---|
| keep | 0.2162 | 0.2133 | 0.281 | 2.89 |
| correct | 0.2251 (**x1.041**) | 0.2187 | 0.370 | 3.02 |

There are 5 correction solves, none failed.

**Where the extra Picard passes go** (per-solve pass log, reduced CPU, 45 cycles):
- Each correction solve takes about 3.5 passes.
- The 2-3 steps after a rotation take about 4 more passes in total than keep.
- The same happens with implicit_predictor_order = 1 (+16 against +19 passes) and with reuse off (276 against 280; keep 261). So neither the order-2 predictor nor the reused build is the cause.
- The trajectory now carries the new set's slope at once, while keep's lags by a step. This is a property of the more accurate start, not a defect.
- The cost scales with 1/dt_rot: vet_gd_rotate_dt = 40 dt would halve it.

The align counter (`rotate_align shortened dt=`) counts NewTimeStep calls, including 2 at start-up, so 122 for 120 cycles is not a bug.

## 4. Verdict / proposal
- correct removes the thin-regime first-order rotation term: rlatt_om5_T reaches the rotate-0 order and error.
- It leaves the stiff and thick ladders, the restart ladders, the X ladders, the steady tests and reduced AG Car A unchanged (pointwise 6e-6 at the photosphere).
- Cost: +4 % average on GPU at rotate_dt = 20 dt; CPU wall is equal.
- Proposal: `vet_gd_rotate_slope = correct` for the vet_gd production decks that use vet_gd_rotate_dt, with keep as the code default for bitwise restarts. Whether to make it the fresh-run default is main's call (defaults-1012 policy).

## 5. Files
- sp/anaF: F's FIN harness, plus the rlat `kNN` token and ordh.py.
- Job lists: sp/jobs_q.txt, jobs_T2.txt, jobs_X.txt, jobs_v1.txt (the jump-count check, with Z's harness in sp/ana).
- Run trees: sp/runs_c (correct + controls), sp/runs (jump check), sp/runs_f -> suiteF runs_fin (keep).
- rst/: rstC.sh, ana/, runs/<case>/FINC_* (F's runs symlinked).
- agc/:
  - inputs: inp/agcA_{FIN,FINC,FINCr0}.athinput
  - runs: runs/{FIN,FINC}, *_c3, *_plog, *_plogp1, FINCr0_plog
  - agcdiff.py
- deltaai/: the TASK; the result is in bsg-files-1009 docs/handover/rotslope-1012/.
- Scripts: poolF.sbatch, agc.sbatch, rst.sbatch, gpu.sbatch, bld*.sbatch.
