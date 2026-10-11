# NOTE segfix (segfix-1011): the short-characteristics segment across an opaque -> thin face

Worktree /raven/u/jinma/ATHENAK/wt-segfix, branch `segfix-1011` (local, NOT pushed), on prodcand-1011 e66ee669.
Work dir /raven/ptmp/jinma/segfix_1011. Binaries: bin/athena_cpu_wip1 (35df5beb, split only), bin/athena_cpu_wip2
(a6a97649: split, emit and the cyl24 pgen).

Commits:
- 35df5beb `vet_seg = split` (vet_sc 4 paths, vet_gd shell sweep + band exit, vet_col_lat)
- a51a8467 copyright header
- f7974448 pgen `cyl` up to 24 objects (the 20 IDORT bubble sources; test-only, as suite's cyl24 build)
- a6a97649 `vet_seg = emit` (the limiter only where S_u > S_0) + line lengths
- 55f7343f comment line lengths only

## Bottom line

- **The bug is real and the fix works for hot emitters.** At a hot opaque slab (B = 1) next to vacuum:
  - with `mean`, J_sc in the first thin cell is 0.28 / 0.050 / 5e-3 / 5e-4 / 3e-5 at tau_cell 1 / 10 / 100 / 1e3 / 1e4;
  - with `split`, it is 0.449 / 0.4996 / 0.5000 / 0.5000 / 0.5000 (exact 0.5).
  - The moment E (REC-be, tJe) goes from 0.39 / 0.22 / 0.19 / 0.19 / 0.19 to 0.54 / 0.51 / 0.501 / 0.500 / 0.500.
- **But `split` is wrong for cold walls heated from the thin side.** In the crooked pipe it makes the walls black:
  - T >= 0.25 keV arrives at probe 3 at t = 1.9 instead of 0.3, and at probe 5 at t = 8.3 instead of 4.0.
  - A resolved 1-D wall (tau_cell 27 -> 0.43) shows `mean` is the right one there: E near the wall is -6.2 % (split) against +1.3 % (mean) at t = 0.2.
- **`vet_seg = emit` takes the best of both.** The split face value is used only where S_u > S_0 (the opaque cell emits into thinner gas); a heated cell keeps `mean`.
  - It equals `split` on every hot-emitter test (slab, hot cylinder, bubbles, vet_gd P1/P2).
  - It equals `mean` on the heated walls: the crooked pipe agrees to < 1e-3 at t = 10, and wall1d is identical.
- **Regressions:** none. Constant-opacity tests are unchanged to solver tolerance. The atmosphere, 1-D ladders, G1/G3/G5, crossing beams, static shadow and non-LTE change by <= 1e-3 relative, with the same orders.
- **Recommendation: make `emit` the default** for new runs of vet_sc / vet_gd / vet_col_lat; a restart keeps `mean`, as vet_bc_diff does. Do not use `split` in production.
  - The default is NOT switched in these commits.
  - `mean` is bitwise identical to e66ee669 (vet_sc and vet_gd + vet_col checked, all output arrays).

## 1. The scheme (src/rad_m1/rad_m1_vetseg.hpp)

The segment from the upwind point u to the cell centre 0 is split at the cell face. The upwind part has optical depth
hu = chi_u L_u and the own part h0 = chi_0 L_0 (vet_sc: L = dx1/(2|mu1|) each; vet_gd / vet_col_lat: the radial face
crossing on the ray, half the chord for a turning ray). S is piecewise linear in tau and continuous at the face value

    S_f = S_u + hu*g,   g = gf + th*(minmod(gb, gf) - gf)
    gf = (S_0 - S_u)/(hu + h0)          the mean segment slope
    gb = S_u - I_u                      the upwind cell's own slope (I - S = -dS/dtau inside an opaque cell)
    th = max(hu - h0, 0)/(hu + h0) * (1 - exp(-hu))   [emit: th = 0 unless S_u > S_0]

and the intensity is the two linear-in-tau sub-segments. The properties:
- **Constant chi:** th = 0, and S_f is exactly the mean scheme's face value, so the result is the old segment to round-off. G5 / non-LTE / crossing beams change by 1e-8 at most (solver tolerance).
- **Smooth fields:** th = O(dx), so the change is O(dx^2). In the 1-D prototype (ana/proto.py), with chi = 2e^{-x}, the mean/split difference falls 6.4e-4 / 1.9e-4 / 4.0e-5 at n = 50 / 100 / 200.
- **Opaque hot body / single hot cell:** gb = 0, or gb has the opposite sign to gf, so S_f = S_u and I ~ S_u.
- **Unresolved diffusion photosphere** (S linear in tau inside the opaque cell): gb = gf, so the mean result is kept. The prototype is exact for S = 1 + 1.5 tau with a cell-to-cell tau ratio up to 100.
- **Positivity:** S_f lies between S_u and S_0, and all weights are >= 0. The scheme is monotone and adds no negative intensities.
- **GPU determinism:** pure per-thread arithmetic on values already read (no new reads, atomics or stencil), so the halo_exact marks are unchanged.

**Keys.** `<rad_m1>/vet_seg = mean | split | emit` (default mean). It is read in VetInit (vet_sc) and in the vet_col branch (vet_gd, vet_col_lat).

**Not patched:** vet_col (the 1-D radial column closure, rad_m1_vetcol.cpp). It has the same mean-chi segment, but its sweep
needs the face fraction per (ray, shell) table; in production it covers only tau > 30 (smooth).

## 2. Thick -> thin tests

**Slab (Cartesian vet_sc, inp/slab.athinput, ana/slab.py).**
- Setup: `cyl` rect hot slab x < 0 (B = 1, held), near-vacuum x > 0 (rho 1e-10), 64 x 8 cells, 64 steps to steady state.
- Exact on the thin side: J = E = 0.5.

| tau_cell | J_sc mean | J_sc split | E REC-be mean | E REC-be split | E cen mean | E cen split |
|---|---|---|---|---|---|---|
| 0.1 | 0.447 | 0.469 | 0.488 | 0.509 | 0.476 | 0.488 |
| 1 | 0.278 | 0.449 | 0.395 | 0.536 | 0.378 | 0.436 |
| 10 | 0.0499 | 0.4996 | 0.221 | 0.513 | 0.108 | 0.114 |
| 100 | 0.0050 | 0.5000 | 0.191 | 0.501 | 0.013 | 0.013 |
| 1e4 | 3e-5 | 0.5000 | 0.1875 | 0.5000 | 1e-4 | 1e-4 |

- emit = split here (tau 1 / 100 / 1e4 rerun).
- **cen (central face flux) does not pick the fix up.** Its E stays at background: the central moment flux across an opaque -> thin face is a separate, moment-side defect. The REC arm (fse/hcs/jse thin switch) follows J_sc.

**Hot cylinder** (suite hcyl, REC-be, 256^2; exact by ray quadrature, ana/hcyl.py; c1000 exact ana/hcylc1000_exact_256.npz):

| case | L1(E) outside, mean -> split | E(R0)/ex | rF/rF_ex (2-3.5) |
|---|---|---|---|
| rho kappa 10 (tau_cell 0.31) | 0.120 -> 0.043 | 0.965 -> 1.010 | 1.03 -> 1.05 |
| rho kappa 1000 (tau_cell 31) | 0.61 -> 0.22 | 0.80 -> 0.99 | 0.58 -> 0.70 |

- f_xx(3) at rho kappa 1000 goes from 0.547 to 0.269 (exact 0.297).
- emit is identical to split.

**vet_gd sphere: P2 thermal channel beam** (sptests pinc, 64 x 48 x 48, 100 cycles; ana/pin_eval.py, PIN_TAG=c):

| arm | L_code/L_ex at r 4.3 / 4.6 / 4.9 | FS L_fs/L_ex | Picard mean/max |
|---|---|---|---|
| Ff (fse+hc) mean | 2.25 / 2.31 / 2.35 | 0.29 | 3.8 / 11 |
| Ff split = emit | **1.50 / 1.53 / 1.56** | 0.53 | 4.5 / 20 |
| Pf (hr) mean | 39.9 / 41.1 / 41.9 | 0.29 | 5.1 / 7 |
| Pf split | 41.2 / 42.4 / 43.2 | 0.53 | 5.1 / 7 |

- **The P2 over-luminosity is partly this bug.** The vet_gd formal solution itself was 3.4x too DIM through the channel (it under-emitted the hot core). The moment solve over-shoots it.
- split/emit takes the FS to 0.53 and fse+hc from 2.3x to 1.5x.
- hr's 40x is a moment-flux problem and is unchanged.

**IDORT multi-source bubbles** (suite bub, 512^2, ring 100, REC, t = 0.1 = 1 t_c; ana/bub.py):

| disc rho kappa | total E/exact at t = 0.005 / 0.025 / 0.1, mean -> split (= emit) | L1 at t = 0.1 |
|---|---|---|
| 1e5 (IDORT, tau_cell 195) | 0.64 / 0.40 / 0.25 -> **0.78 / 0.59 / 0.52** | 0.77 -> 0.67 |
| 1e3 (tau_cell 1.95) | 0.74 / 0.56 / 0.51 -> 0.84 / 0.71 / 0.68 | 0.62 -> 0.76 |

- At the IDORT opacity the total roughly doubles: the discs now emit.
- The test still FAILS. The lit fraction (0.58 against 1.00) is unchanged, and the remaining deficit is the moment front problem of NOTE-suite 3.4.
- At rho kappa 1e3 the total improves but L1 grows (0.62 -> 0.76): E near the discs is closer, the front still lags. This needs a separate look.

**Crooked pipe** (suite pipe, REC, c dt/dx 220, 512 x 256, t = 10; ana/pipe.py). Probe T (keV) and arrival of T >= 0.25 keV:

| arm | T(t_end) probes 1-5 | arrival 1-5 | wall s |
|---|---|---|---|
| mean | 0.469 0.428 0.387 0.344 0.277 | 0.1 0.1 0.3 0.7 4.0 | 2584 |
| split (t_end 8.4, time limit) | 0.458 0.432 0.375 0.342 0.254 | 0.1 0.2 **1.9 4.2 8.3** | >3600 |
| emit | 0.469 0.428 0.387 0.344 0.277 | 0.1 0.1 0.4 0.7 4.0 | 2738 |

- At t = 1, split gives probes 3-5 of 0.093 / 0.053 / 0.050 keV against mean/emit 0.328 / 0.279 / 0.212.
- The walls are tau_cell 27 cold gas heated from the pipe. Split treats their unresolved surface as the cold cell average, so they become black absorbers.

**Which is right for a heated wall: wall1d** (1-D analogue: E = 1 bath | pipe gas 0-1 | wall 1-1.5; pipe opacities and c_v; 108 x f cells, f = 1/4/16/64, so the wall tau_cell goes 27.8 -> 0.43; ana/wall1d.py). E against f64:

| t | f1 mean | f1 split | f1 emit | f4 any |
|---|---|---|---|---|
| 0.2, x = 0.99 | +1.3 % | -6.2 % | +1.3 % | <= 0.1 % |
| 0.5, x = 0.99 | -1.0 % | -3.4 % | -1.0 % | <= 0.05 % |
| 1.0 | -0.4 % | -0.4 % | -0.4 % | |

- The resolved wall skin heats almost at once (its heat capacity per optical depth is tiny), so its surface is near the incident field. That is the mean (linear-in-tau) face value.
- In 1-D the bath dominates. Around the pipe's bends the wall re-emission is the only feed, hence the large pipe effect.

## 3. Regression with split (and emit where run)

All with bin wip1 (wip2 for emit), arms as in the source dirs; values are mean -> split.

- **Grey atmosphere, Cartesian** (deepatm od.py, 32..2048, be). Self-convergence per tau band and Hopf L1 are identical to <= 0.5 % of the error.
  - tJe, tau > 10: 2.36e-3 / 5.95e-4 / 1.49e-4 / 3.72e-5 -> 2.36e-3 / 5.94e-4 / 1.49e-4 / 3.72e-5.
  - Hopf L1 at n 512: 3.632e-4 -> 3.632e-4.
  - Orders unchanged; cen likewise.
- **1-D ladders** (pulse k 0.128 / 12.8 / 1280, rw t10 / t1000; X/T/C, be + hesdirk2): identical, except rw_t10 X be at 512, 3.93e-4 -> 3.43e-4.
  - The Pvs arms are FATAL in both builds (implicit_vimp needs gas_feedback), pre-existing.
- **G1** (atm2d, c dt/dx 1e2 / 1e4) L1 6.996e-4 -> 6.995e-4 (tJe), hr 1.954e-3 identical. **G3** Marshak: PASS, L1 identical.
- **Fixed G5:** identical ratios (constant opacity), all PASS. The old G5 in the battery FAILs at k 12.8/128 in both builds (pre-existing).
- **Crossing beams** cb RECbe r12: identical (L1 E 0.570, J exact).
  - cb REC (hesdirk2) r4 / r12 DIVERGE at cycle 2 with mean too, on this prodcand base (suite's 4c033dcf ran r12). This is pre-existing; please flag it.
- **Static shadow** shd REC r12g2: L1 E 0.011 both; E arrays differ by L1 5e-5.
- **Non-LTE** nltek1-8: max rel E diff <= 1.1e-8; L1 vs exact identical (0.0111 ... 0.0001). ddiff, sdifk, eq1: identical.
- **vet_gd P1** patm_s Ff (4 x 4 lateral cells, N 64 / 128 / 256): E/J-1 at R_ph -6.3e-3 / -1.9e-3 / -1.9e-3 -> -3.9e-3 / -1.35e-3 / -1.7e-3.
  - Bands and orders unchanged; the thin-envelope plateau is unchanged (lateral resolution).
  - emit = split.
- **Cost:** Picard is unchanged on the regressions. On P2 it is 3.8 -> 4.5; on the pipe (emit) 3.01 -> 3.02.

## 4. Files

- ana/proto.py: 1-D prototype.
- ana/slab.py, ana/hcyl.py (HCYL_CHI/HCYL_SFX), ana/wall1d.py, ana/cmp.py, plus the suite's bub.py, pipe.py, shd.py, cb.py, s2.py.
- inp/slab.athinput, inp/wall1d.athinput; jobs_*.txt, run*.sbatch.
- runs/: slab, hcyl, wall1d, bub, pipe, cb, shd, nltek*, ...
- da/: deepatm copy with S arms. cb/: gates (G1/G3). g5/: fixed G5. sp/: the sptests driver with S/E arm suffixes (runs/, runs_base/ = e66ee669).
- Bitwise check: runs/*/*base* and sp/runs_base against mean. All output arrays are identical.
