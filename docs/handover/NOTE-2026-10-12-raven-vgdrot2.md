# NOTE vgdrot2: making the vet_gd quadrature rotation time-accurate

Worker Z, 2026-10-12, MPCDF Raven. Work dir `/raven/ptmp/jinma/vgdrot2_1012`.
Branch **vgdrot2-1012**, from prodfinal-1011 4253af30, in worktree `/raven/u/jinma/ATHENAK/wt-vgdrot2b`. (`wt-vgdrot2` was already in use by the older vgdrot2-1011.) Not pushed.

Commits:
- a78e1977: `vet_gd_rotate_dt`, `vet_gd_rotate_t0`, `vet_gd_rotate_align`
- 3f8c5fd6: `vet_gd_rotate_slope`
- 1cabbd5b: `vet_gd_rotate_stretch`
- ccb7396e: `vet_gd_rotate_slope` defaults to `keep`

The diff is localized:
- rad_m1.cpp / rad_m1.hpp: the keys.
- rad_m1_vetgd.cpp: VetGdAngle / VetGdBlock / VetGdAlignDt / VetGdRotStep / prefetch angle / print.
- mesh.cpp: 2 lines, the dt alignment before the tlim clamp.
- driver.cpp: 2 lines, VetGdRotStep before Time2StartProbe.
- rad_m1_time2.cpp: 2 lines, the probe gate.

Binaries in `bin/`:
- r1 = a78e1977
- r2 = 3f8c5fd6 (probe default)
- r3 = 1cabbd5b
- r4 = ccb7396e (final)

The CPU, he-CPU and A100 binaries are built from each.

## 0. Why the cycle clock is first order
`vet_gd_rotate_every = N` picks the set from `cycle/N`. Two consequences:
- At a given physical time, every dt level samples a different set.
- The number of jumps grows like 1/dt. As dt goes to 0, the solution tends to the angle-averaged closure, and the convergence to that limit is first order.

This is T's finding on the rlat_Tr restart ladder:
- rotate 10: 1.86 / 1.48 / 1.14
- rotate 0: 2.0

## 1. Design evaluation
| design | verdict |
|---|---|
| (a) physical-time clock: block = floor((t + dt/2 - t0)/dt_rot), same golden-ratio angles | **implemented** (`vet_gd_rotate_dt`). The closure is piecewise constant in t and identical on every dt. Restart-safe. |
| (a') + `vet_gd_rotate_align` (default on with rotate_dt) | **implemented**. NewTimeStep reaches each t0 + k dt_rot in ceil(rem/dt) equal steps, so a jump always sits on a step boundary. |
| (b) continuous phi = Omega t | **rejected**. Every build becomes a rotation: rotate_every 1 costs 0.558 s/cycle against 0.394 (+41 %) on 4 A100. The wall maps are mostly nearest-direction (exact 3840 of 72960 on rlat), so the closure is not smooth in alpha anyway. |
| (c) blend old and new closure over a ramp | not implemented. It needs two sets swept per build during the ramp, and the exact-list cache holds one angle (a 0.11 s batch rebuild per switch). This is the only cheap-to-state cure for the FSAL remainder (sect. 3); estimated at about +13-20 % if done as always-two-sets. |
| (d) two staggered sets per build, no rotation | measured as gl_nphi 24 rotate 0 (240 directions): 0.446 s/cycle against 0.372 rotate 0 / 0.394 rotate 10 (**+20 % / +13 %**). Time-accurate by construction. Its ray-effect quality against rotation was not measured. |
| (e) stage consistency | holds. Every build in a step sees the same t and dt (and the same ncycle): predict at stage 1, rebuild at stage 2, in-solve rebuilds, and the probe. So both hesdirk2 stages use one set. |

## 2. Order results (orders over pairs 1-2 / 2-4 / 4-8 / 8-16; finest-pair error in [ ])
Orders are computed with each level's actual mean step, because aligned ladders do not have exact ratios of 2 (`sp/ana/ordh.py`, `ana/ordt.py`).
Arm names:
- rt = dt_rot of 10 level-1 steps, aligned
- rt3 = 3.7 steps, aligned (non-commensurate)
- rt3n = 3.7 steps, not aligned

### 2a. T's rlat_Tr restart ladder (conservative2, slope start)
The r1 binary, cycle clock: RZ reproduces T's errors bitwise.

| arm | dens | E | F1 | F3 |
|---|---|---|---|---|
| RZ (rotate 10, cycle) | 1.91 1.50 1.15 [2.05e-8] | 1.92 1.32 1.04 [1.63e-8] | 1.63 1.21 1.06 [1.70e-6] | 1.62 1.12 0.76 [5.66e-7] |
| RZr0 (rotate 0) | 2.11 2.03 1.57 [8.90e-9] | 2.10 1.83 1.92 [5.36e-9] | 2.01 1.93 1.66 [4.62e-7] | 2.04 2.06 1.97 [8.59e-8] |
| RZt (rt) | 2.08 1.98 1.80 [7.88e-9] | 2.04 1.81 1.95 [5.33e-9] | 1.96 1.94 1.66 [4.61e-7] | 2.04 2.06 1.98 [8.65e-8] |
| RZt3 | 2.12 1.97 1.83 [7.55e-9] | 1.98 1.88 1.75 [5.79e-9] | 2.00 1.87 1.71 [4.43e-7] | 2.02 2.08 1.97 [8.48e-8] |
| RZt3n | 2.12 1.93 1.90 [7.54e-9] | 2.12 1.82 1.71 [6.17e-9] | 2.00 1.92 1.57 [5.01e-7] | 2.03 2.04 1.84 [9.72e-8] |

### 2b. Fresh ladders (X's harness, HFfZ = candidate keys)
- **rlat_k100_T** (rlat at rho kappa_s 100):

  | arm | F1 | E |
  |---|---|---|
  | cycle clock | 1.47 0.59 nan [1.05e-4] | |
  | r0 | 1.94 1.95 1.96 [1.96e-6] | 1.93 1.88 1.92 [1.50e-8] |
  | rt | 1.94 1.94 1.91 [2.04e-6] | 1.92 1.89 1.86 [1.55e-8] |
  | rt3 | 1.50 1.75 1.56 [4.07e-6] | |
  | rt3n | 1.61 1.70 0.54 [9.18e-6] | |

  The leftover error grows with the number of jumps: rt3 has 2.7x as many as rt.
- **rswt_T** (relaxed source): rt = r0, both F1 [2.7e-6]; the cycle clock gives [3.3e-5]. The remaining floor is X's relaxed-source floor.
- With `time2_vet_col_source = state`:
  - r0: 2.0
  - rt: dens 2.40 2.17 1.77, E 2.01 1.98 1.91, F1 2.02 1.97 1.75 [5.4e-7 vs 4.0e-7]
  - rt3: F1 1.38 [7.6e-7]
  - rt3n: F1 1.43 [1.0e-6]
- **rsw_T** (kappa_s 1000): rotation has no effect (all arms 2.0, identical).
- **rlatt_om5_T + state** (thin lateral wave, rho kappa_s 10, rotating gas):

  | arm | E | F1 |
  |---|---|---|
  | r0 | 1.90 1.97 2.03 [1.06e-8] | 1.93 1.97 1.98 [3.5e-6] |
  | cycle clock | 0.69 0.10 0.61 [7.1e-6] | [4.5e-3] |
  | rt | 1.32 1.16 1.08 [1.3e-7] | 1.20 1.10 1.04 [7.0e-5] |
  | rt3 | similar [1.4e-7] | [1.6e-4] |
  | rt3n | 0.65 [3.0e-7] | [5.4e-4] |

## 3. The thin-regime remainder: the FSAL slope
hesdirk2 is FSAL. On the first step of a new block, stage 1 uses K1 = f_old(U^n), the old set's slope. This leaves a local error of O(dt (f_new - f_old)) once per jump. A fixed number of jumps then gives a FIRST-order global term, however well the jumps are aligned. The term is small where the closure jump barely changes f (thick: rlat, rswt), and visible in the thin lateral wave.

Tried `vet_gd_rotate_slope`:
- `probe`: on the first step of a block, drop t2_ok and take the time2_start=slope probe solve with the new set.
- `be`: on the first step of a block, take a backward-Euler step.

Results (r2 binary, `sp/runs2`, `runs/rlat_Tr2`):
- probe:
  - rlatt E is fixed (2.10 2.22 2.35), but F1 and v stay first order (F1 [1.07e-4]).
  - **The stiff cases break:** rlat_k100 rt 2.0 → 1.0 [2.5e-4]; rlat_Tr RZt → 1.0 [4.2e-4].
  - The BE-probe slope (Y_p - U^n)/dt_p is not f(U^n) in cells with lambda dt >> 1.
  - Cost on 4 A100: a rotation cycle goes from 0.516 to 0.931 s.
- be: rlatt F1 0.94 [3.4e-4] (Lie split).

**Default is `keep`.** A real fix needs either:
- the slope corrected by f_new(U^n) - f_old(U^n), which needs an f-evaluation in the implicit operator (rad_m1_implicit, T's area); or
- a closure continuous in t, i.e. (c) / (d).

## 4. Ray effects and AG Car (rt against the cycle clock)
- With a constant dt, rt (dt_rot = 10 dt) rotates at the same instants as rotate_every 10.
- Held-gas sptests at 4253af30 + r1:
  - T6a thin_L16_u 128: E L1 6.16e-3 for both (r0 6.21e-3).
  - T6b core_L16_u 128: exterior +7.85e-3 against +7.92e-3 (r0 +7.83e-3).
  - T7 shd 32: umbra <E>/E_free 0.049 for both, penumbra and lit identical (r0 0.036, lit L1 0.39-0.56).
  - P2 pinc 64: rt against r10 differs by 5e-7 in F1 (r0 against r10: 2e-3).
- So rt keeps exactly the ray-effect behaviour of rotate 10. (These steady tests are not sensitive enough to favour rotation over r0.)
- AG Car A reduced CPU (480x32x32, 120 cycles, A0c candidate deck, job 31092624):
  - Picard: 3.03 (r10), 3.03 (rt), 2.93 (rt20), 2.88 (r0).
  - L_top/L_in = 1.0000x for all arms at the common dump times.
  - The photosphere lateral rms of E (4.05e-3 at t 1.25e5) and of F_r are identical.
  - r0 against r10 differs by 5e-5 in E(R_ph): the reduced run is insensitive to rotation over 120 cycles.
- **Alignment overhead seen:** the AG Car dt drifted from 1252.56 to 1252.4, below dt_rot/10 = 1252.5, so each interval took 11 steps (dt 1138.6). The run reached t = 1.423e5 at cycle 120 against 1.503e5 (-5 %).

- r4 check (job 31093489, A0c; the deck then named slope = probe for rt / rts, so these two ran WITH the probe):
  - `vet_gd_rotate_stretch = 0.02` (rts) removes the overhead: t = 1.5030e5 at cycle 120, the same as r10 (rt: 1.4233e5).
  - Picard: 3.22 with the probe against 3.03 for r10, 0 NON-CONVERGED, stable.
  - E(R_ph) differs from r10 by about 1e-6, and L_top/L_in is the same.
- Regression: r4 (default keep) against r1 is bitwise in every bin dump (data after `<par_end>`; the headers differ only by the new GetOrAdd keys) on rlat_k100_T rt l1/l2 and rlatt_om5_T rt+st l1/l2 (108 files). The cycle clock (RZ) reproduces T's c2order errors exactly.

## 5. Cost (full AG Car A, 4 A100 gpudev, 32 cycles, s/cycle avg)
| arm | job | avg | steady | rotation cycle |
|---|---|---|---|---|
| r10 (cycle clock) | 31092290 | 0.3947 | 0.381 | 0.516 |
| rt (physical, same instants) | 31092290 | 0.3946 | 0.381 | 0.515 |
| r0 | 31092290 / 31092618 | 0.3721 / 0.3713 | | |
| rotate_every 1 (b proxy) | 31092290 | 0.5576 | | |
| gl_nphi 24, rotate 0 (d proxy) | 31092618 | 0.4463 | | |
| rt + probe | 31092904 | 0.4369 | | 0.931 |
| rt20 + probe | 31092904 | 0.3956 (1 rotation in window) | | |

Alignment step overhead, from replaying recorded AG Car dt(t) series (`ana_align.py`; runB, rampA, reatm long A):

| dt_rot | stretch 0 | stretch 0.02 | stretch 0.05 |
|---|---|---|---|
| 10 x median dt | +3.2 to 5.7 % steps | +1.2 to 4.0 % | ~0 to 0.8 % |
| 20 x median dt | +1.3 to 3.5 % | +0.4 to 1.9 % | ≤ 0 |
| 40 x median dt | +1.0 to 2.3 % | | |

`vet_gd_rotate_stretch = s` lets the aligned step be up to s longer than the CFL dt instead of adding one extra step per interval.

## 6. Proposal
- Production: `vet_gd_rotate_dt = 20 x typical dt` (AG Car A: 2.5e4; B: about 20 x its dt), `vet_gd_rotate_align = true`, `vet_gd_rotate_slope = keep`, optionally `vet_gd_rotate_stretch = 0.02`.
- Effects:
  - The rotation frequency per step is half today's, so per-cycle cost goes from 0.394 to about 0.383 s.
  - The step overhead is about 1-2 %.
  - Rotation becomes a function of t, so restarts and dt changes no longer change the closure history.
  - The time order matches rotate 0 in the thick and moderate regimes.
  - In the thin regime a first-order FSAL term remains. It is 65x below the cycle clock and scales with 1/dt_rot.
- Keep `vet_gd_rotate_every` as the fallback default for restarts (bitwise).
- If thin-regime order matters, (d) gl_nphi 24 with rotate 0 is time-accurate at +13 % over rotate 10. Its ray-effect benefit against rotation still needs a measurement on T7 / P2.
