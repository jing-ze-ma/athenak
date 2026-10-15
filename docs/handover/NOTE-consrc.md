# NOTE: conservative gas source of the implicit M1 solve (consrc-1011)

Branch `consrc-1011` (from prodcand-1011 686f6a24). Work dir `/raven/ptmp/jinma/consrc_1011`.
Binaries: `bin/athena_cpu_cs5`, `bin/athena_he_cpu_cs5`. Baseline: prodcand `athena_cpu_558be1c2` (code
parent of 686f6a24).

## 1. What it is

New key `<rad_m1>/implicit_gas_src`:

- `matched` (default). The current term-by-term deposit. It is bitwise identical to prodcand: eq3t
  cen be and tJe h2v+vstage give 704/704 identical tab files against `athena_cpu_558be1c2`, and the
  reduced AG Car A (30 cycles) has identical hst and bin data against `athena_he_cpu_558be1c2`.
- `conservative`. The gas takes minus the radiation's own source of the converged solve.
- `conservative_full`. As `conservative`, but the AP-HLL / blend part of a face is also charged in
  full (section 4). This is not recommended.

Momentum (per cell and direction, code units, radiation momentum F/(c chat)):

    dm = -(F'_c - F^old_c)/(c chat) - 1/2 sum_{interior faces} Z_f
         + sum_{physical x1 faces} [ (F0'_b - F0old_b)/(2 c chat) + w dt k_t F0'_b / c ]

- `F_c` is the stored cell flux (u0). `F^old` is the old vector:
  - BE: u0 at entry.
  - hesdirk2 stage: u0 at stage start plus the stage increment of the cell flux, `csinc = c1 K1_F +
    c2 K2_F`. The slopes are `csk1`/`csk2`, stored by the write-back and combined in
    Time2FormStage.
  - The first solve of a run takes the face state `½ΣF0n + a E^n`, because the problem generator's F
    is not the radiation's state.
- `Z_f` is written by the face kernels of the solve itself: x1 in step (g), x2/x3 in
  ImplicitTransverseTerms. They are the same kernels that set the stored F0, at the same iterate, so
  the identity is exact by construction. This holds for the F0/A decomposition, the sph / vet_gd
  curvature and lateral terms (gr+od, M1SphCurv/Lat/Tan), the cs cross terms, and block/MPI faces,
  which both blocks evaluate from the same halo data. For a face equation
  `F0' = Θ(F0old - chat c dt G - chat dt v g0) + H` (G = gr + od is the face's transport,
  H = the AP-HLL part, ω = Θ(1 + chat dt k_t)), the definition is

      Z_f = -(1 + chat dt k_t) H/(c chat) + ω dt G + (1 - ω) F0old/(c chat)

  This is the face share `-(F0' - F0old)/(c chat) - Z_f`, which equals
  `ω dt/c (k_t F_d + v g0) + dt k_t H/c` plus everything else the face equation removed (the
  implicit_rad_inertia `-Θ dA`, the memory weight, ...).
- On central faces (ω = 1, H = 0), Σ_cells (dm + ΔF_c/(c chat)) = -dt Σ_f G_f. This telescopes, so
  total momentum is conserved to round-off whatever terms the face equation contains.
- A changes of the cell (`a E`, including the velocity change of the kick under vimp/vstage) are in
  ΔF_c, so they are charged too.
- Physical x1 faces have no face equation. Their F0 change is inflow, and they keep the matched
  half-cell drag (a boundary source).
- Energy: the gas total energy changes by `-(c/chat)(SRCR - SRCB E') + W`, with W the exact kinetic
  energy change of dm. E pays `(chat/c) W`, the same structure as matched.
  - A Krylov-residual term, `R = b - A E'`, was tried and dropped. The measured energy defect is
    Picard-level (the lagged row parts), not linear-residual-level.
- F' depends on dm through the clip and the vimp re-forming. There are 3 fixed-point passes, and the
  last dm is the one of the stored F', so conservation stays exact. A clip that binds only after the
  positivity moves gives its momentum to the gas, with the kinetic energy change taken from eint.
- gas_feedback = false: unchanged. Matched path and nofb work. ddn cs = mt to all printed digits.
- In conservative mode, the M1_RCP (implicit_realisable_coupling) scaling, the matched v g0 deposit
  and the a·W/c² fixed point are not used: ΔF_c covers them.

Code: commits b5b92187, d9949457, 5632ff94 (+ the NOTE commit). Files:

- `rad_m1_impl_wb_kernel.hpp`: conservative block, final-clip fix, cell-F slope.
- `rad_m1_implicit.cpp`: key, Z in the face kernels, allocation, the first-solve flag.
- `rad_m1.hpp`: members.
- `rad_m1_time2.cpp`: Time2FormStage, csinc. This is one block before `t2_solve =`, and the predict
  path is not touched.

## 2. Moving equilibrium eq3t (32×32, c = 10, t = 1; domain totals from 32 x2-row tab slices, double)

The target (code frame) is E 1.1326, v 2.9554, P 3.

| arm | matched E / v / max\|ΔP\| | conservative E / v / max\|ΔP\| | max\|ΔE_tot\| (cs) |
|---|---|---|---|
| cen be | 1.13333 / 2.95510 / 1.6e-4 | 1.13296 / 2.95527 / 8.9e-16 | 6.0e-11 |
| cen hesdirk2 | — | 1.13296 / 2.95527 / 8.9e-16 | 3.2e-12 |
| cen ovc | — | 1.17692 / 2.95344 / 1.3e-15 | 8.0e-11 |
| tJe be | — | 1.13296 / 2.95527 / 8.9e-16 | 1.4e-11 |
| tJe h2+vimp+vstage+predict | 1.13329 / 2.95512 / 1.7e-4 | 1.13295 / 2.95527 / 8.9e-16 | 7.3e-11 |
| hr be | 1.12609 / 2.95855 / 3.3e-3 | 1.12585 / 2.95866 / 3.4e-3 | 1.2e-11 |
| hr be, conservative_full | | 1.13394 / 2.95500 / 1.3e-15 | 2.5e-11 |
| cen be, tol 1e-12 | ΔE_tot 2.7e-15 | 8.9e-16 | 2.2e-14 |

Findings:

- Total momentum is conserved to round-off on every central arm, in BE and hesdirk2. Matched leaks
  1.6e-4.
- E and v move toward the analytic target. For E: 1.13296 vs 1.1326 (matched 1.13333).
- Total energy: the defect is the same order as matched and scales with the Picard tolerance (6e-11
  at tol 1e-8, 2e-14 at tol 1e-12). It comes from the Picard-lagged parts of the E row, which are not
  exactly telescoping between passes, and not from the source assembly. It is not round-off at the
  production tolerance.
- hr (blend) keeps its 3.4e-3 leak in `conservative`. The AP-HLL part of a face has no momentum
  equation and keeps the drag-only deposit. `conservative_full` closes it, see section 4.

## 3. Other tests (conservative vs matched, same binary)

- ddn (gas_feedback = false): identical. ddf be / h2v / hvy / hvyh2v, cen and tJe: peak, centroid,
  L1 and D_eff/D are identical to 4 digits (e.g. cen h2v 0.9700 / -0.0001 / 0.0237 / 0.9843).
- nltek1/4/8 and sdifk (static): bitwise.
- eq1: E and eint bitwise. Only noise-level v (1e-16) and F (1e-12) differ.
- Marshak (cen, hr, tJe): agree to 1e-12 relative (dens, eint, E, F1). With the first version, which
  had no boundary drag, the boundary cell v differed by 2 %.
- eq2 (uniform, hot gas, VET noise): total momentum matched -1.5e-7, cs 2e-11. Cell v (6e-3 noise)
  differs by 2 %, dens by 7.5e-5.
  - Side finding: on the hr arm a uniform medium at rest develops a mean F ≈ 27 (E 150, c 100) and a
    gas v 0.02 that grows linearly, in BOTH paths. This is a radiation-side asymmetry of the
    blend/half-range face. It is not this change, but it is worth a look.
- Order ladders: rw_t10a / rw_t1000a × {cen, tJe} be + {cenPvs, tJePvs} hesdirk2, X/T/C, 120 rows
  per path, in `ladder/lad_{mt,cs}.txt`.
  - All orders of E, v, rho and eint are unchanged (largest change 0.03).
  - The one exception is F1 in T mode on rw_t10a hesdirk2: cenPvs 2.29/2.01 becomes 1.46/1.67 at
    the two finest levels (L1 6.7e-6 → 1.3e-5); tJePvs 1.94 becomes 1.66. F1 X/C orders are
    unchanged.
  - The likely source is the radiation-inertia term `d(a(w)E)/dt` that the gas now pays at each stage
    (it is O(E/ρc²)). This is not resolved.
- Reduced CPU AG Car A (480×32×32, 4 ranks, hesdirk2 + time2_vet_col predict + vimp + vstage +
  fse/hc + realisable coupling, 600 cycles to t = 7.44e5 s):

  | arm | dt | Picard mean / max | L_top | budget residual ∫(dEtot - L_in + L_top - P_esrc + P_sponge) |
  |---|---|---|---|---|
  | matched | 1240 s | 3.87 / 5 | 2.8444e38 | -1.19e38 (5.6e-7 of ∫L_top dt) |
  | conservative | 1240 s | 3.87 / 5 | 2.8444e38 | -1.19e38 (identical) |

  - Both arms are stable with the same Picard and Krylov counts (4720/4719 outer passes). Final dumps
    differ by ≤ 1e-7 relative, except v1 at 0.8 % of max|v1| (outer thin region).
  - The first version, without the first-solve rule, charged the problem generator's F at cycle 1: a
    one-time budget jump of 1.7e40. It is fixed and shows that the rule is needed.
  - Radial momentum has gravity and geometric sources, so only the energy budget is checked here.

## 4. conservative_full

It is exact on eq3t hr. On AG Car A it is not usable: dt falls from 1240 s to 137 s by cycle 368,
and the Picard mean rises from 3.3 to 3.9. The HLL part of a thin face carries the front's own
(non-momentum-equation) flux changes, and charging them to tenuous gas gives spurious kicks. It is
kept only as a diagnostic.

## 5. Recommendation

- Make `conservative` the default for gas_feedback = true. On every central / VET-central test it
  conserves total momentum to round-off, it moves eq3t toward the analytic state, and it changes
  nothing measurable in the static, diffusion and Marshak tests or in AG Car.
- Open before switching: the hesdirk2 F1 time order in rw_t10a T mode (1.7 vs 2.0), and the energy
  defect, which is still Picard-level (unchanged from matched).
- moveq patches under `conservative`:
  - `implicit_rad_inertia`: the write-back half (the v g0 deposit and the a·W/c² fixed point) is
    redundant and is skipped. The face half (F0 -= Θ dA) is still radiation dynamics: it decides how
    much of dA stays in F0 (thin) and how much goes to the gas (thick). Conservation no longer
    depends on it, so it can stay or be retired on accuracy grounds only. That was not tested here.
  - `implicit_f0_init`: no longer needed for conservation, because the first-solve rule makes the
    cold start non-creating. It still sets the initial lab F the radiation evolves from (physics, not
    bookkeeping), so keep it.
  - `implicit_nofb_work`: still needed. The gas_feedback = false path is untouched.
- hr/blend arms keep an O(1e-3) momentum leak on eq3t from the AP-HLL part, as in matched. Fixing
  that needs a momentum-consistent thin-face flux, not a source change.
