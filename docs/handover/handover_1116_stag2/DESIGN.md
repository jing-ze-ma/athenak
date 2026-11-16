# stag2: a 2nd-order, realizable, exactly condensable staggered scheme (Orion design, 2026-11-16)

Goal: keep staggered's assets (better smooth accuracy, clean thick time order, AP without Rhie-Chow, scalar-E solve)
and remove its blocker (no realizable low-order fallback, no realizability test that accepts transported fronts).
Status: design + Python prototypes only (clean-cc belongs to viper; nothing goes into C++ before the USER decides).

## 0. Why proto2/stag.py failed, and why this is different
stag.py carries F_f as a free unknown with its own cell "flux" Psi_i = c^2 D_i E_i - alpha (a/2)[F]_i. That has two
consequences. (i) Its low-order form upwinds neither E +- F/c nor anything else that is monotone, so even the all-low
step is not realizable (3883/6644 steps). (ii) The F-F coupling forbids condensation.
stag2 never gives F an F-dissipation. Every face flux is an explicit function of E (plus lagged/old data), so the stage
system is ALWAYS a scalar E system, as in production and the C++ staggered path.

## 1. Unknowns and stage equations (H-ESDIRK2 / SDIRK2 stage, h = gamma dt)
Cell:  E_i (1 + h c kE_i) + (h/V_i) sum_f A_f F_f = Et_i + h c kP_i a T_i^4     (Et = explicit stage part)
Face flux, blended per face f:  F_f = F^L_f + b_f (F^H_f - F^L_f),  b_f in [0, 1]
- High order (2nd order, AP, keeps F memory):
  F^H_f = theta_f [ G_f - h c^2 (grad_n (D E))_f + C_f ],  theta_f = 1/(1 + h c chi_f)
  G_f = explicit part of F (old F + stage history). In the thick limit this is -(c/chi_f) grad(DE), the compact diffusion,
  with no RC and no dissipation; in the thin limit it is the staggered wave operator.
- Low order (kinetic, positive, no memory):
  F^L_f = c (h+_L E_L + h-_R E_R),  h+ >= 0, h- <= 0 from the rays' half ranges over J (isotropic: +-1/4).
  It is the first-order upwind discrete-ordinates flux taken in moments.
- C_f = consistency bracket (lagged): F_ray,f - [F^H_f evaluated on the rays' own J]. The rays' face flux comes from the
  sweep's face intensities. It is RHS-only, so condensation is kept.
- Gas: momentum source from face F_f (face-centred -chi F/c, split to the two cells). v.F and the Ma2025 beta terms are
  lagged in the gas fixed point, which also keeps condensation.

## 2. Realizability definition (staggered-native; accepts transported fronts)
- R1: E_i >= 0 for every cell.
- R2: donor bound on every face. If F_f > 0 then F_f <= c E_L; if F_f < 0 then -F_f <= c E_R. The flux leaving a cell
  cannot exceed c times that cell's energy density.
- There is NO co-located cell-F test. Nothing in the scheme needs a cell F, except diagnostics and the gas, and the gas
  uses face F. That test was what flagged every beam front in proto2.

## 3. Positivity guarantee (the missing piece in proto2)
- With b = 0 on every face and backward Euler (Et = E^n >= 0), the E matrix is an M-matrix. The diagonal is
  1 + h c kE + (h c/V) sum A h+_own, and the off-diagonals are -(h c/V) A h+_nbr <= 0. So E >= 0 for any dt.
- R2 also holds for b = 0. F^L is bounded by c times the donor E only if h+ <= 1 and |h-| <= 1, which holds by construction.
- Safety-net chain:
  1. stage with b = 1 - w(C1) (the smooth ray weight; w -> 1 in beams);
  2. if R1/R2 fail: set b = 0 on the bad faces and re-solve, at most 3 rounds;
  3. if R1/R2 still fail: BE with the marks;
  4. if R1/R2 still fail: BE all-low. This step is guaranteed.
- Arm FCT: instead of step 2, a linearised implicit FCT. Solve BE-L, then add a Zalesak-limited antidiffusive flux
  b_f (F^H - F^L) with bounds from R1/R2. Compare the limiter work against marking.

## 4. Multi-D
- x faces carry F_x, y faces carry F_y. The off-diagonal stress (D_xy E) in the x-face flux needs corner E, which gives
  a 9-point stencil that is not an M-matrix.
  - Arm A: lag D_off in the Picard loop (RHS).
  - Arm B: D_off in the matrix, limited against the diagonal (design doc 2.2a).
  The L flux uses only normal half ranges, so it stays an M-matrix in any dimension.
- Cubed sphere (later): face-normal F needs one frame per face, which is natural for staggered.

## 5. Tests (pass/fail before anything goes further)
1-D (prototype stag2.py in the proto2 framework):
- the proto2 table (thin/thick waves, time order thin and thick, radiating wave ka 1 / 100, odd-even, conservation,
  176-case realizability sweep, beam front width, absorber leak) vs cell A3 and the old stag.py;
- limiter work;
- linear-system size and iterations.
2-D (beamreal framework lf2d.py / raycl.py):
- narrow normal beam, crossing 45 deg beams, shadow;
- R1/R2 violations per stage and after the chain, BE fraction, all-low fraction;
- score vs rays;
- D_off arms A/B, bracket on/off.
Targets:
- smooth 1-D orders >= 1.9 in space and in time (thick included), with accuracy at stag.py level;
- 0 R1/R2 violations after the chain in every case;
- beam limiter work <= cell A3;
- every solve a scalar E system.
