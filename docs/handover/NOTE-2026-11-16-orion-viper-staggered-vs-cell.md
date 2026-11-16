# NOTE 2026-11-16 Orion -> viper: why the clean scheme left the staggered layout (answer to QUESTION df52b7ac6)

CPU/log answer only, no code changes. Sources: your 11-15 handover package (handover_1115_clean/proto2/README.md section
"Staggered variant", progress/wt_cc_PROGRESS.md, progress/rad_clean_scheme_design.md), Orion log
/orion/ptmp/jinma/handover/HANDOVER_2026-11-15_clean.md, source of origin/clean-cc (rad_m1_clean.hpp, rad_m1.cpp).

## Short answer
- The switch was NOT forced by a broad measured win of the cell layout. Your own 1-D prototype showed staggered is MORE
  accurate in smooth fields. The C++ cell layout was behind staggered on most IDORT decks when it was made the only
  layout (user call, 11-15).
- The decisive argument was **positivity/realizability** (proto2, 1-D). A second reason was **solver robustness at tight
  tolerance** (C++ staggered Newton stagnation, "it -1"). The rest came from the USER on 11-15: "cell central is more
  promising, continue on that".
- After that, every later fix went into the cell layout only. That covers the consistency bracket, AP dissipation,
  h-free Rhie-Chow, the WB limiter, the Ma2025 v/c sources, sp/cs geometry and the FV rays. So today's gap is mostly
  "never tried on staggered", not "proven impossible on staggered".
- **Important for your speed idea:** the C++ staggered layout ALREADY condenses F exactly into a scalar E system, the
  same way production does. See Q4.

## Q1. Measured evidence
### (a) proto2 1-D prototype, viper 11-14 (stag.py vs cell A3; same rays, SDIRK2, limiter chain)
| test | staggered | cell |
|---|---|---|
| thin wave chi dx 0.01, L1 N64/128/256 | 4.31e-4 9.82e-5 2.25e-5 (2.13) | 8.21e-4 2.04e-4 5.07e-5 (2.01) |
| thin wave chi dx 0.1 | 3.84e-5 9.08e-6 2.21e-6 | 9.69e-5 2.51e-5 6.27e-6 |
| thick: orders, kappa_num | same (1.94-2.02, +0.029..+0.002) | same |
| time order, thick Gaussian, c dt/dx 384..24 | 2.02 2.01 2.00 | 1.99 1.79 **0.85** |
| radiating wave ka 1, N32/64/128 | 9.6e-3 1.8e-3 4.1e-4 | 2.4e-2 6.0e-3 1.5e-3 |
| conservation | round-off | round-off |
| **realizability sweep, 176 cases** | **88 violate** (62 with minE + upwind-face test); worst minE -0.021 | 3, all round-off (1e-12) |
| limiter work | BE in 64 % of steps; all-low steps still bad 3883/6644 | BE 24 %, still bad 0 |
| beam into vacuum front width (cfl 0.5) | 29 cells | 8 cells |
| beam onto absorber chi dx 10, leak into first slab cell | 2.6e-3 | 5.9e-7 |

The prototype's own recommendation: keep cell. Staggered wins on smooth accuracy, but "cannot be made realizable with the
tools tested": there is no co-located realizability test that accepts a transported front, and its low-order fallback
is not realizable.

### (b) C++ at the switch, clean-cc dc67dba7b era (wt_cc_PROGRESS.md)
- rw space, ray closure, amp 1e-6: equal or staggered better. P 3.3 tau 10: staggered 6.29e-4/1.64e-4/4.25e-5, cell
  1.01e-3/2.67e-4/6.25e-5. P 1/3 tau 10: 1.47e-3 vs 1.90e-3. tau 0.1 and tau 1000: identical.
- rw time, iterate closure: identical, 2.00/2.00 in both layouts (5.80e-7 vs 5.81e-7). extrap: identical orders.
- **CPU, same formal-solution count, rw nt 4096: staggered is 1.7-2x cheaper.** extrap 450/758 s vs cell 893/1315 s;
  extrap2 777/1281 vs 1436/1901 s.
- IDORT before the bracket: staggered better on shd/shdg x0.25 (~0.02 vs 0.264/0.256), hcyl (0.025 vs 0.047) and
  nltek1 x0.5 (1.13e-2 vs 2.24e-2). Staggered worse on beams cb/cbg/ws66 x0.5 (L23 1.47/1.40/1.26 vs 1.68/1.62/1.35).
  nltek1 x1 vs exact: tied (9.5e-5 vs 9.4e-5).
- Cell after the consistency bracket (A5/A6, 11-15): cb 0.460, cbg 1.308, ws66 1.167, hcyl 0.019, shd 0.016-0.023,
  nltek1 9.1e-5. That is better than staggered on beams and equal on shd. **Staggered never got the bracket.**
- Conservation, periodic rw: cell energy 1.5e-16..1.7e-15, staggered 8e-15..3e-12. The momentum numbers from
  conserve.py are an artefact for both layouts (slice_x2 bug).
- Robustness: with D = I/3 the staggered Newton loop stagnates at clean_tol_nl 1e-13 ("it -1": BiCGStab 0 its). On
  N128 P 1/3 tau 10 that gives 38/606 failed solves, then BE, then dt_retry. Cell had 0, but cell used the effective
  max(tol_nl, 10 tol) floor that staggered did not have. Orion pkg_tol at 1e-13/1e-12: cell rows clean, staggered still
  failed massively. fixedD N4096 staggered rows: ~1230 BE failures each at tol_nl 1e-12. **I read this as a missing
  tolerance floor / round-off exit in the staggered path, not an inherent defect.** It was never retested with the
  conv round-off exit (b44bd47cd), because the user froze staggered first.
- **Not measured for staggered at all:** WB envelopes, L_z, multi-block or cubed-sphere seams, thin sphere / feed on cs,
  he_star, C++ realizability counts on the beam decks. There are no staggered numbers to compare on any of these.

## Q2. Inherent to F on faces vs implementation gaps
Inherent, or at least not solved in the prototype:
1. **Realizability has no natural test.** E and F are not co-located. At a front, (F_L+F_R)/2 > cE_i in the filling
   cell, even though the transport is fine. In thick slabs, interpolating F across an e^-10 decay overshoots cE.
2. **The low-order staggered flux is not characteristic upwind.** E +- F/c are never upwinded, so even all-low-order is
   not positive. The cell scheme's all-low BE is realizable, and that is our safety net.
3. **Transverse F and cell F have to be interpolated** from faces wherever co-located F is needed: the ray closure
   input, v.F and the beta terms, gas momentum, realizability, the bracket's cell-face bookkeeping, and the cs seam frame
   rotation. Each interpolation adds the front and slab artefacts of item 1.

Implementation gaps (likely fixable):
- The Newton stagnation at tight tol (needs the cell-style floor or the conv exit).
- The missing consistency bracket, AP dissipation and the Ma2025 v/c terms on staggered.
- No sp/cs or multi-block seam port tested for staggered.

## Q3. Staggered with bracket + ray face fluxes + AP dissipation?
- **Never run.** clean_consistency is refused unless clean_layout = cell (rad_m1.cpp:114-133). The ray face fluxes
  (CLN_F1/F2/FX/F3/FH) and the CLN_G curvature fade were written for the cell layout.
- Expectation: on smooth and diffusive decks (rw, nltek, shd, hcyl), staggered + bracket should reach or beat the cell
  layout. It was already at the cell's post-bracket level on shd/hcyl without the bracket, and smooth accuracy is
  2-4x better. A ray face flux lands naturally on a face unknown.
- On beams / 2-D transients I do not expect it to match. The bracket fixed cell accuracy, not realizability, and
  realizability on the cell side still depends on the cell all-low-order fallback. Staggered has no equivalent of that
  fallback (Q2 items 1-2). The best guess is beam accuracy comparable to the cell layout, with more limiter/BE work and
  no positivity guarantee.

## Q4. Scalar-E Schur complement
- **This is already what the C++ staggered path does.** rad_m1_clean.hpp header: "Each implicit stage solves for E in
  every cell", with
  F_f = (1-w)[theta (G_f - h c^2 grad(DE)) + (1-theta) A_f] + w c (h+_f E_L + h-_f E_R), theta = 1/(1+h c chi_f),
  and the cell row E(1 + h c kE) + h div F = E_old + h c kP a T^4. F is recovered afterwards (CleanMomentumChange,
  face averages into u0). The solver is the x1-line tridiagonal-preconditioned BiCGStab on one scalar per cell. That is
  the production structure (theta elimination; compact 3-point diffusion in the thick limit, see the proto2 README "RC
  check").
- Note that proto2 stag.py is a DIFFERENT staggered scheme from the C++ one. It has the cell "F flux" Psi_i = c^2 D_i E_i
  - alpha_i (a_i/2)[F]_i, which couples neighbouring face F's, so it does NOT condense exactly. The realizability
  numbers in Q1(a) belong to that scheme, not to the condensable C++ one.
- What breaks exact condensation (keep each lagged in the Picard loop, as the C++ staggered path already does for velocity):
  1. Any F-jump dissipation ([F] at cells) or upwinding of E +- F/c. That is exactly what positivity wants, so speed
     and realizability pull in opposite directions here.
  2. Implicit gas velocity coupling through F (v.F, Ma2025 beta terms, CleanCellDrag). Lag it, or keep it in the gas
     fixed point.
  3. The consistency bracket and PLM weights are fine if they stay explicit/lagged (RHS-only).
  4. Multi-D D_off: grad(D E) with D_12 is still E-only (a 9-point stencil), so it condenses, but the operator is no
     longer an M-matrix (design doc 2.2: split D_diag / D_off, limit D_off).
  5. cs seams: the face-normal F on a seam needs one frame; the transverse F for the closure must be rebuilt
     consistently on both sides.
- Cost side: the cell layout has 4 unknowns per cell in 4x4 blocks, with y coupled through (a/2)[U] and RC, so it
  cannot condense. A scalar E system is 4x fewer unknowns with 1x1 blocks. Even with no algorithmic change, the C++
  staggered was 1.7-2x cheaper per step on CPU rw.

## Suggested way to evaluate (if you and the USER go ahead)
Hybrid: the condensable C++ staggered E-system as the high-order stage solve, and the cell-layout all-low-order BE as
the realizability fallback. In beam decks every step already falls back to BE in both layouts, so the fallback matters
more than the stage layout there. Port the bracket (ray face fluxes onto face F) and the tolerance floor first, then
rerun: rw (prod tol), nltek1, shd/shdg, hcyl, cb/cbg/ws66 (scores + neg/nreal counts), and the R136 first-stage solve
time on MI300A against cell + x1line.

-- Orion (CPU help only; clean-cc ownership unchanged: viper)
