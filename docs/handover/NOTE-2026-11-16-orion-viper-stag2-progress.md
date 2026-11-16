# NOTE 2026-11-16 Orion -> viper: stag2 (local staggered study) progress report, rounds 1-13

FYI only; no action needed. clean-cc is untouched and Orion pushed nothing there.
- Local branch clean-cc-stag2: 23 commits on d6b6ed4e2. It rebases cleanly onto push 15 e3f4881dc.
- With every new key off, the full reg harness is bitwise to d6b6ed4e2 after every round (last check: job 230649).
- Write-ups (Orion): /orion/ptmp/jinma/clean_1115/stag2/cxx/RESULTS.md .. RESULTS13.md

## Recommended staggered key set
clean_layout = staggered + clean_stag_split (default on for staggered) + clean_stag_fsolve + clean_stag_fkin + clean_stag_fwork
+ clean_consistency + clean_rays = fv1.

## What the set does
1. **Solve and energy.** Every solve is scalar-E: E is solved with the energy flux Phi, where Phi = the theta-eliminated high-order flux
   blended by w(C1) with the rays' kinetic half-range flux. The consistency bracket is applied on Phi.
2. **Momentum.** A separate face-momentum F solve, banded per direction, runs inside the Newton loop. It has cc-style LF dissipation, and
   in streaming cells the rays' kinetic momentum flux, including (D - s^2) min(J_ray, E). It is in flux form, so momentum and L_z are
   conserved to round-off.
3. **Gas.** The gas gets only physical quantities: the force h chi (F - A)/c and the exchange on the solution. The linear-solve residual
   stays in E.

## Results vs the cell layout (d6b6ed4e2), production tolerance, 0 CLEAN WARNING
- **Cost per cycle.** Staggered is 3-15x cheaper on cb/shd/hcyl (about 1 solve per stage, 1 formal solution per cycle) and 2-4x cheaper
  on sp. The extra F solve costs +12%.
- **Beams with chain off.** cb 0.460 / FWHM 1.0, cbg 1.100, ws66 0.921, shd 0.0118, shdg 0.0295. That equals cell.
- **Beams with the positivity chain on.** The bracket is stripped: cb 1.448, FWHM about 11.
- **Negativity with chain off.** Final states are about as negative as cell's (cb min E/max E -7.9e-7 vs -8.9e-7; shd -5.9e-5 vs
  -8.7e-5). A user decision is pending.
- **Smooth rw.** tau >= 10: orders 2.0, errors 0.8-0.9x cell. Thin rows need clean_tol about 1e-2 x the cell rule (|B| is dominated
  by the old-flux divergence). One point remains: P 1/3, tau 0.01, N512 is still at order 0.98.
- **nltek1.** 9.105e-5 (cell 9.10e-5).
- **Spherical polar.** 1-D, 2-D and 3-D. The uniform isotropic state is exact; closed-box energy and L_z drift are 1e-15.
  - flz q: 0 warnings, t 32.88 in 25 cycles (cell 41), wall v_r inward like cell.
  - flz Q = cell.
  - sp3 + force_reference = clean, 1 vs 2 blocks: Er 1.1e-9 (cell 2.4e-9).
- **Moving medium (Ma2025 terms ported).** b^4 behaviour like cell (-6.65e-9 at beta 0.01).
- **Multi-MeshBlock and restarts.** Decomposition-invariant to 1e-9..1e-11; restarts bitwise.

## Findings that may matter for the cell layout
- Cell shdfb: the Newton de denominator (NOTE 3a5363e6b).
- In 1-D proto2 the cell outflow-copy boundary lets energy in from vacuum at cfl 1000.
- The F solve had a ghost-slot convergence bug. That was staggered-only, but check any line solve whose test includes the line's
  last face slot.

## Missing
Cubed sphere (planned, ~250-350 lines), the r-face kinetic transverse stress on sp, and the bracket positivity decision.
