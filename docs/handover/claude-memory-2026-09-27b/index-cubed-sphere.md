---
name: index-cubed-sphere
description: Memory links for the cubed-sphere grid: seams, MHD, validation, GPU/MPI history.
metadata:
  type: reference
---

## Cubed sphere — seam halo (from RULE ZERO block)

- **[cs SEAM HALO WRONG WITH >1 BLOCK/PANEL (cc fixed 94c7165d on cs-seam-4x4; 2x2 was 11x off too); bvals_fc (MHD) STILL HAS BOTH DEFECTS](cs-seam-multiblock-halo-fix.md)**

## Cubed sphere — start here

- **[RT SEMI-IMPLICIT source](rt-semi-implicit-changes-dhj-answer.md) — 048dff30; explicit overcools, semi-implicit is right**
- **[IN FLIGHT, STOPPED 09-11](inflight-2026-09-09-viper.md) — cs collapse closed; c26b01ac pushed**
- **[GPU BUILD RACE: flag never compiled](gpu-build-race-flag-not-compiled.md) — verify flags with `strings` plus a startup print**
- **[cs RADIAL = sp RADIAL](cs-radial-unification.md) — 979edada; earlier cs baselines stale**
- **[cs STRETCHED-grid SOURCE TERM off 0.57-2.15x](cs-stretched-source-term-bug.md) — fixed 13a97399, index-space dr**
- **[cs STRETCHED-grid RESISTIVITY was ANTI-DIFFUSIVE](cs-stretched-resistive-rcm-bug.md) — fixed c5c85e3b, r_cm missed the stretch**
- **[cs seam CO-LOCATION found and fixed](cs-seam-colocation-fixed.md) — the clamped-window cubic, 5-11x**
- **[cs seam order: 2nd order confirmed](cs-seam-order-limiter.md) — global L1 2.01, done**
- **[cs seam de-staggering fix REFUTED](cs-seam-destag-refuted.md) — superseded**
- **[resistive seam 1st order: closed](cs-resistive-seam-order.md) — operator innocent, halo inputs guilty**
- **[cs cross-level SEAM halo: closed](cs-crosslevel-seam-halo-first-order.md) — a01ace75 + 180a9b3e**
- **[cs CUBE-VERTEX corner halo x radial ghost: fixed](cs-cube-vertex-corner-radial-ghost.md) — 397b4ad3**
- **[cs seam conservation closed](cubed-sphere-seam-conservation.md) — 985faa22, exact**
- **[cs seam EMF closed](cubed-sphere-seam-emf.md) — 42323a66; div B is the wrong gate**
- **[cs RESISTIVITY](cubed-sphere-resistivity.md) — 4bfacdd8; two-pass curl, 2nd order**
- **[cs SMR: refined MHD converges, MPI-clean](cubed-sphere-smr.md) — 74cbc8df; rank dependence retracted**
- **[cs MHD convergence](cubed-sphere-mhd-convergence.md) — e63b571a; residual is radial-BC phase lag**
- **[SHOCKS through cs seams](cs-shocks-through-seams.md) — resample clamped; FOFC not needed**
- **[cs blast vs a CARTESIAN grid](cs-blast-vs-cartesian.md) — seam costs nothing, vertex ~2x worse**

## Cubed sphere — validation and limits

- **[cs_test FF decay + rot_axis](cs-test-ffdecay-rotaxis.md) — 4990eb41; no 1.5-order region**
- **[cs PURE HYDRO validation](cs-hydro-validation.md) — machine precision under a shock; L1 2.3-2.4**
- **[cs MHD + RESISTIVE validation by region](cs-mhd-validation.md) — 2nd order; energy drift is physical Ohmic heating**
- **[cs ANGULAR MOMENTUM: sp exact, cs cannot be](cs-angular-momentum.md) — 6.9e-4/rot at nx2=32**
- **[cs: radiation + srcterms now REFUSED](cs-unsupported-physics-guards.md) — b283ad3a**
- **[cs NARROW-BLOCK resample degeneracy](cs-narrow-block-resample-degeneracy.md) — stencil inverts below 3 cells**
- **[cs GENERAL/TABULATED EOS stale cache: fixed](cs-general-eos-stale-cache.md) — cached pre-correction, 6-15% wrong**
- **[cs MHD blast from a VECTOR POTENTIAL](cs-mhd-blast.md) — b283ad3a; gate with two numbers**
- **[cube-vertex corner fill promoted, default ON](cs-wire-fill-wip.md) — `<mesh>/cs_vertex_fill`**
- **[cs cube-vertex REAL fill prototyped, not built](cs-cube-vertex-real-exchange.md) — sampling beats extrapolation 3.4-5.5x**
- [cs determinism test on orion](cs-determinism-test-orion.md) — DONE 193377: CPU restarts are BIT-IDENTICAL, even 16x7 vs 28x4; the viper 5e-6 divergence is GPU-specific or a binary difference
- [cs orthogonal-KE AUDIT: red_giant ghost fills + IC + inner passes FIXED in the working tree 09-09 (T14: correct, not the killer); latent sites listed in docs/handover/NOTE-2026-09-09-cs-covariant-basis.md (pushed 40aef0ad)](cs-orthogonal-ke-audit-2026-09-09.md)

## Cubed sphere — history

- **[Cubed sphere ON GPU](gpu-this-capture-device-lambda.md) — af539941; constant-memory closure limit**
- **[cs GPU multi-block fault: fixed](cs-gpu-multiblock-fault.md) — 4302a008**
- **[Cubed sphere + MPI: fixed](cubed-sphere-mpi-hang.md) — 1c2e29d6, unmatched receive at a vertex**
- **[cs MHD seam fixed](cubed-sphere-mhd-seam.md) — 4f19a244; transform always exact**
- **[Cubed sphere: MHD](cubed-sphere-mhd.md) — e0c74357; dxedge was zero**
- **[Cubed sphere committed](cubed-sphere-committed.md) — 9492a946**
- **[cs along-seam resample, 2nd order](cubed-sphere-seam-interp.md) — flat-residual claim retracted**
- **[cs seam basis transform fixed](cubed-sphere-seam-basis.md) — tangent-basis transform**
- [Cubed sphere panel frames](cubed-sphere-panel-frames.md) — duplicate copy had 3/4 swapped
- [Cubed-sphere hydro fix](cubed-sphere-hydro-fix.md) — halo axis bug fixed
- [Cubed sphere: x1 as radial](cubed-sphere-x1-radial.md) — done
- [Cubed sphere: hydro state](cubed-sphere-hydro-state.md) — 4 bugs fixed; rigid-rotation part retracted
- [Cubed sphere for the hot Jupiter](cubed-sphere-for-hot-jupiter.md) — dt argument dead
