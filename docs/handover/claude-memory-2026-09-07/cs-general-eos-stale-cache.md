---
name: cs-general-eos-stale-cache
description: FIXED -- on the cubed sphere a GENERAL/tabulated EOS cached p, Gamma_1 and T from BEFORE the gnomonic metric correction, leaving them 6.3% wrong at the midlines and 15.4% at panel corners. Prerequisite for any cubed-sphere hot-Jupiter port
metadata:
  type: project
---

## The bug

On the cubed sphere the order of operations was:

1. `peos->ConsToPrim()` (`general_{mhd,hyd}.cpp`) root-finds T and p from the conserved
   state, subtracting an **orthonormal** kinetic energy `0.5*(m1^2+m2^2+m3^2)/d` -- wrong
   on the non-orthogonal gnomonic basis by the metric cross term -- and for MHD a magnetic
   energy built from the NON-ORTHOGONAL face-normal triple.  It caches the results in
   `wder(IDPR)`, `wder(IDG1)`, `wtemp`.
2. `Coordinates::GnomonicEquiangleRaiseVel{,MHD}()` then corrects `w0(IEN)` with the
   metric and rebuilds `bcc` in the orthonormal frame.
3. **Nothing refreshed the cache.**  `wder`/`wtemp` are written NOWHERE ELSE.

They are read by `mhd_fluxes`, `mhd_newdt`, `hydro_fofc`, `prolong_prims`, and by the
gnomonic geometric source term itself (`coordinates.cpp` reads `wder(IDPR)` when `gen_`).

## Measured -- much larger than the "few percent" I guessed

New gate `CSTestEosCacheCheck` recomputes p and T from the state `w0` actually holds and
compares with the cache, split by `|cos_cell|` so the corner cells cannot be diluted:

```
  WITHOUT the fix   midline |cos|<=0.2  dp = 6.330e-02   corner |cos|>0.2  dp = 1.543e-01
  WITH the fix                          dp = 0.000e+00                     dp = 0.000e+00
```

**6.3% at the midlines and 15.4% at the panel corners**, exactly zero after.  The corner
excess confirms the O(cos_cell) scaling; the midline floor is the MAGNETIC energy rebuild,
which is large here because the test field is O(1).

## The fix

`GnomonicEquiangleRaiseVel{,MHD}` now take the `EOS_Data` plus `wder`/`wtemp` (the same
argument-passing precedent `SrcTermsGnomonicEquiangle` already uses, because a two-fluid
run has both modules' arrays) and, when `eos_data.IsGeneral()`, re-solve via
`eos.TemperaturePressureGamma1(d, eint_corrected, tguess = wtemp, ...)` -- warm-started on
the temperature ConsToPrim already found, so the root find is cheap.

**IDEAL EOS IS BITWISE UNCHANGED** (guarded by `IsGeneral()`), verified against a HEAD
build; MPI np=1/2 agree; HIP compiles.

## THE GATE NEEDS A FLOW -- it is vacuous on the static problem

The discrepancy is O(cos_cell) x tangential KINETIC energy, so on `iprob = 11` (v = 0) it
is identically zero and the gate proves nothing.  It is called from `CSTestConvErrors`
(iprob = 9, rigid rotation, tangential v is O(1)) as well as the resist check.  Test
config: `cubed_sphere_mhd_conv.athinput` with `eos = general`, `general_eos = gamma`
(no table file needed) plus a `<units>` block -- a general EOS refuses to run without one.

## GPU TRAP HIT ON THE WAY (and it was MY earlier commit)

`static const Real x = getenv(...)` at function scope is a HOST variable; a device lambda
referencing one fails to compile on HIP with "reference to __host__ variable in
__host__ __device__ function".  The shear-correction commit introduced exactly this and I
missed it by only CPU-building.  Read the env into a file-scope static, then take a LOCAL
`const` copy, which is captured BY VALUE.  Same trap as f493c209.
**GPU-COMPILE ANY COMMIT THAT TOUCHES A DEVICE LAMBDA.**
