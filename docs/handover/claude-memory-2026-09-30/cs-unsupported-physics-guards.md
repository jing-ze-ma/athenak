---
name: cs-unsupported-physics-guards
description: b283ad3a -- radiation and the geometry-dependent source terms are now a startup FATAL on the cubed sphere; records the two guard-writing traps that produce a gate which can never fire
metadata:
  type: project
---

**b283ad3a.** Cubed-sphere support was already a startup FATAL for fofc, viscosity,
conduction and AMR (9f79a031), but `<radiation>` and the source terms were **never
gated**, so they RAN TO COMPLETION and returned a wrong answer -- the same silent
failure mode the resistive curl had before 4bfacdd8.

`src/radiation/` contains not one mention of `use_cubed_sphere`, and not even the
`use_spherical_polar` branch other modules carry. It is not untested, it is **unwritten**.

REFUSED: `<radiation>`; `rel_cooling` (forms the Lorentz factor as a NORM
`sqrt(1+ux^2+...)`, which needs the metric on a non-orthonormal basis); `rad_beam`
(Cartesian positions + angular momentum components); `const_accel` along x2/x3
(those basis vectors are neither orthogonal nor of unit physical length).

STILL ALLOWED, because their expressions carry **no geometry** -- not because anyone
measured them: `const_accel` with dir=1 (RADIAL), `ism_cooling` (drag is -k*rho*v in
any basis), the EOS (pointwise).

**TWO TRAPS, both of which produce a guard that looks fine and can NEVER FIRE:**
1. The block is **`<hydro_srcterms>` / `<mhd_srcterms>`**, NOT `<hydro>` / `<mhd>`.
   Probing the latter is a dead guard, and a test that puts the parameter in the same
   wrong block AGREES WITH IT and reads as a pass.
2. These switches are read with **`GetOrAddBoolean`**, so they are legitimately
   PRESENT AND FALSE in an input that spells out what it is not using. Testing
   `DoesParameterExist` alone would **REFUSE such a run**. Test the VALUE.

Note `turb_driving` needs no guard: `<turb_driving>` is not in
`ParameterInput::CheckBlockNames`' whitelist, so that block is rejected on EVERY grid
before the gate runs -- a separate pre-existing oddity, given
`srcterms/turb_driver.cpp` reads exactly that name.

**How to apply:** verify a guard with NEGATIVE CONTROLS (the switch set false, the
"off" case) as well as positive ones. **Command-line overrides cannot create a block**
(`Block name 'hydro_srcterms' on command line not found`), so these must be exercised
with real input files -- a `-i ... hydro_srcterms/x=true` test dies identically in every
case and reads as a uniform pass. Verified over 8 cases on CPU, MPI (np=2) and GPU, all
three agreeing. See [[cubed-sphere-resistivity]], [[validate-the-instrument]].
