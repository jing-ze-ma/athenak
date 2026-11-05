# TASK 2026-11-05 Caltech -> ORION (via viper if needed), USER: does Arepo IDORT make the Athena++ Doppler-velocity mistake?

**User:** "we are using a similar scheme in arepo and I want to check if we are making that same mistake there."

## The Athena++ issue (found 11-04, details: NOTE-2026-11-05-caltech-athpp-wall-collapse-cause.md, docs/handover/athpp-toy-1104/REPORT.md)
- Athena++'s implicit DO radiation (`GetTgasVel`, src/nr_radiation/implicit) hands the implicit radiation solve the velocity v* of the
  HYDRO-UPDATED state of the stage, i.e. after the hydro/gravity update but BEFORE the same step's radiation momentum source.
- Where radiation supports the gas (optically thick, near-Eddington: Gamma ~ 0.3-1), the hydro update pulls the gas by dt*a_rad each step and
  the radiation kick gives it back. v* is then off by dt*a_rad (BSG base: v* = -2.95 km/s vs the true end-of-step -0.63 km/s).
- The O(v/c) Doppler/advective flux (4/3) v* E/c then exceeds the true diffusive flux by Lambda ~ (4/3) kappa E dt / c (~65 at the BSG base).
  At an open/bath boundary this drains energy (the paper's wall collapse and deflation); at a closed wall it drives a cell-scale hot/cold instability.
- Signatures: the error is proportional to dt; it vanishes if the solve uses the start-of-step velocity v^n or a velocity that is
  self-consistent with the radiation force (AthenaK M1 does that); frozen-gas tests are clean; angle count / face scheme / resolution don't matter.

## What to check in Arepo (IDORT: src/modules/idort/, esp. idort_coupling.h, the implicit solve and its velocity / comoving-frame terms)
1. Which velocity the implicit radiation solve uses in the O(v/c) terms (comoving transform, Doppler/aberration, advective flux, work term):
   v^n, the post-hydro/gravity pre-radiation-kick v*, or a velocity updated with the radiation force inside the solve?
2. The order of operations per step (and per sub-step/time-bin with local time-stepping): hydro + gravity kicks, then the radiation solve,
   then the radiation momentum kick? If v* is used, estimate Lambda = (4/3) kappa E_rad dt / c in your production stars (RSG g25 envelope,
   He giant, ...), at the deepest thick layers and the inner boundary. Lambda >~ 1 means the error matters.
3. Cheap test (as our toy): a 1-D (or thin-column) thick, radiation-supported static atmosphere with an open (bath) inner boundary:
   run at dt, dt/2, dt/4. If the boundary energy drain or the wall-cell temperature drift scales with dt, the scheme has the issue.
   Then rerun with v^n in the velocity terms: the drain should vanish.
4. Report what you find (code path, velocity used, Lambda estimate, test result) in a NOTE on this branch.

Context: Caltech's local Arepo2 copy (arepo2_gpu, 10-24) is being checked in parallel; Orion's current IDORT may differ, so please check your version.

## Addendum (Caltech code read of the local 10-24 copies arepo2 6cef4d1d6 / arepo2_opt c663a54c5 / arepo2_nltree 31a5b39ff)
- Default `IdortUseVelocity 3` (module_idort.h:108, module_idort.cc:68; fix 7a7dc2c61 on 10-16) uses v^n = `SphP.IdortVel0`, stored in
  onFindTimestepsBeforeGravity (module_idort.cc:1277-1280), for all O(v/c) terms (beta at module_idort.cc:1362-1369). So no Athena++ error by default.
  The work term v.dp uses the post-hydro momentum, which is the correct energy bookkeeping.
- Error still possible: (i) mode 3 FALLS BACK to P.Vel = v* when none of SELFGRAVITY / EXTERNALGRAVITY / EXACT_GRAVITY_FOR_PARTICLE_TYPE
  is compiled (module_idort.cc:1215-1229): suggest using IdortVel0 unconditionally (one line); (ii) param files with IdortUseVelocity 1 or 2
  (mode 2 = exactly the Athena++ case). Lambda in a deep RSG/He-giant envelope ~ 3e2 if it happens.
- Please confirm for the PRODUCTION runs (RSG g25 etc.): the Config.sh gravity flag and the IdortUseVelocity value, and that Orion's current
  IDORT still has mode 3 as above.
