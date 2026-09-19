---
name: cs-raisevel-missing-floor
description: The cubed-sphere dhj NaN was GnomonicEquiangleRaiseVel re-deriving the internal energy WITHOUT re-applying the floors, feeding a non-positive eint to the tabulated EOS; fixed, spherical polar bitwise unchanged
metadata:
  type: project
---

**THE CAUSE OF THE CUBED-SPHERE HOT-JUPITER NaN** (found 2026-09-01, see
[[cs-dhj-long-run]] for the 2x2 that isolated it to the grid).

`Coordinates::GnomonicEquiangleRaiseVel` (and its MHD twin) runs immediately after
`ConsToPrim` on the cubed sphere and **overwrites `w0(IEN)`** with an internal energy
re-derived using the non-orthogonal metric:

    eint = u0(IEN) - 0.5*(m1*v1 + m2*v2 + m3*v3)

It then handed that value straight to `eos_.TemperaturePressureGamma1`.  **No floor, and
no positivity guard.**  `SingleC2P_GeneralHyd` guards the very same call with
`if (e_positive)` and follows it with the pressure and temperature floors -- precisely
because the tabulated inversion is undefined below the table.  The routine's own doc
comment claimed "the floors ... are left exactly as ConsToPrim left them", which is the
one thing it did NOT do: ConsToPrim floored an energy computed with an ORTHONORMAL kinetic
energy, and the metric cross term MOVES that energy, so the floor has to be retested.

**Signature, and how it was found.**  One cell goes NaN, then a symmetric +/-2 STENCIL
CROSS appears around it (not a filled cube) -- the fingerprint of a single bad POINTWISE
value spread by axis-aligned reconstruction.  Whole grid NaN ~100 cycles later.  Localised
by bisecting `tlim` on a restart (a dump is always written at tlim, which is the way to get
dense output from a restart -- overriding `output*/dt` on a restart does NOT work, the
schedule comes from the restart file).

**Refuted along the way, all cleanly:** the panel seam (layer profiles across every seam
are smooth, L0/L4 = 0.75-1.12, NaN sites not outliers); the CFL / non-orthogonal timestep
(CFL/4 still dies identically); global mass conservation (-0.63%, and stable);
`etotgrav` (its routines are geometry-agnostic, phi depends on r only); the density drain
itself (**spherical polar does the SAME drain to a LOWER minimum, 6.4e-14 vs 9.9e-14, and
survives** -- it is a normal shared transient, not the cause).

**Fix** (uncommitted at time of writing, on `cs-wire-wip`): re-apply the floors to the
corrected energy, mirroring `SingleC2P_GeneralHyd` step for step, and update `u0(IEN)` in
step so the conserved energy cannot keep sinking and re-trip the floor every cycle.
`u0` is now non-const in both signatures.  **The PRESSURE floor alone is not enough under
a tabulated EOS**: at upper-atmosphere densities e(d,pfloor) sits at T ~ 0.05 K, far below
the table's 31.6 K floor, so it is the TEMPERATURE floor that keeps the lookup in range.

**Evidence.**  Identical restart at t=1.90e4: unfixed -> first NaN t=1.930e4, 100% NaN by
t=2.05e4, dt frozen.  Fixed -> clean to t=3.0e4, zero non-finite cells, zero NaN history
rows, min dens exactly on dfloor.  Regression, all clean:
  * **spherical polar BITWISE unchanged** -- all 44 dumps identical old vs new binary;
  * **every ideal-EOS cubed-sphere test BITWISE unchanged** -- cubed_sphere_mhd,
    _rigidrot and _smr histories identical, _blast and _resist gate output identical.

**Multi-rotation confirmation.**  Cubed-sphere hydro FROM SCRATCH now runs 25100 cycles
to t = 5.81e5 = **1.91 rotations** with zero NaN rows and a healthy, steady dt (22-24).
The unfixed binary died at t = 1.93e4 = 0.063 rotations, so this is 30x past the failure.
Mass drift settles at -0.65% and is FLAT (-0.630% at 0.3 rot -> -0.647% at 1.91 rot);
energy -0.97% and creeping, which is the floors plus radiative loss, not a runaway.

**Two rounding traps hit while proving that regression, both worth remembering.**
(1) Folding the MHD `(E - KE) - ME` into `E - (KE + ME)` RE-ASSOCIATES the arithmetic and
moved every symmetric-zero history column in the last bits.  Keep the subtractions
separate and in the original order.  (2) Even merely HOISTING `const Real e_k = ...` above
the subtraction perturbed the result, by changing FMA contraction -- compute it only
inside the branch that writes it back.  Mass, energy and the non-zero columns were exactly
identical throughout; only quantities that are ZERO BY SYMMETRY moved, at 1e-17 and below,
which is how to tell round-off from a behavioural change.
**And: `-d <dir>` APPENDS to an existing .hst.**  Re-running an A/B into dirs left over
from a previous build compares accumulated output from BOTH binaries and reports a
spurious DIFFERS.  Clear the directories.  [[validate-the-instrument]].

This is also very likely the "negative internal energy at panel corners" left open in
[[cubed-sphere-mhd-seam]] -- the ideal-EOS path had the same missing floor, it just
degrades silently instead of returning NaN.
