---
name: cs-dhj-long-run
description: RESOLVED -- the 2x2 {cubed sphere, spherical polar} x {MHD, hydro} says the dhj NaN is the GRID, not the physics - cs dies in 1/10 rotation with hydro alone, sp is clean for tens of rotations
metadata:
  type: project
---

The first long cubed-sphere hot-Jupiter run went NaN inside one rotation, silently, with a
temperature-floor storm.  A 2x2 of {cubed sphere, spherical polar} x {MHD, hydro} at
nx1 = 128 (jobs 11301436 / 11301510 / 11301511 / 11301512, launched 2026-08-31 23:55,
read 2026-09-01) **isolates the cause to the GRID**:

| | MHD (eta 1e12) | hydro |
| --- | --- | --- |
| **cubed sphere** | NaN at t = 3.05e4, the FIRST hst row after t=0 | NaN by t = 3.8e4 |
| **spherical polar** | clean to t = 1.0e7 (33 rot) | clean to t = 2.0e7 (66 rot) |

Rotation = 3.05e5 s, so the cubed sphere dies in **under 1/10 of a rotation** with the
magnetic field and the resistivity entirely absent.  The two inputs differ ONLY in
grid parameters (nx2/nx3, the x2/x3 BCs, use_cubed_sphere vs use_spherical_polar +
theta stretch) -- verified by diff.  Physics, EOS, RT and the radial stretch are identical.

**Why the runs reported success.**  `sacct` says COMPLETED, exit 0, and the cs hydro job
took 38 SECONDS.  Once the state is NaN the timestep is NaN too, `time` races to
`tlim = 8.64e7` in ~170 cycles and the code exits normally.  The cs MHD job instead froze
dt at exactly 1.04466e+02 and ground out 827k cycles over 8 h, every hst row NaN.
**A job that COMPLETES with exit 0 is not a job that ran** -- see [[validate-the-instrument]].

Dense-output reproduction is cheap: the failure is ~700 cycles / ~26 s of wall on 2 GPUs.
The t=0 binary dump is the only clean one in the archived runs -- the bin interval
(6.10e5) is far longer than the time to failure, so nothing useful was captured.

Related: [[dhj-cubed-sphere-port]] (the port itself works and tracks sp on a static setup),
[[cs-hydro-validation]], [[cs-angular-momentum]].

**RESOLVED, same day.**  The cause was `GnomonicEquiangleRaiseVel` re-deriving the
internal energy with the gnomonic metric and never re-applying the floors, so a
non-positive energy reached the tabulated EOS.  Full account, evidence and the list of
refuted hypotheses in [[cs-raisevel-missing-floor]] -- **go there**.  This file is kept
for the 2x2 itself, which is the measurement that pointed at the grid.
