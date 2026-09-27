---
name: cubed-sphere-mhd-convergence
description: The cubed-sphere MHD convergence test (cs_test iprob=9) -- the flat error was FOUND and FIXED (corner EMF, e63b571a); the residual is an O(dt) rotation phase lag GENERATED AT THE RADIAL BOUNDARY, now reduced to one law delta = -0.906*omega*(frac+0.41)*dt, with only the offset -0.41 unexplained
metadata:
  type: project
---

Opened 2026-08-28 (WIP fe39529d), worked 2026-08-28/29. The non-converging error is
found, explained and **FIXED in e63b571a**. The smaller residual left behind is traced to
the radial boundary (12bd3c67, 6ebef4f3); one narrow question about it remains.
**HEAD is 282d6e31, `src/` clean, everything committed.** Follows [[cubed-sphere-seam-emf]] and [[cubed-sphere-mpi-hang]].

## THE BUG: `mhd_corner_e.cpp` had NO cubed-sphere handling

GS07 builds each edge EMF from the four face EMFs **plus an upwind correction made from
CELL-CENTRED EMFs**, and `e_cc_3d` forms those with the flat-space cross product
`v x B`. On the cubed sphere that is invalid two ways at once: the gnomonic tangent basis
is non-orthogonal (so a cross product needs the metric), and the two operands are in
**different frames** -- `w0`'s velocity is CONTRAVARIANT (`GnomonicEquiangleRaiseVel`)
while `bcc0` is in the ORTHONORMAL frame (`GnomonicEquiangleRaiseVelMHD`). The face EMFs
were fine; `GnomonicEquiangleEmfX1` rotates them. Only the correction was wrong, and it is
wrong by O(cos_cell) -- **O(1), not O(h)** -- so the scheme was INCONSISTENT, not merely
inaccurate.

**Fix:** a `use_cubed_sphere` branch in `CornerE`'s 3D path using the plain four-face
average (Balsara-Spicer), which is frame-correct. The proper fix is a gnomonic
cell-centred EMF; the branch says so.

## The measurements (all `cs_test` iprob=9, `inputs/tests/cubed_sphere_mhd_conv.athinput`)

Everything below was needed; **the ruling-out is what made the diagnosis safe**:

* **iprob=3, no field** -- CONVERGES (L1(v) 4.40e-4 / 1.02e-4). Hydro is clean.
* **omega=0** (field, no flow) -- CONVERGES. Needs BOTH B and v, i.e. `v x B`.
* **cfl halved** -- L1(B) changes 0.3%. **Purely spatial**; kills the radial-BC
  substage-time candidate and any temporal-order story.
* **b0c=0.01 (passive field)** -- hydro returns EXACTLY to the iprob=3 numbers while
  `L1(B)/b0c` stays flat. Isolates it to the induction; the back-reaction was a symptom.
* **Best-fit uniform B** -- only ~10% of the error is a systematic rotation-angle error;
  90% is structured. (The fit recovers B0 to 1e-13 at t=0, which validates it.)
* **`tlim=0.05` then `nlim=1`** -- the killer. Fast waves move <1 cell, so the deep
  interior cannot be seam error advected in. Local truncation error per unit time in the
  deepest bin: **2.691e-3 / 2.371e-3 / 2.302e-3** at nx2 = 16/32/64. FLAT = zeroth order.

**Result of the fix**, passive field, tlim=0.05: L1(B) 2.350e-5 / 6.415e-6 / 1.777e-6,
ratios **3.66, 3.61 -- second order**, and 124x smaller at nx2=32. Linf stops growing.
At tlim=1, b0c=1, refining all three directions: 5.535e-4 / 2.092e-4 / 9.357e-5, ratios
2.65 / 2.24 -- converging (it was FLAT at 6.2e-3), but only ~1.2 order.

## The residual ~1.2 order at t=1: TRACED TO THE RADIAL BOUNDARY

Committed as **12bd3c67** (diagnostics only; `problem/bc_time_frac` defaults to 0 and
reproduces the previous run BITWISE).

**CAREFUL, this was initially misread.** The radial BC IS the source (see ORIGIN below),
but the naive framing -- "the ghosts hold `pm->time` instead of the substage time, shift
the clock and it goes away" -- is WRONG, and shifting the clock makes it worse. The first
thing that misled me was the profile against distance from the RADIAL boundary: the
residual is **UNIFORM**, 2.12e-4 at the boundary vs 2.03e-4 mid-shell, and equally flat
against seam distance. That is NOT evidence against a boundary origin -- the error is
generated at the boundary and then spread through the volume, which the 1/L scaling below
proves directly. Do not re-run the profile expecting it to localise this.

**What it is:** the best-fit uniform field is a near-perfect RIGID ROTATION of the exact
one -- magnitude preserved to 2e-5 -- lagging in PHASE by 4.65e-4 rad, i.e. **0.23% of the
omega*t = 0.2 the field should have turned through**. So the scheme rotates B almost
exactly, just slightly behind.

**It is not the integrator.** At fixed grid the systematic part scales as dt: rk2 gives
1.471e-4 / 8.128e-5 / 4.821e-5 as the CFL is halved, and **rk3 gives 1.034e-4 / 5.927e-5 /
3.716e-5** -- still ~O(dt), though rk3's own truncation error is O(dt^3). Only a boundary
term can be O(dt) here, but see below.

**`problem/bc_time_frac`** scans the ghost time over `t^n + frac*dt`. If the ghosts simply
held the wrong time, `frac = 1` would fix it -- SSP-RK2 (Heun, confirmed in driver.cpp)
puts BOTH stage outputs at t^n + dt. **It does not:** L1(B) rises
2.09e-4 / 3.30e-4 / 4.58e-4 / 5.77e-4 / 6.77e-4 over frac = 0 .. 1, and at frac = 1 the
error stops being a pure rotation. **frac = 0 -- the existing code -- is the best of that
family.** Do not "fix" this by advancing the ghost time; it is 3.2x worse. (This half of
the scan looked monotonic and was read that way for a session; it is a V -- see **The
ghost clock, RESOLVED** below, which supersedes this paragraph.)

### What the lag IS -- the fingerprint (6ebef4f3)

Nine measurements, all at nx1=16, nx2=nx3=32, tlim=1 unless stated. 1-7 say what it is
NOT; 8-9 identify it:

1. **A rigid-rotation phase lag.** |Bfit|/|Bex| = 1.0000194, so the field keeps its
   magnitude; it is just turned 4.65e-4 rad short of where it should be.
2. **Accumulates linearly in time** -- delta/t = -3.48e-4 / -3.93e-4 / -4.29e-4 / -4.65e-4
   at t = 0.125 / 0.25 / 0.5 / 1. A constant RATE deficit: the field turns ~0.2% slower
   than the flow.
3. **Proportional to dt** -- systematic part 1.471e-4 / 8.128e-5 / 4.821e-5 as the CFL is
   halved; `A + K*dt` fits it to <1%, with A = 1.51e-5 (the genuine spatial part).
4. **NOT the flow.** `omega_fit/omega - 1` <= 1e-4, accounts for only 2.9e-5 rad of the
   4.65e-4, and has the wrong sign at two of three timesteps. In pure hydro it is
   1.000033350 vs 1.000033353 at CFL 0.3 vs 0.075 -- dt-independent to EIGHT DIGITS.
5. **NOT the time-integrator order.** delta/dt = -0.075 (rk2), -0.052 (rk3), -0.078 (rk1).
   No explicit RK method has an O(dt) phase error on a linear rotation anyway.
6. **NOT the reconstruction order.** Donor cell gives -1.4598e-4 against PLM's -1.4702e-4
   at CFL 0.075 -- the same to 1%. First and second order spatial reconstruction produce
   the SAME lag.
7. **NOT the grid.** At essentially fixed dt (6.51e-3 vs 6.20e-3), doubling the resolution
   from nx2=16 to 32 changes delta by only 18%.

### ORIGIN: the RADIAL BOUNDARY, diluted over the shell as 1/L

Two further measurements identify it, and the second is quantitative:

8. **The deficit is INDEPENDENT OF omega.** `deficit/dt` = 0.378 / 0.376 / 0.375 / 0.374
   at omega = 0.05 / 0.1 / 0.2 / 0.4 -- constant to 1% over an 8x range. So the field
   turns at **omega*(1 - 0.375*dt)**: the fractional rate deficit is a pure function of dt.
   This is exactly the signature of a ghost-zone PHASE error: the ghosts hold the exact
   field, which lags the interior by omega*dt, and dragging the interior back by a phase
   proportional to omega gives a FRACTIONAL rate deficit with the omega cancelled out.
9. **The coefficient scales as 1/L**, L = shell thickness. At FIXED dt (6.203e-3) and
   FIXED dr (0.0625), varying only the radial extent: **0.375 / 0.186 / 0.124** for
   r in [1,2] / [1,3] / [1,4], against 1/L = 0.375 / 0.1875 / 0.125. **Agreement to 1%.**
   A boundary-generated error diluted over the volume, and nothing else, does that.

Supporting: **no `dt` appears anywhere** in `mhd_fluxes.cpp`, `src/reconstruct/*.hpp` or
`hlld_mhd.hpp`, so for a given state the interior EMF cannot be dt-dependent at all; the
ghost zones are the only operator input that changes with dt. And the earlier "not h, not
rk1/2/3, not dc/plm" results all follow: a Dirichlet boundary error is not a truncation
error of any order.

### The ghost clock, RESOLVED to one law (282d6e31, 2026-08-29)

The "frac = 1 is worse" puzzle was an artefact of scanning only frac in [0,1]. Extending
to frac = -1 (nx1=16, nx2=nx3=32, tlim=1) gives, with the phase read off the best-fit
uniform B (+ = ahead of exact):

    frac    -1.0    -0.75   -0.5    -0.25    0      +0.5    +1.0
    L1(B)  6.73e-4 5.15e-4 3.61e-4 2.27e-4 2.09e-4 4.58e-4 6.77e-4
    phase  +2.8e-4 +2.4e-4 +0.7e-4 -1.8e-4 -4.7e-4 -9.9e-4 -1.06e-3

* **L1(B) is a V with its minimum at frac = 0** -- any offset opens a ghost/interior jump
  that adds non-rotational structure. `frac = 0` stays the default. The old
  "monotonic over [0,1]" reading is **RETRACTED** (in the source comment too).
* **The PHASE is linear in frac** over |frac| <~ 0.5 (slope -1.01/-1.13/-1.04e-3 per unit
  at the three intervals there), zero at **frac = -0.41**, slope **-0.906 per unit of
  ghost phase offset** omega*frac*dt. It flattens at |frac| >= 0.75, but there the error
  is no longer a rotation at all, so the extracted "phase" is not meaningful -- do not
  read the flattening as saturation of a physical response. Slope of order one = the shell's mean field LOCKS to
  the phase the boundary holds (Alfven crossing L/v_A = 1 = tlim) rather than accumulating
  from it. One law covers the family:  **delta = -0.906*omega*(frac + 0.41)*dt**, and the
  frac = 0 lag is that law at frac = 0 -- which is why it is O(dt) and dilutes as 1/L.
* So **the clock cannot both null the phase and keep the field a rigid rotation.** That,
  not a sign error, is why "advance the ghosts to t^n + dt" fails.

**REFUTED suspect (do not re-propose).** ConToPrim runs over the ghost zones
(bcs -> prol -> c2p) and REBUILDS bcc from the faces, so for the first ghost cell the
radial average 0.5*(x1f(is-1)+x1f(is)) straddles the boundary: bc_time_frac moves the
ghost cell's angular components by the full omega*frac*dt and its radial one by only half.
Plausible, and wrong. `problem/bc_bcc_match` sets the ghost x1f so that average IS exact,
and it changes nothing (2.16e-4 vs 2.09e-4 at frac 0; 7.21e-4 vs 6.77e-4 at frac 1).
The other old suspect, overwriting the ACTIVE faces x1f(is)/x1f(ie+1), was already fixed.

**STILL OPEN: the -0.41.** `problem/bc_probe` fits the interior's own phase in the first
ACTIVE radial layer (two-parameter fit of the face fields to Rz(theta)B0) and measures
theta = omega*(t^n + dt) in BOTH rk2 stages, to 1-7% -- exactly where SSP-RK2 puts its
stage outputs, and confirming `pmesh->time` advances only at driver.cpp:428. So the ghost
that nulls the lag sits **1.41*dt BEHIND the interior it faces**, not level with it, and
no RK-boundary argument yet explains that number. A worthwhile next check: hold the ghost
phase offset FIXED while halving the CFL -- the law says the -0.906 coefficient is
dt-independent, so the same absolute offset should shift delta by the same amount.

Practically: this is a 0.2% rotation-rate error that shrinks with dt and with a thicker
shell, and it is a property of THIS TEST's boundary, not of the cubed-sphere scheme.

## Audit of the rest of `src/mhd/` for missing cubed-sphere support (user asked)

* `mhd_fofc.cpp` -- **a genuine gap, LATENT.** It recomputes face fluxes AND face EMFs with
  `SingleStateLLF_*` applying none of the three gnomonic transforms, and its trial update
  divides by Cartesian `size.dx1/2/3` instead of area/dxedge. `mhd/fofc` defaults to false
  and no cubed-sphere input sets it, so nothing measured has been affected. NOT fixed.
* `mhd_corner_e_uct.cpp` -- `CornerE_UCT` is called from NOWHERE in `src/`. Unreachable.
* `mhd_etotgrav.cpp` -- coordinate-free (adds phi x mass flux). Correctly needs nothing.
* `mhd_wellbalance.cpp` -- runs BEFORE the gnomonic rotation by construction; the call
  order in `mhd_fluxes.cpp` is deliberate.
* `rsolvers/*` -- by design see an orthonormal frame the caller prepares. Correct as-is.
* `mhd.cpp` -- allocation and parameter reading only.

## Harness -- what `CSTestConvErrors` now prints, and which tool answered what

Runs FIELD-FREE via `problem/conv_errors = 1` (how the iprob=3 control was measured);
`problem/bc_time_frac` sets the ghost time to `t^n + frac*dt` (default 0 = unchanged),
`problem/bc_probe` prints the INTERIOR's phase at the first active radial layer every N
cycles, and `problem/bc_bcc_match` matches the ghost x1f to the bcc ConToPrim rebuilds.
All three default to off and reproduce the plain run BITWISE.
Five diagnostics, each of which decided something:

* **L1(B) budget by region** -- first hint the seam was enriched.
* **profile of L1(B) vs distance from the SEAM, in cells** -- separated "generated here"
  from "advected in from the seam". Combined with `nlim=1` (no wave can travel a cell)
  this is what proved the zeroth-order interior error.
* **profile vs distance from the RADIAL boundary** -- showed the residual is uniform.
* **best-fit uniform B** -- splits systematic (a wrong rotation ANGLE) from structured;
  recovers B0 to 1e-13 at t=0, which validates it. The phase angle comes from this.
* **`omega_fit/omega`** -- the flow's own rotation rate, basis-free. Ruled out "the field
  is fine, the flow is slow".

**Reproducing.** Build `cmake -S . -B <dir> -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release`
(plain bare-login g++, ~1 min at -j48). A run at nx1=16, nx2=nx3=32, tlim=1 takes ~3-6 min
on one core. NOTE: a command-line override can only change a parameter that ALREADY EXISTS
in the input file -- `problem/bc_time_frac` must be added to the athinput before it can be
scanned. This session's runs are in the session-60b484ad scratchpad (`c9/`, `ab/`, `res/`,
`alpha/`, `wscan/`, `thick/`, `flow/`, `bisect/`).

Still open, needs the user: a cs_test regression test cannot join CI until `cs_test` is
promoted to a built-in pgen (the suite builds one binary with no `-D PROBLEM=`). Not
started; ask first.
