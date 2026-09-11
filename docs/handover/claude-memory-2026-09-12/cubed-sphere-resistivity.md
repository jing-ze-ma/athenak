---
name: cubed-sphere-resistivity
description: Resistivity now works on the cubed sphere (4bfacdd8) -- the current density was falling through to the CARTESIAN curl; new resistivity_gnomonic.cpp takes a two-pass curl on the non-orthogonal basis, plus the Poynting flux and the diffusive dt, validated by the new cs_test iprob=11
metadata:
  type: project
---

Done 2026-08-29, committed as **4bfacdd8**. Asked for after the support audit in
[[cubed-sphere-mhd-convergence]] found resistivity half-implemented.

## What was broken

`current_density.hpp` branches on `use_spherical_polar` ONLY, so the cubed sphere fell
through to the **Cartesian** curl dividing by `size.dx1/dx2/dx3` -- index spacings, not
lengths (x2/x3 are equiangular coords on [-1,1]). O(1) wrong, and unguarded, even though
`resistivity_ct.cpp` and `resistivity_update.cpp` already had `use_cubed_sphere` branches.
Its curvilinear branch could not be reused either: **the tangent basis is not orthogonal**,
`b0.x2f` holds B.nhat_xi, and Stokes needs B.e_xi = (B.nhat_xi + c B.nhat_eta)/s -- an
O(c) correction, O(1) at a panel corner. Same trap as `mhd_corner_e.cpp` (e63b571a).

## The fix

* **`src/diffusion/resistivity_gnomonic.cpp`** (new). Two passes, which cannot be fused
  because the rotation mixes J at x2 edges with J at x3 edges (different points):
  pass 1 = Stokes on the loop through the four cell CENTRES around each edge -> J in the
  FACE-NORMAL frame; pass 2 = rotate into the edge frame (E.that, what `mhd_ct.cpp`
  consumes and what `GnomonicEquiangleEmfX1` produces) and multiply by eta. Entered from
  a branch at the top of `AddEMFGeneralResist`. Needs `jnorm`, an extra edge field
  allocated only for the cubed sphere.
* **`CoordGnomonicEquiangle`** now fills the DUAL mesh -- `Coordinates::dxface` and
  `::areaedge`, which existed but were spherical-polar only. The x1 loop is in the tangent
  surface so its area carries sin(basis angle); the x2/x3 loops contain rhat and do not.
* **Poynting flux** crossed covariant edge EMFs with orthonormal `bcc`. Only x1 and x2
  needed work: nhat_eta IS the third orthonormal axis, so the x3 flux was already right.
* **Diffusive dt** used `mb_size.dx*` (index spacings); now `Coordinates::dx{1,2,3}`.

## The test: `cs_test` iprob = 11, `inputs/tests/cubed_sphere_resist.athinput`

B = b0c*(-y,x,0) + a uniform tilted field. curl B = **2*b0c*zhat exactly**, so eta*J is
known at every edge INCLUDING its projection on that edge's own unit tangent -- which is
the thing a metric-free current density gets wrong by O(1). The uniform part adds nothing
to the curl but makes every Stokes term nonzero, so a dropped term cannot hide. The state
is static and heats uniformly at eta*|J|^2 -- the only quantitative check on the Poynting
flux. Measured, nx1=8, nx2=nx3=16/32/64, EMF L1 relative to |eta*J|, zero-dt limit:

    x1e  1.450e-3 / 3.733e-4 / 9.491e-5   ratios 3.88, 3.93
    x2e  1.186e-4 / 2.906e-5 / 7.229e-6   ratios 4.08, 4.02
    x3e  1.904e-4 / 4.614e-5 / 1.161e-5   ratios 4.13, 3.97

x1e splits into panel INTERIOR (8.38e-4 / 2.14e-4 / 5.41e-5, ratios 3.91/3.96 = clean 2nd
order) and the 2-cell panel-edge RING (2.31e-3 / 9.13e-4 / 3.96e-4, ratios 2.53/2.31 =
FIRST order), worst point always a cube VERTEX. That is inherent: differencing the seam
halo's own O(h^2) interpolation costs one order. Domain L1 stays 2nd order because the
ring shrinks. Ohmic heating at t=0.2 vs (gam-1)*eta*J^2*t: **0.99893 (nx2=16), 0.99996
(nx2=32)**; with b0c -> 0 the spurious heating is 4.6e-8, 2e-5 of that scale.

## Three traps the TEST had -- all cost real time, all now in the file

1. **Measure the EMF at a tiny CFL** (nlim=1, cfl~0.003). Otherwise the one-step evolution
   of the field swamps the EMF truncation error, amplified by 1/dr -- it looked like a
   first-order error that got WORSE when the radial grid was refined.
2. **The initial energy must use the ANALYTIC |B|^2/2**, not the face-averaged `bcc`
   (which is unset at that point in the pgen anyway). Otherwise the state is not in force
   balance and the shell rings, with a non-converging drift 2.7x the heating signal.
3. **The ghosts must heat with the interior** (`dp_ohm` in `CSTestRadialBC`). A static
   Dirichlet pressure lets the heated shell push against a cold boundary; the leak grows
   like t^2 and dilutes as 1/L, and read exactly like a 10% heating DEFICIT. Diagnosed by
   the [[localise-by-dilution]] rule -- vary the shell thickness at fixed dx.

## Still NOT supported -- but now FATAL at startup (9f79a031)

`hydro_fofc.cpp` / `mhd_fofc.cpp`, `viscosity.cpp`, `conduction.cpp` and SMR/AMR
(`mesh_refinement.cpp` / `prolongation.cpp` have no cubed-sphere path). Each used to run
to completion and return a wrong answer; the guard in `Mesh`'s ctor now names what is
missing and exits. All four verified to fire, a supported cs run untouched, AMR off the
cubed sphere unaffected. Resistivity is the model for the work any of them needs.

**TRAP when writing such a guard:** `pin->GetOrAdd*` CREATES the block it is asked about,
and `AddPhysics` builds a module for every block that exists -- probing `hydro/fofc` that
way made a pure-MHD input die on a missing `hydro/eos`. Use `DoesParameterExist` first.
Related: a command-line override can only change a parameter that ALREADY EXISTS in the
input file, so testing these guards needs edited input files, not `-` overrides.

**RKG super-time-stepping (`mhd/use_rkg_sts`) is VALIDATED too** (7e290c59, no code
change). Same iprob=11 problem, nx2=32, tlim=0.2, eta=0.5: heating ratio 0.996605 vs the
direct path's 1.000567, EMF errors agreeing to <1%, for **61x fewer MeshBlock-cycles and
10x less CPU** (138 vs 8406 cycles, 24.6 s vs 247 s). The 0.34% is the super-stepping's
OWN temporal error, not geometry: at fixed grid it falls 0.99661 / 0.99862 / 0.99965 as
the CFL is halved twice = second order in dt. Both null tests match the direct path to
the digit.

**The super-stepping is FIRST ORDER IN TIME, and only in the ENERGY** (f6792ac4; the
"second order" in 7e290c59 is RETRACTED). Scanning dt at FIXED grid over 16x, cfl
0.3/0.15/0.075/0.0375/0.01875, the heating ratio is
0.996605 / 0.998619 / 0.999654 / 1.000177 / 1.000440; against the self-consistent dt -> 0
asymptote 1.000706 the errors halve each time, implied orders 0.975/0.988/0.992/0.992.
Over the same scan the EMF errors and max|v| are dt-INDEPENDENT and converge to the direct
path's values, so induction and momentum are fine -- it is the Ohmic HEATING that is
applied to first order. Suspect the split between the ideal flux (frozen, weight f0) and
the resistive flux (weight fjm1, zero at stage 1) in `resistivity_update.cpp`; NOT
verified, and the user called a halt before it was chased -- **do not reopen unasked**.

Two things were established before stopping, both worth keeping:
* **The dt error is a pure function of dt.** Isolating it as D = ratio(cfl 0.075) -
  ratio(cfl 0.3) (the spatial error cancels in the difference), D = 3.020e-3 / 3.033e-3 /
  3.049e-3 / 3.083e-3 for eta = 0.125 / 0.25 / 0.5 / 1.0 -- constant to 2% while eta
  varies 8x and the RKG substage count s (~ sqrt(eta)) varies ~3x. So it does NOT track
  the sub-integrator's stage count.
* **Strang splitting would probably not fix it, and this test could not tell.** Lie and
  Strang differ only through the commutator [A,B]. Here the state is an exact equilibrium
  of the ideal operator A, and what the resistive operator adds is a SPATIALLY UNIFORM
  pressure, which A also annihilates (only grad p enters) -- so [A,B]u is ~0 in this
  problem and no splitting arrangement should change the measured number. By the same
  token this test is BLIND to splitting error, including in B (the exact solution is
  resistively static, curl(eta J) = 0), so the "energy only" result above must not be
  read as proof that the induction has no splitting error.

**METHOD TRAP, this is what fooled me:** do NOT measure the temporal order against the
ANALYTIC value. The grid's own dt-independent spatial error (+5.7e-4 at nx2=32) then sits
in the residual and makes a first-order sequence look second order. Solve for the dt -> 0
asymptote from the sequence itself, or difference against the dt-converged direct path.

`mhd/use_rkg_sts = true` is enabled in `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`
(spherical polar), so this matters for production: confirm it there before relying on
Ohmic heating in an STS run.

## RESISTIVITY ACROSS A RADIAL BLOCK BOUNDARY -- the alarm was MY GATE. RETRACTED.

Asked whether SMR + resistivity works. First answer was "no, x2e/x3e are 2x off in L1 and
44x in max as soon as x1 is split". **That is RETRACTED.** It was a measurement error, the
third of this shape in this project.

**`CSTestResistCheck` reads `pres->efld_resist`, which is PRE-EXCHANGE.** The gnomonic
curl fills it inside `AddFluxGeneralResist`; `AddEMFDirect` then merges it into `efld`,
and only after that do `SendE`/`RecvE` AVERAGE the two sides of a shared block face. A
block-face value in that array is *supposed* to be two-valued -- reconciling it is what
the averaging is for -- so its block-face maximum is not evidence of anything.

**The gate that means something, now added to `CSTestResistCheck` (EVOLVED FIELD error).**
For iprob=11 eta*J is a CONSTANT vector, so curl(eta*J) = 0 and the resistive term cannot
move B at all; B drifts only through -v x B with v ~ 1e-4. So B(t) = B(0) to that
accuracy, and |B - B_exact| over all three face components is a PHYSICAL norm with no
reference to any block -- it must be decomposition-independent. Measured, iprob=11,
nx2=16, refinement = none, changing only `<meshblock>/nx1`:

| | L1 | Linf |
|---|---|---|
| 1 radial block | 4.0658e-05 | 8.074e-04 |
| 2 radial blocks | 4.0400e-05 | 9.903e-04 |

**Agreement.** The averaging does reconcile the pre-exchange discrepancy, and resistivity
across a same-level radial block boundary is fine.

**THE LESSON, and it is the SAME one as the rank-local norms:** before believing a gate,
find out WHICH ARRAY it reads and at what point in the task list. `efld_resist` is
pre-exchange; `u0`/`b0` at the end of a run are post-everything. Prefer a norm against an
exact solution over any per-block intermediate.

**Structural audit (asked for explicitly).** The resistivity files DO mirror the MHD ones:
`AddEMFGeneralResist` dispatches to `AddEMFGnomonicResist`; the Poynting flux in
`AddFluxGeneralResist` has a cubed-sphere branch that rotates the edge EMF into the
orthonormal frame before the cross product; `resistivity_ct.cpp` and
`resistivity_update.cpp` carry area/dxedge forms. The two routines with no gnomonic form
(`AddEMFConstantResist`, `AddFluxConstantGridResist`) are DEAD CODE with explicit FATAL
tripwires. The resistive EMF also rides the SAME `efld` buffers as the ideal one, so the
x1x2/x3x1 edge-buffer fix ([[cubed-sphere-smr]], 74cbc8df) covers it too -- which is why
that fix moved the pre-exchange x1e back to its one-block value exactly.

**RULED OUT: stencil width.** nghost=3 gives x2e L1 3.4516e-04 against nghost=2's
3.4515e-04. Do not re-try.


## AND THE RETRACTION NEEDS ITS OWN CORRECTION: there IS a small residual bug. **CLOSED 397b4ad3 -- see [[cs-cube-vertex-corner-radial-ghost]]. It was the CUBE-VERTEX corner halo, whose extrapolation ran only over the ACTIVE radial range, so a radial split left its radial-ghost layers written by nothing. Not an EMF defect at all; the resistive curl was merely the only stencil that reads that block. Linf and the heating ratio are now decomposition-independent at 1/2/4 radial blocks. Everything below is the record of the hunt.**

Having retracted the O(1) alarm above, I then called 1-block vs 2-block "essentially
identical" and moved on. **The user caught that, and was right.** The bar is the one that
certified ideal MHD: a radial split reproduces L1(B) to ALL SEVEN PRINTED DIGITS
(5.277894e-04 both). Resistivity does not:

| EVOLVED FIELD error, iprob=11 nx2=16 | L1 | Linf |
|---|---|---|
| 1 radial block | 4.0658e-05 | 8.074e-04 |
| 2 radial blocks | 4.0400e-05 | 9.903e-04 |
| refined (level boundary at r=1.5) | **2.9483e-05** | 2.482e-03 |

**Linf moves 23% on a plain radial split. That is not round-off, and a correct scheme
cannot care how the mesh is cut into blocks.** Note the L1 comparison is CONFOUNDED --
the face counts differ (39936 vs 41472) because a shared block face is counted twice in
the 2-block mesh, so the two L1 norms are not over the same set. **Linf is the clean
number**: a max is immune to double counting. Fix the L1 double count if a sharper number
is wanted.

**Character of what is left**: small, not structural. Linf 2.5e-3 refined is 0.25% of
|B|, nowhere near the O(1) the pre-exchange array implied, and L1 IMPROVES 1.38x with
refinement, which is what extra inner-shell resolution should buy. So this reads as a
consistency error at the radial interface, not a missing transform.

**LOCALISED, first run (790d02c0).** The gate reports where the max sits and tags i == 0.
iprob=11, nx2=16, tlim=0.02, changing ONLY `<meshblock>/nx1`:

    1 radial block    Linf 5.885e-04   on x2f at (i,j,k) = (5,0,15)          mid-block in r
    2 radial blocks   Linf 7.441e-04   on x2f at (i,j,k) = (0,0,15), x1=1.5..2
                                       **ON the block's INNER RADIAL face**

**So the error the split ADDS sits exactly on the radial block interface, on the
TANGENTIAL component x2f, with j = 0 and k = 15 -- both on tangential block edges.** That
is the same RING and the same component family as the edge-buffer defect fixed in
74cbc8df: the radial interface where it meets the tangential boundary. **START HERE next
session**: the x1-FACE flux buffer (n = 0/4) carries x2e and x3e on that face and applies
its seam transform only when `lev == mblev` and the panels differ; the x1x2/x3x1 EDGE
buffers carry them on the ring itself. Read what the resistive contribution does
differently from the ideal one at exactly (i = is, j = js, k = ke) -- and remember the
resistive EMF is ADDED to efld before SendE, so it rides the same buffers.

**Do NOT re-run these, they are done:** nghost width (ruled out), the structural audit of
the resistivity files against the MHD ones (they DO mirror them -- see above), and the
pre-exchange `efld_resist` numbers (they measure the wrong thing).

**CLOSED 397b4ad3** (2026-08-30). The residual radial-block-interface bug above was the
cube-vertex corner halo, not the flux buffers -- see [[cs-cube-vertex-corner-radial-ghost]]
for the account. **c4735f18** then made resistivity + REFINEMENT a supported combination
(it had been a startup fatal). The "START HERE next session" note above is spent.
