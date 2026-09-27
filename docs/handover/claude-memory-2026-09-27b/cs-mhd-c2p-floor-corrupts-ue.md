---
name: cs-mhd-c2p-floor-corrupts-ue
description: A REAL cs MHD bug, found, proven and FIXED -- the C2P floor writes back a magnetic energy built from the non-orthogonal triple. Restores 2nd order at the vertex. But it does NOT save the dhj production run (0.2023 -> 0.2129 rot), so it is NOT the cause of that blow-up -- on the cubed sphere ConsToPrim computes the magnetic energy from the NON-ORTHOGONAL face-normal triple; at low beta that drives the internal energy negative, a floor fires and OVERWRITES the conserved energy u.e, and RaiseVelMHD's later correction cannot undo it
metadata:
  type: project
---

**2026-09-03. This is the cause of the cubed-sphere MHD blow-up.** Proven, with the fix
direction demonstrated. See [[cs-mhd-low-beta-divergent]] for the reproducer and the
thirteen mechanisms eliminated on the way.

## The chain

1. `ConsToPrim` (`SingleC2P_IdealMHD`, `src/eos/ideal_c2p_mhd.hpp`) computes
   `e_m = 0.5*(bx^2 + by^2 + bz^2)` from `u.b`, which on the cubed sphere is the
   **NON-ORTHOGONAL face-normal triple** (B.rhat, B.nhat_xi, B.nhat_eta).
   `nhat_xi . nhat_eta = -c`, so the sum of squares is NOT |B|^2. The error is O(1) and
   largest where |cos_cell| is largest -- **the cube vertex**.
2. At low beta that wrong (too large) magnetic energy drives `w.e = u.e - e_k - e_m`
   NEGATIVE, so a floor fires.
3. **The floor OVERWRITES the conserved energy**: `u.e = efloor + e_k + e_m`, using the
   WRONG e_m. The temperature floor does the same; disabling one just makes another fire.
4. `Coordinates::GnomonicEquiangleRaiseVelMHD` then rebuilds bcc ORTHONORMALLY and
   subtracts the CORRECT magnetic energy -- from the already-corrupted `u.e`. The residual
   is exactly (e_m_wrong - e_m_correct): **O(1), non-convergent, at the vertex.**

The code documents the hazard without recognising its severity, in RaiseVelMHD's own
header: *"the floors are left exactly as ConsToPrim applied them -- it tested them against
an internal energy that lacked both the metric cross term and the corrected magnetic
energy."*

## The proof

`cs_test` iprob=8 seeds u0(IEN) with the ANALYTIC 0.5*b0c^2, so any error in the
subtracted magnetic energy shows directly in the recovered pressure. Since bcc is LINEAR
in b0c, dp/p MUST scale as b0c^2 -- any departure is a branch. Measured at nx=32:

    b0c   beta     dp/p        dp/p / b0c^2      eos_efloor
    1     2        2.008e-04     2.008e-04            0
    2     0.5      8.030e-04     2.008e-04            0        <- clean b0^2
    4     0.125    1.339e+00     8.370e-02         7 700       <- floor starts firing
    16    0.008    3.643e+01     1.423e-01        29 948

Disabling only the u.e OVERWRITE (clamp w.e, leave u.e alone) restores it exactly:

    b0c   dp/p default   dp/p keep u.e   / b0c^2
    2      8.030e-04      8.030e-04     2.0076e-04
    4      1.339e+00      3.212e-03     2.0076e-04
    16     3.643e+01      5.139e-02     2.0076e-04     <- 709x better, IDENTICAL constant

And the real gate, per-step vertex convergence order (nlim=1):

    beta     default      keep u.e
    0.1       -0.04        **+1.97**
    0.01      -0.02        **+2.00**

**Second order restored.** The scheme is consistent again.

## The fix, COMMITTED as a46b75e0

**`GnomonicEquiangleRaiseVelMHD` ALREADY re-applies the floors correctly** to the corrected
internal energy and writes `u0(IEN)` back, for both the ideal and the general EOS -- the
machinery was there all along. The only problem was that ConsToPrim had already corrupted
`u.e` before it ran. So the fix is small:

  * `EOS_Data::defer_cons_floors` (eos.hpp), set in the EOS ctor to
    `pp->pmesh->use_cubed_sphere`;
  * the four Newtonian-MHD conserved-energy writes inside the floor branches --
    two in `ideal_c2p_mhd.hpp`, two in `general_c2p_mhd.hpp` -- guarded by
    `if (!eos.defer_cons_floors)`.

The PRIMITIVE is still clamped exactly as before; only the CONSERVED write is deferred.
Off for Cartesian and spherical polar, so those answers are bit-identical. The HYDRO half is guarded too (the user asked for it), for structural uniformity --
**not** for a demonstrated accuracy gain. MEASURED with the hydro gate below: with 15904
floor firings in a single step the cube-vertex order is +1.08 against the interior's
+1.31, i.e. NO consistency loss, and the guarded result is BIT-IDENTICAL. That identity
is expected: when a cell ends up floored either way, deferring the write lands on the same
value; MHD differed only because its two magnetic energies genuinely disagree. So
`cs_prod_hyd` answers should be unchanged -- worth confirming on a restart before relying
on it.

**A hydro gate now exists**: `problem/conv_errors = 1` makes `CSTestConvErrors` run for
iprob=3 (the exactly-steady rigid rotation), giving region-split L1 for hydro. Lower `p0`
to raise KE/e_int -- the hydro counterpart of lowering beta. That is what turned "I think
hydro matters too" into a measurement; my first recommendation (fix it for accuracy) was
WRONG and was withdrawn on the strength of it.

RESULTS after the fix:

    per-step VERTEX order   beta=2 +2.03 | 0.5 +1.96 | 0.1 +1.97 | 0.01 +2.00
    many-step vertex order  beta=0.1: -2.0  ->  +0.98   (divergent -> convergent)
    dp/p / b0c^2            2.0076e-04 IDENTICALLY at b0c = 2, 4, 16

**Gate any fix with:** `cs_test` iprob=9, `nlim=1`, per-step VERTEX order at beta = 0.1
and 0.01 must be ~+2, not ~0. The defect is invisible at the default beta = 2, and a
many-step run muddies the order with accumulation.

## Why it took so long, and the lesson

Every earlier suspect was in the flux/EMF/halo path, and **all of them were innocent** --
the error is present at `nlim=0`, before a single flux or EMF is computed. That one fact
would have redirected the search immediately. **When a defect shows up in a state that has
only been INITIALISED, the transport operators cannot be responsible; look at the state
representation and the conversions.** Cf. [[validate-the-instrument]].


## CRITICAL CAVEAT: this does NOT fix the dhj blow-up

Validated on the production case (bench `cs_floorfix`, cs + MHD + bbot = 3, the arm that
died at 0.2023 rot). **With the fix it still dies, at 0.2129 rot -- 5% later.**

So: the defect is REAL, the unit gates prove it is FIXED (per-step vertex order -0.04 ->
+1.97, many-step -2.0 -> +0.98, dp/p back to a clean b0c^2), and the fix is worth keeping
on correctness grounds. **But it is not the cause of the cubed-sphere hot-Jupiter
blow-up**, which survives it almost unchanged. Do not describe this commit as fixing that.

That is the [[measure-impact-before-claiming]] lesson again, and it nearly bit hard: every
unit measurement said "found it", and only the production A/B showed otherwise. **A
mechanism proven in a unit test is not thereby the cause of a production failure -- run
the failing case before claiming it.**

The dhj trajectory is essentially UNCHANGED up to the failure (dt 6.52 vs 6.56 at 0.130
rot), so whatever kills that run is something else, and it is not gated by beta at the
vertex in the way iprob=9 is. [[cs-mhd-dhj-blowup]] stays OPEN.
