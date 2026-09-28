---
name: cs-gs07-emf-gate
description: cs_gs07_emf MEASURED on the cs_test iprob=9 convergence gate (2026-09-05) -- same order as the plain average but 1.8x LARGER L1(B) on a smooth field at every resolution and beta; its value is robustness (low-beta arm matrix), not accuracy. The old "worth 30x" is NOT what this gate shows
metadata:
  type: project
---

Recipe: build_cs (cs_test), inputs/tests/cubed_sphere_mhd_conv.athinput + `cs_gs07_emf` declared
in the file (command-line overrides need the parameter present), iprob=9, p0=1, raw HLLD
(cs_lowbeta_fallback=0), nx2=nx3 = 16/32/64, b0c = 1.0 (beta 2) and 4.472 (beta 0.1); per-step
(nlim=1) and finite-time (tlim=0.02236). Sweep script: scratchpad gs07_gate/one.sh.

    L1(B) total at n64        per-step (order 16->32/32->64)     finite-time
    beta 2   plain average    1.95e-7  (2.71/2.63)               1.72e-6 (1.60/1.38)
    beta 2   GS07             3.44e-7  (2.65/2.62)               2.98e-6 (1.50/1.41)
    beta 0.1 plain average    1.89e-7  (2.79/2.74)               5.88e-6 (1.44/1.25)
    beta 0.1 GS07             3.80e-7  (2.79/2.73)               7.06e-6 (1.75/1.45)
    cube-vertex region: same picture, GS07 1.8x larger, orders ~2.0 per step for both.

So GS07 is MORE DIFFUSIVE on a smooth field (an upwind EMF, as expected) with identical
order. The finite-time order below 2 is the same for both -- the radial-BC phase lag of
[[cubed-sphere-mhd-convergence]], not the EMF. Its justification is the low-beta arm matrix in
[[cs-dhj-production-retry]] (+50 % time to death alone; wb+emf never died), i.e. ROBUSTNESS.
The "worth 30x on its own" line in older notes is not reproduced by this gate and should be
read as referring to the low-beta per-step residual at the vertex, not to smooth-field accuracy.

**Why it stays a switch (2026-09-05):** the GS05/GS07 upwind terms are exactly the operator
that grew the sp polar face-field mode ([[sp-polar-field-blowup]]), and the cs deep zonal sheet
([[cs-deep-toroidal-sheet]]) is a thin face-field structure in a run that has GS07 ON. Do not
flip the default until bench/cs_deep_ablate/nogs07 has read.

**SWITCH REMOVED (b71ef392, 2026-09-05, user decision): GS07 is unconditional on the cubed
sphere.** `<mhd>/cs_gs07_emf` no longer exists (an input line with it is silently ignored);
the plain average is reachable only as the `<mhd>/bs_emf` diagnostic. Verified: no flag ->
L1(B) 2.112e-06 (GS07), bs_emf=true -> 1.207e-06 (average) on the nx32 beta-2 per-step gate.
Binaries BEFORE b71ef392 (e.g. the production 97bf06ab) still have the switch, default off.
