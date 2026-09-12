---
name: restart-bc-reads-w0-bug
description: 2026-09-12 RESTARTS WERE NOT BITWISE because the user physical BCs read PRIMITIVES (w0 / wtemp) -- Driver::InitBoundaryValuesAndPrimitives calls ApplyPhysicalBCs BEFORE the first ConToPrim after a restart, when w0 is still zero (the pgen returns early on restart), while a running stage sees the previous stage's w0; u0 (incl. ghosts) is restored bitwise. red_giant fixed (3fd3a836: gates + state_i read u0, bc_use_cons removed); deep_hot_jupiter_rt has the SAME defect and worse (divides by w0 density, TGuess from wtemp) -> every chained dhj production run restarts non-bitwise with inf/NaN-prone outer ghosts for one stage
metadata:
  type: project
---

Found with the P/E/F0/F one-cycle test on the GPU (bench/RG_fofc_long/bitwise: continuous vs restarted, dumps
with ghost_zones=true, per-class comparison): active cells identical at the restart cycle, ALL outer radial ghosts
differ (rel up to inf), next cycle the top active rows differ. Every physics switch left it (bisect round 1+2).
Rule: a user BC must read u0 (+ b0), never w0/wtemp/bcc; convert locally (utils/eint_from_cons.hpp subtracts KE,
rho*Phi, magnetic energy). Removing the read also removes a one-stage lag in the running boundary (small change
to running results; measure it). Test recipe: bench/RG_fofc_long/bitwise/fixtest/submit.sh (red giant),
bench/cs_mhd_prod3/bitwise/ (dhj, being staged). See [[rg-fofc-gpu-survived]] for the trail.

**Status 2026-09-13 ~01:30.** red_giant: 3fd3a836 (BC gates/state read u0, bc_use_cons removed) + ec5e1f6d (driver: on
restart KEEP the file's ghosts, run only ConToPrim; <time>/restart_refill_ghosts=true = old; AMR path unchanged) --
in worktree, being validated against the suite's restart tests before merge. dhj: af2f87d2 -> e4b318b3 PUSHED
(outer fill: rho/e_int/T from u0 via EintFromCons + eos inversion, no wtemp warm start, no active-column c2p inside
the BC that edited u0(IEN,ie) on every call; inner: e_mag from ghost b0 faces before the reflection; trajectory impact
<= 2.7e-6 over 11 cycles; ck test passes). The dhj outer fill did NOT divide by a zero w0 (it ran c2p on the ie column
first) -- the earlier claim was wrong. RESIDUAL in both problems: ~1e-7, = the general-EOS c2p WARM START from
wtemp (not in the restart file; the root find stops at a tolerance). Fix = converge the inversion to machine
precision or write wtemp into the restart. ideal-EOS dhj restarts are bitwise with the fix.
**MERGED + PUSHED 2026-09-13 ~02:30:** 3fd3a836 (red_giant BC reads u0), ec5e1f6d (driver: restart keeps the file's
ghosts, runs only ConToPrim; <time>/restart_refill_ghosts=true = old; AMR unchanged), 8a665c40 (general-EOS wtemp
cache written to / read from the restart; ideal-EOS files unchanged; legacy files load with a warning + one cold
start). Residual after all three (fixtest round 3, 11654134): ONE outer radial ghost eint cell at 1 float32 ulp at
the restart cycle (BC writes w0 ghosts directly; the restart recomputes w0 = c2p(u0) -- EOS round trip), which the
implicit radial conduction spreads to ~1e-6 in one cycle. Left OPEN. CPU tabulated-EOS linear wave: bitwise.
The suite has NO restart tests -- add one.
