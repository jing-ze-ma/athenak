---
name: test-state-use-prod4-restarts
description: User rule 09-23 - dhj tests (implicit ck, RT, cost/timing arms) start from prod4 restarts (bench/cs_mhd_prod4/rst, read-only), not prod3's (odd-even columns); hydro tests from cs_hyd4_prod
metadata:
  type: feedback
---

User 2026-09-23: "the prod3 restart is not very good. maybe try prod4 restart from now on".
**Why:** every prod3 state after rot 0 has mildly flagged odd-even columns (180-340), and the rot-283
restart dhj.00567.rst has 9 strong ones (a night-side corner column, MeshBlock 13) that the implicit ck
solve cannot converge and where T3/T4 differ by up to 56 % (tests_ck_implicit/README_T1_stall.md,
README_T3.md, README_T4.md). prod3 also predates the prod4 physics changes (tm sweep, spherical beam,
polytropic WB, pfloor 1e-5).
**How to apply:** use the newest prod4 restart in bench/cs_mhd_prod4/rst/ (read-only; copy the file into
the test dir, never write into the production dir) with the prod4 input; state its rotation. Early prod4
states (rot < ~20) are still spinning up (dt ~4 s, ME growing), so for "developed state" questions say so
and prefer a later restart when one exists. Old timing numbers from the prod3 rot-283 restart are not
directly comparable with new prod4-restart numbers: rerun the reference arm. Hydro tests: bench/cs_hyd4_prod/rst.
Also check a new restart for odd-even columns before trusting it (scan in bench/impl_t3_0923/scan).
See [[cs-mhd-prod4-run]], [[cs-hyd4-prod-run]].
