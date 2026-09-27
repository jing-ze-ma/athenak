---
name: overnight-plan-0926-caltech
description: User asleep from 2026-09-26 ~00:00 PDT -- overnight plan on Caltech (merge m1-port if gated, ck experiments report-only, Kokkos 4.6.02 pin bump + re-gate, mg_gf stays held, one morning summary)
metadata:
  type: project
---

The user went to bed on 2026-09-26 at about 00:00 PDT ("can i leave you to it"). The plan I stated and the user accepted:

1. **m1-port:** merge into rt-integration and push, but only if every gate passes: CPU bitwise, GPU vs CPU values, 1 vs 2 GPUs, GPU restart. Otherwise leave it unmerged and write it up.
2. **ck-h200:** experiments only, so do NOT merge. Write up the ranked candidates for >= 25 % for the user to decide.
3. **Kokkos pin** ([[kokkos-pin-4400-not-4602]]): after both agents finish, bump to 4.6.02, do a full rebuild, then re-gate. Push only if clean; otherwise hold and report.
4. **mg_gf** ([[default-flips-0925]]): stays HELD. Report the H200 numbers only.
5. Fetch viper pushes; act only within these tasks. No new side threads.

**Limits:** testing only (short jobs, <= 2 GPUs, timing on H200, correctness may use H100), no production.

**Deliver:** ONE morning summary with results, merges, anything held, and the decisions waiting for the user.
