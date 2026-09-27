---
name: ck-cadence-retest-0926
description: 09-26 dhj ck retest under the p>1e-6 bar rule - xstep<=4 is a no-op under every 4; xstep 8 -11 % OK; e8x16 -28 % borderline; production itself has 2 % nonconv at the top/T-floor cells
metadata:
  type: project
---
WASP-121b 1x, 3 rotations from sparc_w121x1/base rot 8.14 (/viper/ptmp2/jinma/ckcad_0926/RESUME.md, ana.txt). Pn = P + 1e-14 restart kick (the noise reference).
- ck_impl_xstep counts cycles. Under ck_impl_every 4, xstep <= 4 is bitwise production, so production's xstep 2 does nothing.
- X8: -11 % ms/cycle. Within noise at 1-10 bar; 1.2-1.5x the single-noise spread at 1e-3..1 bar. Acceptable.
- E8X16: -28 %. 1.4-1.9x the noise at 1e-3..1 bar, and dt is set from the top 5 % of the time (P 1.2 %). Borderline: needs a 2nd noise member.
- The old b2 noise (job 11980229) confirms that the earlier 3-rotation deviations of x16 and e8x16 were chaotic noise.
- Production P itself has NOT-CONVERGED in 2 % of calls, mostly the top 3 cells (3e-8..1.6e-6 bar), many at the 200 K T-floor. Possibly the same KKT projection defect as the 10x stall (ck_impl_kkt_row, branch ck-newton10x-b).
- Diagnostics merged 78de21dd: problem/ck_impl_ncloc, seed_restart, dtloc_every (default off). A restart needs -i overlay.athinput to accept the new keys.
- ck_impl_kkt_row (merged 427f9f58, local until the build check): 1x production NOT-CONVERGED 20/2113 -> 0 at the same cost (fresh 9 -> 0). 10x 164 -> 85/200 at +1.6 %; the rest are linearly contracting stiff day-side columns (would need Broyden / banded Jacobian). RECOMMEND kkt_row = true in production inputs.
