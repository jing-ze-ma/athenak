---
name: orion-merge-roundoff-8da093f5
description: 2026-09-11 orion merge a10e367d changes dhj results at 1-ULP level; bisected to 8da093f5 (EOS floors commit), codegen in GnomonicEquiangleRaiseVel, all new switches off
metadata:
  type: project
---

After pulling the orion fixes (c26b01ac -> a10e367d, HIP build fix cdd7d2a5 on top), the dhj RT run is
NOT bitwise reproduced on either cs or sp: cs first differs at cycle 12 by one float32 ULP in vely in one
cell, sp by 1e-12 at cycle 100; both grow chaotically to 1e-3..1e-2 by cycle 400; conserved quantities
stay within 1e-10. Old binary vs itself is bitwise identical (GPU determinism holds).

Bisect (cs 20-cycle gate, per-cycle dumps, bench/ab_orion/cs/bis_*): culprit is 8da093f5 "EOS floors:
pay floors from kinetic energy, density-floor momentum scaling, velocity ceiling". Its new behaviours
(efloor_from_ekin, dfloor_keep_velocity, vceil, eos_floor_consistent) are switch-gated and OFF in the dhj
inputs; the 1-ULP entry is codegen: GnomonicEquiangleRaiseVel (src/coordinates/coordinates.cpp) was
rewritten from par_for to a parallel_reduce with a hoisted `ekin`, letting hipcc contract/reassociate
the m.v products differently. Commits 1-4, 6-9 and the merge resolution are innocent.

**Why:** any A/B against pre-a10e367d binaries (sp_mhd_prod3 runs 9a9396f7) is now round-off-seeded;
a restart of sp_mhd_prod3 on the new binary will diverge chaotically but not physically.

**How to apply:** treat the new binary as equivalent; do not chase 1e-3 differences at cycle 400 in
cs/sp comparisons across this merge. NOT independently verified that the new counters never fire in the
run (an event-log check would close that). Worktrees bench/wt_{k3,k6,nEOS,nRT} can be removed.
Related: [[cs-ulp-amplification]], [[dhj-run-to-run-nondeterminism]].

## GATED, caad9247 (2026-09-11)
`EOS_Data::floors_legacy` (true iff dfloor_keep_velocity, vceil, eos_floor_consistent, efloor_from_ekin all off)
selects the pre-8da093f5 kernels verbatim: GnomonicEquiangleRaiseVel, general/ideal ConsToPrim (new switch-aware
launches moved to src/eos/{general,ideal}_hyd_floors.cpp), SingleC2P_*Legacy, prolong_prims. The hunk that
actually mattered was the GENERAL-EOS c2p, not cs_raisev. Trap: a legacy kernel in the SAME translation unit as
the new one still differed (hipcc inlines both); the kernels must live in separate .cpp files. Gates: cs and sp
20-cycle ndiff 0 vs c26b01ac. Default path is now bitwise the old one.
