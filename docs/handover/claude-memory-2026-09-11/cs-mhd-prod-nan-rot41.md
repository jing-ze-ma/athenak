---
name: cs-mhd-prod-nan-rot41
description: cs_mhd_prod (11516277, 551cd78a binary, the old cs MHD production) went NaN EVERYWHERE between rot 41.5 and 41.6 with NO precursor in the history; kept running on NaN to rot 54+ (the code does not stop on NaN). NaN-hunt restart from the rot-40 rst with 0.05-rot dumps submitted as 11529359 (pending behind the 09-07 maintenance)
metadata:
  type: project
---

**Facts (09-07 09:00).** History: KE/ME/mass flat and normal through rot 41.5 (KE3 1.18e34,
ME3 2.17e33, M/M0-1 -4.9e-3); the rot-41.6 row is NaN in every column; the dump at rot 42
is NaN in all 786432 cells. The run continued with dt 10.4 s (constant) to rot 54.6 --
AthenaK does not abort on NaN, so a job that "runs" can be dead. The event log's int32
counters have OVERFLOWED (eos_efloor = -1.13e9), so they say nothing about the moment.
Restart files every 10 rot: 00004 = rot 40 (last good), 00005 = rot 50 (NaN). Chain jobs
11516227/8 would restart from 00005 -> asked the user to scancel 11516277 11516227 11516228.

**Last good dump (rot 40):** v max 21 km/s at i=66 (upper atmosphere, normal), rho at the
floor 5e-14 in the top layers of most blocks (normal top pinning), |B| max 1.09 kG at the
BOTTOM (i=9, block 6), 402 cells > 300 G all at i < 32. That bottom field is NOT a
runaway: max|B| i<16 = 478 G at rot 10, 964 at 20, 929 at 30, 1090 at 40 (saturated), and
cs_mhd_prod2 has the same 443 G at rot 10. beta there ~ 5e3 (p 225 bar), dynamically weak.
So the rot-41.6 event has no visible precursor in any global or deep quantity.

**NaN hunt:** bench/cs_mhd_prod_nan (copy of rst/dhj.00004.rst, same input and binary),
job 11529359: restart at rot 40, output3/dt = 1.525e4 (0.05 rot), tlim = 1.281e7 (rot 42).
Read the first NaN dump's cell set (block, i, k, j) -> vertex? bottom? pole row of a
panel? floor cell at the top? Then rerun with ncycle_out to the step. Also worth checking:
whether cs_mhd_prod2 (new defaults) passes rot 41.6 -- it is at rot ~11 and continues
after the maintenance.

**HUNT RESULT (09:20):** the restart from rot 40 (job 11529359, 0.05-rot dumps) passed
rot 41.6 with NO NaN (Bmax bursts 1.0 -> 1.5 kG in block 8, i=8, over rot 41.0-41.65,
the only lively deep feature). The rerun's history departs from the original at 5e-6 in
the first 0.1 rot and by 1-4 % in a rotation: a RESTART IS NOT BIT-REPRODUCIBLE for this
cs MHD path (nondeterministic reductions/atomics, or restart state not exact) and the
trajectory is chaotic, so the event was either a rare state-dependent failure this
trajectory did not revisit, or a transient hardware/MPI fault. Cause NOT found.
Follow-ups: (1) determinism test on GPU after the maintenance (two restarts from the same
rst, 200 cycles, bitwise compare) -- if they differ, find the nondeterministic kernel
(cs seam exchange? resistivity? WB cache?); (2) time/nan_check_cycles guard (default 100)
added to Driver::Execute so a NaN run aborts instead of burning the allocation.
