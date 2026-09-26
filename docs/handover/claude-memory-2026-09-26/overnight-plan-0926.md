---
name: overnight-plan-0926
description: 09-26 overnight (user asleep) - finish m1-wedge + He-box cfl, update handover, merge only gated work, fetch before push, no new threads, morning summary
metadata:
  type: project
---
User went to bed 09-26. Standing permission: continue the running work only.
- **m1-wedge agent** (Caltech TASK-2026-09-26-m1-real-wedge.md; /viper/ptmp2/jinma/sprhd_0926, branch m1-wedge): let it finish. Merge into rt-integration only if its gates pass (keys off bitwise, restart bitwise, existing problems bitwise vs its own base).
- **He-box cfl convergence** (hebox_cfl_0926, R15/H6): put the verdict in the handover.
- Every result goes into docs/handover/HANDOVER-2026-09-26.md (Caltech reads it).
- `git fetch fork` before every push, merge Caltech's commits, NEVER force-push.
- Viper GPU builds must link bench/wt_rgbox/kokkos (4.6.2): the submodule is 4.4.00.
- No new side threads, no production runs. Morning summary only; ignore interim notifications.
Related: [[thread-triage-0925]], [[no-interim-agent-acks]].
- 09-26 night: m1-wedge MERGED (d64b96db pushed). sph_wedge test pgen; gates pass; handover section 4c. Open: ringing, convection onset, job 11984805 unanalysed. He-box cfl verdict still pending.
