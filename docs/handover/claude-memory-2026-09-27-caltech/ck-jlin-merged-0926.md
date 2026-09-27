---
name: ck-jlin-merged-0926
description: ck-jlin (implicit ck speed-up for H200) MERGED 587825ef, pushed 581596f7 (2026-09-26) -- about 21-22 % of ck cost (-15.3 % total on 1 H200, about -13 % on 2); CPU bitwise, GPU round-off, noise 0.89x; the final 2-GPU timing of a861ebb6 is queued (jobs 3524816/7)
metadata:
  type: project
---

User: "merge it now".

**The branch** (steps on f8f7232f):
1. c0c4d2a6: cache the jlin Jacobian coefficients. Reuse calls skip jlin. The biggest win: -9 % total on 2 H200.
2. 2946c3f6: rt_apply maps j-fastest.
3. 269948bf: jlin pass 1 fused into ck_lin_build.
4. b6de2e66: store-kernel launch bound. Reverted in a861ebb6, because it lost on both H100 and H200.
5. 368415b4: direct ck_c0/ci/co reads.

**Result.** The store kernel is now the largest ck kernel (~32 %, 164 registers). It is the next lever: register pressure, splitting its beam/JAC branches. Two smaller items: shrink lP from 9 to 6 slots, and fuse lin_sum into rt_apply.

**Pending**
- The final 2-GPU timing (tf2/tf1 in /resnick/scratch/jingze/ckjlin, jobs 3524816 and 3524817). Read it with `python3 ana2.py tf2 40 400`.
- Viper: HIP rebuild and gate. Noted in the handover.

**FINAL TIMING CONFIRMED (hpc-sm-01-17, 3 interleaved reps, median 8-cycle windows):**

| setup | k46 | a861ebb6 | change |
|---|---|---|---|
| 2 H200 | 32.94 ms/cycle | 28.47 ms/cycle | **-13.6 %** |
| 1 H200 | 59.30 ms/cycle | 49.75 ms/cycle | **-16.1 %** |

- About 21 % of ck.
- Node-to-node base differs: 32.94 on 01-17 vs 35.12 on 02-10, so compare only within one job.
