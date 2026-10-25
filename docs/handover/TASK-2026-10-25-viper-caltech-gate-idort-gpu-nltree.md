# TASK viper -> caltech (USER): gate the combined branch idort-gpu-nltree NOW, in parallel with Raven (speed)

The user wants the latest Arepo2-IDORT gated as fast as possible. Raven's gpudev runs one of our jobs at a time (~2-3 h of queue).
Please run the same gates on any free GPU slice at Caltech. Whichever machine finishes first counts.
No Arepo source in git, as before.

**Code: NOT yet gated.** The user copies the bundle to Caltech:
- Raven path: `/raven/ptmp/jinma/arepo2_gpu/caltech_bundle_1025/arepo_idort_gpu_nltree_1025_caltech.bundle`, md5 ff1e715df343496cfcfded67c12ab52b, 4.5 kB.
- Prerequisites, all in your clone: nltree-g25 1108335a71, idort-gpu-s1 6cef4d1d69, idort-gpu-caltech-opt c663a54c5.
- `git fetch <bundle> idort-gpu-nltree:idort-gpu-nltree` gives e93504ccef.

**What the branch contains:**
1. 08491b0559: merge of nltree-g25 (621fbe3665) into idort-gpu-s1.
2. 3405e3ab39: the g25 fix. `gravity_direct` skips NO_SELFGRAVITY_TYPE sources and zeroes NO_GRAVITY_TYPE sinks, as the tree walk does.
   Under ALLOW_DIRECT_SUMMATION, gas used to get the core's gravity twice.
3. e93504ccef: merge of your idort-gpu-caltech-opt.

**Gates, as in your results 1 and 2:**
1. The 10 IDORT gates, IdortGpu 0 vs 1: all ALL SAME.
2. nltree_multibin_3d(_cuda) and nltree_mns_3d(_cuda): PASS.
3. g25_red, 2 top steps: c1 (CPU node-local tree) vs g2 (full GPU recipe). Compare idort_solver.txt and the final state; report whether they are bitwise.
   Note: the direct-summation fix changes g25 results against the pre-fix code, by design.

**Reply** with a NOTE as soon as gates 1 and 2 are done; gate 3 can follow.
- Do not wait for a full node.
- Your slices are fine (32 ranks).

## WITHDRAWN (USER 10-25): do not run this TASK
The user does not want Caltech to run these gates. Raven will gate the branch.
- If you have started, stop and discard.
- No reply is needed beyond a one-line ack.
