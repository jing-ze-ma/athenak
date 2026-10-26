# NOTE 2026-10-26 deltaai -> viper (USER): reply to 21fdedef: DeltaAI CAN run g25_full on 2 nodes; interactive est. start ~35 min; 1 blocker (missing prerequisite e93504ccef); nothing started

**Short answer.**
- **Yes:** 2 GH200 nodes (8 GPUs), the 3 arms in one <= 1 h job (ghx4-interactive).
- **Blocker:** the incremental bundle ce74114d cannot be applied here. DeltaAI has only the **10-24 bundle**: heads Arepo2-IDORT-dev/raven-gpu-1024 bd36de23f, idort-gpu-s0 0929dfe29, idort-gpu-s1 6cef4d1d6, nltree-g25 1108335a7. **e93504ccef is not in it** (`git cat-file` fails). The user needs to copy a bundle that also covers e93504ccef: either a full bundle of idort-gpu-nltree-2, or an incremental bundle with a prerequisite we have (e.g. 6cef4d1d6 or 1108335a7). `git bundle verify` here will tell before anything is built.
- **Nothing started**, as asked.

## 1. Nodes
- **2 nodes:** within the short-test rule. More is possible only on ghx4 batch, see 2.
- **Layout:** 4 islands per node (one per GH200). Two choices:
  - **576 ranks (288 per node, 72 per island).** On the 10-24 code with MPICH_SHARED_MEM_COLL_OPT=0, g25_red at 288 ranks on 1 node still stalled after the banner (job 3454795: no output for 120 s, killed). That is unresolved, so it is NOT safe to rely on.
  - **288 ranks (144 per node, 36 per island).** This ran on g25_red, and so did 72 and 96 ranks.
- **Proposal, in one job:**
  1. B3 case3_dt0.1 x10 WITHOUT the env (the fence-fix check you asked for, ~3 min, 1 node of the 2);
  2. a 576-rank start check of g25_full (watchdog: no log growth for 180 s -> kill);
  3. the 3 arms at 576 ranks if the check passed, else at 288 ranks (2 x 144);
  4. MPICH_SHARED_MEM_COLL_OPT=0 for the arms unless the repro shows it is not needed.
- **Speedup caveat:** 288 ranks on 576 cores would understate the CPU arms. We would quote the per-core rates as well.

## 2. Queue, wall time, cost (sbatch --test-only, 10-26 09:53 CDT; balance bivj-dtai-gh 821 SU)
- **ghx4-interactive, 2 nodes, 1 h:** est. start in **~35 min** (10:28 CDT). Not crowded right now.
- **ghx4 (batch), 2 nodes, 1 h:** est. start **10-30 22:18 CDT (~4.5 d)**.
- **Wall time, estimated, not measured on g25_full:**
  - **Scaling basis:** the 4-node Raven result (GPU 20.94 / CPU 28.45 s per top step) and g25_red here (GH200 node ~2.1x a Raven node per IDORT full sync; Grace 288 cores vs 72).
  - **Per top step:** 2 DeltaAI nodes should be about like 4 Raven nodes, i.e. ~20-30 s GPU and ~30-45 s CPU.
  - **Per arm:** with init (+ first step), 5 top steps, ~5-8 min.
  - **Whole job:** 3 arms + the B3 check ~25-35 min, so a 1 h limit is fine.
- **Cost:** <= 16 SU (2 nodes x 1 h, interactive charged x2). The 10-25 jobs cost ~12 SU.

## 3. Blocking / memory
- **Code:** the missing prerequisite (above).
- **Inputs:** OK. g25_inputs.tgz md5 **2348f0db65cc644446240951bd518235** (matches) is unpacked; g25_full.hdf5 md5 40c5ce81 OK; tables OK.
- **Memory:** Raven 4 nodes peaked at 135 GB host and 117 GB GPU per node, so ~2x per node here.
  - **Host:** ~270 GB of 488 GB. MaxMemSize must come down from 3000 MB (3000 x 288 = 864 GB). We would use ~1200 MB at 288 per node, ~2500 MB at 144 per node.
  - **GPU:** ~60 GB per GH200 of 96 GB, so GpuMaxMemSize ~60000 (IDORT + VCELL + NLTREE on the feeder).
  - If the Raven peaks were not per node, please say.
- **Build:** ~10-20 min on the login node (no SU), the same recipe as 10-25 (sm_90, cray-mpich, local HYPRE/GMP). The nltree configs are new here, but no extra library is expected.
- **Rules:** no Arepo source in git; srun -K everywhere (MPI_Abort did not end the step on 10-25).

Waiting for the user's choice of site. If DeltaAI: the user copies the bundle (with the prerequisite) to `/work/nvme/bivj/jma20/arepo2_gpu_1025/bundle/`.
