# NOTE 2026-10-20 deltaai -> raven: RSG design E queued on 2 ghx4 nodes (8 GH200), race copy (TASK cc6acfc6)

- Build: origin/fin51-r-fix @ 0fbfa31f (82275582 was not on the branch at build time), he_star_m1, nofma recipe
  (CUDA, Kokkos_ARCH_HOPPER90, Release, MPI ON, CMAKE_CXX_FLAGS=-ffp-contract=off; CMakeCache checked), cray-mpich + GTL.
  Binary /u/jma20/ATHENAK/builds/bin/athena_hesnf_gpu_0fbfa31f52d5, md5 cd328571a83eaede115bea3162e47c96.
- Inputs from bsg-files-1009 (cc6acfc6): deck rsg_E_start1020r.athinput.in md5 393d443a (only the 3 @FILES@ lines replaced),
  IC ic_rsg_g25_deep.txt 2d7de3c8, Rosseland fa0f4f42, Planck a9a8d6d3 (all match). FRESH start, tlim 1e9.
- Run dir /work/hdd/bivj/jma20/rsg_E_1020 (run/ = output, rst pruned in-job to 3 complete; link.sub, JOBS, logs/).
- No pre-start check job. Chain: 3423244 (link 1, 48 h) -> 3423245 (link 2, afterok). 2 nodes x 4 tasks x 4 GPUs, -A bivj-dtai-gh.
  Budget cap: each link trims its own TimeLimit so the chain spends <= 767 SU (867 now; ~100 SU kept) = ~96 h of 8 GH200.
  If I win, the run stops cleanly at the end of link 2 at an rst and I say so. STOP only on FATAL / DIVERGED / NON-CONVERGED=[1-9] / NaN / rc.
- Start estimate: sbatch --test-only for 2 nodes x 48 h says 2026-10-24 06:39 CDT (~3.7 d out; ~2390 jobs pending on ghx4);
  squeue --start N/A so far. DeltaAI is the long shot vs Raven's ~19 h. Tell me when another copy has started and I cancel (0 SU if pending).
- ghx4-interactive (not touched): MaxTime 2 h (default 30 min), MaxNodes 2 (= 8 GH200 per job), QoS qos_ghx4int MaxJobsPU 1 /
  MaxSubmitPU 1, billed x2 (GPU weight 2000 vs 1000 on ghx4). ghx4: MaxTime 2 d.
