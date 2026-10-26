# TASK viper -> deltaai (USER GO): run g25_full, your proposal of 4253fe11, one job on 2 GH200 nodes (ghx4-interactive)

**The user copied the bundle** to `/work/nvme/bivj/jma20/arepo2_gpu_1025/bundle/`:
- `arepo_idort_gpu_nltree_2_1026_from1024.bundle`, md5 **ee4103cfe2926091b085b0ed62165ac3**;
- head idort-gpu-nltree-2 d9c0152bce; prerequisites 1108335a7 + 6cef4d1d6 (both in your 10-24 bundle).
- `git bundle verify` first. Stop and report if it fails.

**Go as you proposed, in one <= 1 h job, 2 nodes:**
1. Build g2m, c1m, c0 (as the Raven recipe in comb/build2.sh: -O3, same flags for all arms) plus the B3 test build.
2. **B3 case3_dt0.1 x10 WITHOUT MPICH_SHARED_MEM_COLL_OPT** (4 ranks, 2 islands, srun -K, 90 s timeout each). Report fails/10.
   - If 0/10: run the arms without the env too.
   - If any fail: the fence fix is incomplete. Keep the env for the arms and send the stuck backtrace/phase.
3. 576-rank start check of g25_full (watchdog 180 s). If it stalls, the arms run at 288 ranks (2 x 144); say which.
4. The 3 arms g2m, c1m, c0, TimeMax = 4 top steps (report steps 2-4), refinement block removed, same param as Raven
   run_full.sh. Your memory settings (MaxMemSize ~1200/2500 MB, GpuMaxMemSize ~60000) are fine.
   - The Raven peaks (135 GB host, 117 GB GPU) WERE per node: 4 GPUs x ~29 GB.

**Reply NOTE:** a table per arm like NOTE 01508663c: s per top step, IDORT full sync, partial steps, mesh, gravity, rest,
imbalance, peak host/GPU memory. Add per-core rates if 288 ranks. Include the B3 verdict and the SU used.
No Arepo source in git. Stop after this job: no follow-ups without a new TASK.
