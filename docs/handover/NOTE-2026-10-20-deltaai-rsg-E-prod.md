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

## Update 10-20 12:37 CDT: USER override of step 5 -> ghx4-interactive chain (USER: "yes")
- ghx4-interactive test-only said 2 nodes x 2 h start in ~25 min (vs ~3.7 d on ghx4). USER asked to use it despite step 5
  ("do not touch ghx4-interactive"): there is no BSG chain of jma20 on DeltaAI (squeue shows nothing else), and the slot is per user.
- chain.sh (pid on gh-login03) submits 2 h interactive links one at a time (2 nodes, 8 GH200, x2 charge = 32 SU per link), each
  restarting from the newest rst (final rst written at the -t limit). First interactive link: 3423285 (pending, Resources).
- ghx4 links resubmitted with the run lock (old ones 3423244/45 cancelled pending, 0 SU): 3423286 -> 3423287 (afterok). When one
  starts the interactive chain stops; the ghx4 link waits on the lock (<= 130 min) and trims itself to the budget.
- Budget unchanged: all links together <= 767 SU (~100 SU kept). Interactive only = ~48 h of 8 GH200.
- I will post a start NOTE (s/cycle, dt, memory). Tell me if another copy started first and I cancel everything.

## Update 10-20 13:45 CDT (reply to 3188fe2e): switched to prodmerge-1019b 3896d529 BEFORE the first link started
- Built 3896d529 he_star_m1 GH200 nofma (same worktree/recipe, -ffp-contract=off checked): athena_hesnf_gpu_3896d52937d8,
  md5 698680999418c7886ec6015977278b7f. Every link from now on reads it (BINARY file); deck unchanged; still a FRESH start.
- First interactive try 3423285 (12:37) died at launch: node gh070 "Error configuring interconnect", srun hung 15 min (3.9 SU).
  link.sub now has --kill-on-bad-exit + a launch watchdog (bad node -> BADNODES, excluded on resubmit). Retry 3423395 (1 h limit,
  excludes gh070) is first in line, est start ~13:50 CDT; ghx4 links 3423396 -> 3423397 pending (est ~10-21 15:30 CDT).
- Race rule noted (DeltaAI start = head start, others keep running) pending USER.

## Update 10-20 13:50 CDT: 8 GH200 do NOT fit 3 blocks/GPU -> GPU OOM at setup; everything stopped, nothing queued
- Interactive link 3423395 (3896d529, gh[084,152]) launched fine, then died 45 s into setup:
  "Kokkos ERROR: Cuda memory space failed to allocate 621.3 MiB", right after
  "vet_gd: ragged per-shell band intensity array 60.7004 GB on rank 0". DeltaAI GH200 = 96 GB HBM (gres gh200_120gb);
  3 blocks x ~25-30 GB (TASK caltech estimate) does not fit; it fits H200 141 GB, not GH200. 0.3 SU; total spent 4.2 SU.
- ghx4-interactive caps a job at 2 nodes (8 GPUs), so the interactive route is closed for this layout. Interactive chain stopped;
  ghx4 links 3423396/97 cancelled pending (they would OOM the same way).
- Would fit: 3 ghx4 nodes (12 GPUs, 2 blocks/GPU, 12 SU/h -> ~64 h within budget) or 6 nodes (24 GPUs, 1 block/GPU, 24 SU/h -> ~32 h).
  test-only start for either: 2026-10-25 ~09:00 CDT (~4.8 d out). Waiting for the USER's call; my recommendation: DeltaAI stands down.
- Run dir kept (fresh, no rst): /work/hdd/bivj/jma20/rsg_E_1020.
