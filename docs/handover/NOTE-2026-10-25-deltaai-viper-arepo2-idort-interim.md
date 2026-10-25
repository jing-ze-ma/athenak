# NOTE 2026-10-25 deltaai -> viper (USER): Arepo2-IDORT GPU on DeltaAI, INTERIM (TASK fc5cf8f2): build OK, 6/10 gates ALL SAME, island parallel sort fails on aarch64

**Verdict so far.**
- Build: OK.
- Node facts: OK.
- Gates: 6 of 10 ALL SAME.
- g25_red smoke at 288 ranks: **no timings**. Both arms hang in init in the island parallel sort.
- B3 (4 ranks, 2 islands): the same sort TERMINATEs intermittently, and no other gate fails.
- **Request to viper:** the sort is Arepo code (not IDORT_GPU), and nothing was patched here. Please look at it. Details below.
- **Not in git:** no Arepo source, diffs or patches were committed anywhere.
- **Job:** 3454194 (ghx4-interactive, 1 node, 56 min). The follow-up 3454468 is queued, est. start 10-25 ~18:00 CDT.

## 1. Node facts (gh090, inside the job)
- **Hardware:** 1 node = 4 x GH200 (96 GB HBM each) + 288 Neoverse-V2 cores (4 x 72, one Grace per GPU); 488 GB host.
- **Software:** cray-mpich 9.0.1 (PrgEnv-gnu, gcc 14.2), CUDA 12.9, driver from nvidia-smi in the job log.
- **Peer access:** cudaDeviceCanAccessPeer = 1 for all pairs.
- **CUDA IPC between processes WORKS:** 1 rank per GPU, cudaIpcOpenMemHandle of a rank-0 buffer, 127 GB/s copy GPU0 -> GPUk.
- **Bundle:** caltech_bundle_1024. MD5SUMS all OK; bundle md5 905ec3cb6ad069aa1e491658956b3da8. The g25 MD5SUMS (incl. g25_full.hdf5) are also OK.

## 2. Build adaptations (idort-gpu-s1 6cef4d1d6, all 10 gate configs + g25g)
- **Makefile.systype** (untracked): SYSTYPE Raven-OpenMPI with `override` lines:
  - CUDA `-arch=compute_90 -code=sm_90` in the compile and the nvcc link line;
  - MPICHLIB `-ldl -lz`: -lbfd dropped (BACKTRACE is unused);
  - local HYPRE and GMP paths.
- **Modules:** cray-mpich through the mpicxx wrapper; cray-hdf5 1.14.3 (serial), gsl 2.8, cray-fftw 3.3.10.
- **Built here:**
  - HYPRE 2.32.0 (CPU, static, from recipes.tgz);
  - **GMP 6.3.0**: DeltaAI has no gmp.h (runtime only). Static.
- **Flags:** OPTIMIZE `-g -march=native -O2 -flto` for the gate builds, `-O3` for g25g. -march=native is fine: the login node is Neoverse-V2 like the compute nodes. -ffp-contract=off and -fmad=false are in-tree and untouched.
- **Python:** venv with Python 3.11, numpy 1.26.4, h5py 3.16, scipy 1.17.
- **mpiexec:** DeltaAI has none. test.py's mpiexec is a shim to `srun --overlap -K -n N -c 1 --cpu-bind=cores`.
- **Test copies:** the idort_s0 copies (A4_shadow_vcell_m4, I4_star_lite_fhyb_levels_m4) come from port/make_s0_tests.sh.
  - A7_table_photosphere_rv: copied with the table paths set to the bundle tables. create.py and check.py get a prepended matter_tables path header; OpalEos5's default argument is patched too.
- test.py's make in a run dir only regenerates version.cc and relinks (as on Raven).
- **No defaults, solver settings or preconditioner were changed.**

## 3. Gates (IdortGpu 0 vs 1, same binary, 2 islands = CUDA_VISIBLE_DEVICES 0,1; I2a 1 island; AREPO_IDORT_GPU_VERIFY=1)
cmp_g1.py: ALL SAME for the 6 tests that ran.

| test | snapshots bitwise | solver rows | matvecs (g0/g1) | budget |
|---|---|---|---|---|
| A4_shadow_vcell_m4 | 2/2 | equal | 455/455 (= Raven) | identical |
| A1_crossing_beams_vcell | 2/2 | equal | 47/47 | identical |
| A4_shadow_vcell | 2/2 | equal | 53/53 | identical |
| A6_photosphere | 2/2 | equal | 44/44 | identical |
| A2_sphere_3d | 2/2 | equal | 78/78 | identical |
| I2a_rank_invariance (np1/np4/np8) | 6/6 | equal | 44, 32, 47 | identical |

- **GPU path ran:** the g1 logs say "IDORT: GPU island solve: 2 islands, 4 to 4 ranks per island".
- **VERIFY lines are 0,** as on Caltech: verify is only on the IdortGpu 3 path.
- check.py passes in both arms.

**Not done.**
- **B3_thermal_equilibrium:** the g0 arm stopped in check.py's sub-run case3_dt0.1 (see 4b), so there is no g1 arm.
- **A7_table_photosphere_rv, I4_star_lite_fhyb_levels_m4, I1_restart_bitwise:** out of time. They are in job 3454468.

## 4. Island parallel sort (mycxxsort_parallel_shared, src/sort/parallel_sort.h) on DeltaAI
**a) g25_red, 288 ranks = 4 islands x 72, both arms.**
- **Where:** init -> test_id_uniqueness -> mycxxsort_parallel -> mycxxsort_parallel_shared, after "PARALLEL SORT: Starting global sort.".
- **Hang:** no progress for 15 min (CPU arm, killed by timeout) and >2 min (GPU arm, cancelled). The ranks spin in MPI.
- **gdb backtraces (every 12th rank):**
  - ranks 0 and 72 (island masters): in MPI_Sendrecv of the per-ngrp header exchange (parallel_sort.h:1573);
  - their other island ranks: in NodeBarrier (:1583);
  - island 2/3 ranks: in the data MPI_Sendrecv (:1683) or NodeBarrier (:1554).
  - So the islands are in different ngrp rounds, and some pairs wait for partners that skipped the exchange.

**b) B3 case3_dt0.1 (1024 cells, 4 ranks, 2 islands, IdortGpu 0).**
- **Failure:** `TERMINATE ... mycxxsort_parallel_shared(), parallel_sort.h line 1710: local_count=1, ...` in test_id_uniqueness.
- **Intermittent:** the other three B3 sub-runs passed the same sort 5 times each.
- **Second problem:** MPI_Abort did not end the srun step. The other 3 ranks stayed until the job ended (24 min). srun -K is now used everywhere.

**Reading.**
- **Same code:** the sort code is the same on Arepo2-IDORT-dev. Raven (x86, 18 ranks per island) and Caltech (x86, 8 per island) pass B3 and g25_red.
- **Hypothesis 1 (not proven):** aarch64 weak memory ordering. The sort exchanges its counters through mymallocshared windows (FreeElementsNode, SendNode, ...), separated only by NodeBarrier = MPI_Barrier. Cray MPICH runs that barrier in shared memory (MPIDI_Cray_shared_mem_coll_bcast in the backtrace) and may not act as a full fence for the window memory on Grace. On x86 (TSO) the same pattern happens to work. If so, an MPI_Win_sync (or a full fence) around the NodeBarriers of the shared sort would fix it.
- **Hypothesis 2 (aside, not shown to cause this):** the data tag is formed from ngrp times ptask, not ngrp times 2^ptask, so tags can repeat across rounds. No same-pair collision was found by reading.

**Job 3454468 (queued) tests the first hypothesis without a code change.**
- **A/B:** B3 case3 x 10 per MPICH env: default, `MPICH_SHARED_MEM_COLL_OPT=0`, `MPICH_COLL_OPT_OFF=1`. It takes the first clean one.
- **Then, with the chosen env:**
  - g25_red at 288 ranks (GPU then CPU arm);
  - the 4 remaining gates;
  - a g25_red layout scan (144 / 72 / 96 ranks).
- **Questions for viper:**
  1. Is there a known issue with the island sort on ARM / Cray MPICH?
  2. Should a Win_sync / fence fix go into the bundle for stage 2?
  3. Or should DeltaAI run with a non-shared-memory barrier env, if 3454468 shows that one works?

## 5. Cost
Job 3454194: 56 min on 1 node in ghx4-interactive (charged x2): ~8 SU. Balance before it: 833 SU.
Files: /work/nvme/bivj/jma20/arepo2_gpu_1025 (job1.sub, job2.sub, logs/, runs/3454194/).
