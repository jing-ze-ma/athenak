# NOTE 2026-10-25 deltaai -> viper (USER): Arepo2-IDORT island sort on aarch64: MPICH_SHARED_MEM_COLL_OPT=0 fixes it; 10/10 gates ALL SAME; DeltaAI stopped (exception a1d0a3c2)

**Verdict.**
- **Cause:** the failures come from Cray MPICH's shared-memory collectives on GH200/aarch64.
- **Fix:** with `MPICH_SHARED_MEM_COLL_OPT=0` (or `MPICH_COLL_OPT_OFF=1`), the B3 case3 repro is clean, and all 10 gates are ALL SAME.
- **Not in git:** no Arepo source was changed or committed.
- **DeltaAI is stopped:** nothing is running or queued, no stage 2, no new jobs.
- **Job:** 3454795 (resubmit of 3454468, which I had cancelled on the withdrawal 598febf4 before the exception), 31 min, 1 node, ghx4-interactive.
- **Previous NOTE:** 0fd42724 (interim).

## 1. A/B repro: B3 case3_dt0.1, 4 ranks, 2 islands (CUDA_VISIBLE_DEVICES 0,1), IdortGpu 0, 10 runs per env, srun -K, 90 s timeout
| env | failed |
|---|---|
| default | **2 / 10**: both hang in `VCELL: create mesh (direct cell-by-cell construction)` (killed at 90 s) |
| `MPICH_SHARED_MEM_COLL_OPT=0` | 0 / 10 |
| `MPICH_COLL_OPT_OFF=1` | 0 / 10 |

**Reading.**
- **Two code paths fail:** with the default env, the failures are not only the island sort (job 3454194: the sort hang at 288 ranks and the TERMINATE at parallel_sort.h:1710). The VCELL mesh build fails too, and both use mymallocshared windows synchronised by NodeBarrier.
- **Cause:** turning off Cray's shared-memory collective implementations (MPI_Barrier/Bcast/... over the node's shared memory) removes both failures.
- **Consistent with the aarch64 ordering hypothesis of 0fd42724:** Cray's shared-memory barrier does not order the ranks' plain stores to the windows on Grace, while the point-to-point barrier does.
- **Suggestions for the code:**
  - a full fence or MPI_Win_sync at NodeBarrier in SHAREDNODEMEMORY builds;
  - or document `MPICH_SHARED_MEM_COLL_OPT=0` as required on Cray/aarch64.
- **Hypothesis 2 (aside):** the tag-formula remark of 0fd42724 was not needed to explain this.

## 2. Gates with MPICH_SHARED_MEM_COLL_OPT=0 (IdortGpu 0 vs 1, same binary): the 4 left from 3454194, all ALL SAME
| test | snapshots bitwise | solver rows equal | matvecs (g0 = g1) | device solves |
|---|---|---|---|---|
| B3_thermal_equilibrium (6 runs) | 102/102 | 816/816 | 73, 73, 466, 496, 2962, 1051 | 560 |
| A7_table_photosphere_rv | 4/4 | 8/8 | 694 (= Raven S0) | (VERIFY path only) |
| I1_restart_bitwise | 10/10 | 16+16 | 855 / 855 restart | 16 |
| I4_star_lite_fhyb_levels_m4 (3 runs) | 30/30 | 272/272 | 4519 / 11332 / 6208 | 232 |

- **Budgets:** identical in every run, incl. tsolve_its.
- **check.py:** PASS in both arms.
- **Full set:** with 0fd42724 (A1, A4, A4_m4, A6, A2, I2a), **10/10 ALL SAME on DeltaAI**.
- **I4_m4 matvec totals:** 4519 / 11332 / 6208, vs Caltech 4539 / 11368 / 6241. This is cross-machine round-off, as Caltech noted.

## 3. g25_red, GPU arm only (IdortGpu 1, IdortGpuMinCells 1e5, 4 islands, CPUsUsedPerGPU 1, 2 top steps, MaxMemSize 1000)
Ranks are spread evenly over the 4 Grace sockets (`-c 288/N`). Env: MPICH_SHARED_MEM_COLL_OPT=0.

| ranks (per island) | Code run for | 2 full syncs, step | per full sync | sweeps | gpu io | 44 partial steps | cpu.txt misc / treegrav / voronoi |
|---|---|---|---|---|---|---|---|
| 72 (18) | **76.7 s** | 3.68 s | **1.84 s** | 1.88 | 1.29 | 3.91 s | 60.8 / 4.3 / 7.5 |
| 96 (24) | 88.2 s | 3.43 s | 1.71 s | 1.76 | 1.22 | 3.43 s | 74.1 / 3.5 / 6.9 |
| 144 (36) | 119.2 s | 3.67 s | 1.83 s | 1.46 | 1.71 | 2.67 s | 106.2 / 2.5 / 6.3 |
| 288 (72) | - | - | - | - | - | - | stalled: no output after the banner for 120 s, killed by the watchdog |

**References (GPU arm, per full sync):** Raven 72 ranks 3.90 s; Caltech 32 ranks 4.17 s.
- **IDORT full sync:** a GH200 node does it in ~1.7-1.8 s, about 2.1x faster than Raven per node.
- **Wall time grows with rank count:** misc rises from 61 to 106 s, although the solver rows of the full syncs and partial steps are flat or falling. The extra cost is outside the timed solver phases (most likely setup/collectives that MPICH_SHARED_MEM_COLL_OPT=0 slows down at more ranks). That was not measured further.
- **No CPU arm:** none ran in this job (scan design: GPU arm first; the TASK was withdrawn before a CPU follow-up). So there is no DeltaAI CPU/GPU speedup number.
- **288 ranks:** not settled. It is either a slow start-up with the env, or another hang before the parameter read.

## 4. Cost and state
- **Jobs:** 3454194 (56 min) + 3454795 (31 min), interactive x2 charge: ~12 SU in total.
- **Nothing running or queued.**
- **Files (DeltaAI):** /work/nvme/bivj/jma20/arepo2_gpu_1025: builds/ (11 binaries), env_dtai.sh, Makefile.systype overrides, job1.sub, job2.sub, runs/3454194, runs/3454795, logs/.
