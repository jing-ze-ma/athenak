# NOTE Caltech -> Raven (2026-10-18): 3-D accretor production copy (TASK-2026-10-18-caltech-acc3d-prod), status 1

All jobs SUBMITTED (Caltech, account carnegie_poc, --exclude=hpc-sm-01-09). Link 1 has NOT started yet; the START time
follows in the next NOTE update.

- Build: accretor-3d-1017 @ 186306cb, PROBLEM=ry_per_accretor, CUDA HOPPER90, MPI, Release, -ffp-contract=off (nofma),
  nvcc_wrapper, gcc 13.2 + cuda 12.9 + hpcx 2.17.1 (git-archive snapshot build, bsg_1009/build_inc_nofma.sh,
  BUILD_ROOT builds_accr). Two build jobs (4970290 expansion, 4970319 gpu, whichever runs first; lock-protected).
  Binary: /resnick/groups/carnegie_poc/jingze/builds_accr/bin/athena_gpu_ry_per_accretor_186306cb_nofma (md5 follows).
- Files: /resnick/groups/carnegie_poc/jingze/acc3d_prod_1018/files, all 5 md5s match the TASK (e9dab743, d81b7a46,
  99d668f6, 5b813886, 1ad77aed); deck acc3d_prod_mk.athinput (md5 9748aeaa) = .in with only the 3 @FILES@ paths substituted.
- Layout: 1 node x 4 H200, 4 ranks (srun --mpi=pmix -n 4 -c 8 --cpu-bind=cores gpuwrap.sh, CUDA_VISIBLE_DEVICES=SLURM_LOCALID,
  UCX_TLS=^xpmem), deck MeshBlocks unchanged (8 MBs, 2 per GPU).
- Pre-check 4970356 (A nlim 10 / B restart to 20 / C fresh 20, sect. 4) afterok on the build; it exits 1 unless
  B vs C bins differ in <= 1 byte, hydro.hst rows identical and user.hst cols 1-6,16-20,23 identical.
- Chain: link 1 4970357 afterok:4970356, link 2 4970358 afterok:4970357, both --kill-on-invalid-dep=yes,
  walltime 7-00:00:00 (normal QOS max; WT = 603600 s), acc3d_chain.sh unchanged; job rc = chain rc.
- Run path (rsync-able): /resnick/groups/carnegie_poc/jingze/acc3d_prod_1018 (static/, stream/, chain.txt, starts.txt).

## Update 2 (2026-10-18)
- Build 4970290 COMPLETED (5 min): binary md5 **9283c619823898f3d0196a3abca6129f**
  (/resnick/groups/carnegie_poc/jingze/builds_accr/bin/athena_gpu_ry_per_accretor_186306cb_nofma); 4970319 cancelled (redundant).
- Pre-check 4970356 TimeLimit cut to 40 min (backfill). Not started yet.
- USER GO: parallel **2-H200 variant** (DEVIATES from the TASK 4-GPU layout: 2 ranks, 4 MeshBlocks per GPU, same binary,
  same deck md5 9748aeaa) in /resnick/groups/carnegie_poc/jingze/acc3d_prod_1018/g2: pre-check 4972809 (40 min),
  link 1 4972811, link 2 4972812 (7 d, afterok + kill-on-invalid-dep). Whichever variant's link 1 STARTS first runs;
  the other set is cancelled at once. Start time + layout follow in the next update.

## Update 3 (2026-10-18): pre-check results; links cancelled per HOLD b592b459, waiting for the new TASK
All 4 links (4970357, 4970358, 4972811, 4972812) were cancelled per HOLD b592b459 (no production on 186306cb); none had started.
The 4-GPU pre-check 4970356 had not started when the 2-H200 one finished, so it was cancelled. Only the 2-H200 result exists.

**2-H200 pre-check 4972809** (hpc-sm-01-05, 2 ranks, 4 MBs/GPU, binary md5 9283c619, deck 9748aeaa), run dir
/resnick/groups/carnegie_poc/jingze/acc3d_prod_1018/g2/pre:
- rc 0 for A/B/C; FATAL 0, summary NON-CONVERGED=[1-9] 0, Picard DIVERGED 0, DT_RETRY 0, NaN 0.
- dt 1.419254e-05 (t = 2.838468e-04 at cycle 20). Init ~44 s. Loop time A 13.6 s / 10 cycles, C 18.5 s / 20 cycles
  -> ~0.49 s/cycle after warm-up (0.93 s/cycle averaged over 20 cycles incl. the first ones), static phase, stream off.
- GPU peak memory 14661 MiB on each of the 2 GPUs (of 141 GB), so 4 MBs per H200 fits easily.
- B vs C: final ryper.hydro_w.00001.bin and ryper.m1.00001.bin match byte for byte (0 bytes differ). hydro.hst B rows (20)
  match C rows 2-21; user.hst cols 1-6,16-20,23 match too.
- **Caveat, the restart test is WEAKER than intended:** A wrote 2 rst files here, ryper.00000.rst at t=0 and ryper.00001.rst
  at its end (cycle 10). B used ryper.00000.rst as the TASK says, so it restarted at cycle 0, not at cycle 10. The job
  script reported fail=1 only because B does not rewrite the t=0 bins (my "MISSING" check, a script artifact). A proper
  cycle-10 restart (B from ryper.00001.rst) is a ~3-min 2-H200 job; I can run it when a new TASK asks for it.

## Update 4 (2026-10-18): PRODUCTION v3 SUBMITTED (TASK fe694be5)
- Build: accretor-3d-stream-1018 @ **a6624af6** (verified = origin), PROBLEM=ry_per_accretor, HOPPER90, -ffp-contract=off,
  job **4977640** (expansion); binary /resnick/groups/carnegie_poc/jingze/builds_accr/bin/athena_gpu_ry_per_accretor_a6624af6_nofma
  (md5 follows when built).
- Run dir /resnick/groups/carnegie_poc/jingze/acc3d_v3_1018: deck v10 5e342048 (@FILES@ sed only), chain 1ad77aed, rosseland
  d81b7a46, planck 99d668f6, column 5b813886 (all verified). Layout: 2 H200, 2 ranks, 4 MB/GPU, gpuwrap (as g2).
- Links (12 h, WT 41999, afterok + kill-on-invalid-dep; no pre-check): link 1 **4977641** (afterok:4977640 build),
  links 2-6 **4977642 4977643 4977644 4977645 4977646**. Job rc = chain rc. Start time + startup lines follow.

## Update 5 (2026-10-18): binary built
- Expansion queue estimate was ~10 h, so a backup build **4977664** (gpu partition, P100 node, compile only, HOPPER90 target)
  ran instead: COMPLETED, 2 files recompiled, binary md5 **76f5201569a523d7071d0525ee24c34f**
  (/resnick/groups/carnegie_poc/jingze/builds_accr/bin/athena_gpu_ry_per_accretor_a6624af6_nofma). Expansion build 4977640 cancelled.
- Link 1 4977641 dependency satisfied, pending on resources (Slurm estimate in ~6 h; backfill may be earlier).

## Update 6 (2026-10-18): LINK 1 STARTED
- Link 1 **4977641 START Sun Oct 18 13:48:19 PDT 2026 (= 20:48:19 UTC = 22:48:19 CEST)** on hpc-sm-02-10, 2 H200; static phase began 13:48:22 PDT.
- Startup lines (static/run.log):
  - `ry_per_accretor: stream_sector_light surface (back-sink sector radiation)`
  - `ry_per_accretor: env_ic_vz: von Zeipel IC T x (g/g_col)^(w/4), taper depth 1131 (psi), factor 0.94676 .. 1.00383 (rank 0)`
  - NOTE: the low end 0.947 is below the TASK's expected ~0.97-1.03 (rank 0 only is printed). Not a stop condition; flagged for review.
- Static phase healthy so far: cycle 250, t 3.53e-3, dt 1.407e-05, ~0.43 s/cycle; FATAL/NaN/DIVERGED 0.

## Update 7 (2026-10-18): QUESTION FOR RAVEN - env_ic_vz factor low end 0.947
- Caltech link 1 startup line: `ry_per_accretor: env_ic_vz: von Zeipel IC T x (g/g_col)^(w/4), taper depth 1131 (psi), factor 0.94676 .. 1.00383 (rank 0)`.
  The low end 0.947 is below the TASK's expected ~0.97-1.03.
- **Raven: please compare with your A100 v10 startup line.** Our layout is 2 ranks (4 MeshBlocks per rank), so "rank 0" covers
  a different, larger part of the mesh than your rank 0 (8 ranks x 1 MB); the min/max may differ for that reason alone. Please say
  whether we should continue. **We keep running unless told otherwise.**
- Static phase status: t 0.0216 at cycle 1550, ~0.44 s/cycle, dt 1.23e-5, 0 DT_RETRY, GPU memory 14661 MiB per GPU; FATAL/NaN/DIVERGED 0.

## Update 8 (2026-10-18): QUESTION FOR RAVEN - static-phase dt falls: mid-latitude infall, not breathing decay
- Static-phase dt (static/run.log dt_report) is NOT a breathing decay. A mid-latitude infall in the thin atmosphere grows:

  | t | dt | v_r at limiting cell |
  |---|---|---|
  | 0 | 1.42e-5 | |
  | 0.0214 | 1.22e-5 | -93 km/s |
  | 0.0392 | 1.01e-5 | -118 km/s |
  | 0.0440 | 7.1e-6 | -164 km/s |
  | 0.0480 | 6.0e-6 | -221 km/s |

- Limiting cells: r 9.07-9.3 (i 124-174), rho 2e-7..2e-6, theta 70-76 and 102-109 deg (both hemispheres, mostly rank 1),
  e/rho ~2000 (km/s)^2, flat. 0 DT_RETRY, NON-CONVERGED 0.
- Hypothesis (NOT verified): the cooler von Zeipel off-equator atmosphere (env_ic_vz factor down to 0.947) loses pressure support and sags.
- **Raven: please post your A100 v10 static dt history and limiting cells. Do you see the same infall? Should we continue?
  We keep running unless told otherwise.**
- Stream-on status: not yet. Static phase at cycle 4450, t 0.04932, dt 5.90e-6 (elapsed 1927 s); stream phase starts at t 0.05.
