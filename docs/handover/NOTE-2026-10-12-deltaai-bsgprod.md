# NOTE DeltaAI: REAL BSG true reproduction, final method, chained 1 h links (TASK-2026-10-12-deltaai-bsgprod) -- RUNNING (link 2 queued)

**Picked up** 11:16 CDT (TASK 40558c5d, user GO). bsgcand-1012 **d3bf42d7** built: `athena_hes_gpu_d3bf42d7dc27`,
md5 63528bf7f170366e4fd5cbc05f4b74b6 (deck keys implicit_thin_js_face, vet_gd_rotate_slope, vet_gd_halo_ipc, vet_gd_twin_lowmem present).
- Files: the 10-09 BSG files (md5 ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805). Deck bsgP_seed.athinput.in
  (md5 e193df54, after the 53995cf9 fix) with @FILES@ resolved only.
- Layout 1 node x 4 GH200, 4 ranks x 4 MeshBlocks, plain srun, MPICH_GPU_SUPPORT_ENABLED=1.
- Link 1 = job **3367085** (ghx4-interactive, 1 h): 25-cycle smoke (nlim 25, ndiag 1, no rst; rc / FATAL / NaN / NON-CONVERGED -> STOP),
  then FRESH production with `-t` = wall left - 5 min. ic0check run on the smoke dumps from the login node while production runs.
- Later links: newest rst + rst_info.py outputN/last_time keys; STOP on rc != 0 / FATAL / NaN / DIVERGED / NON-CONVERGED > 0; DONE at t >= 4.96e6.
  A login-node chain submits each link when the previous one ends (1 job per user on ghx4-interactive); it stops if a bsgprod
  CANCEL / WITHDRAWN NOTE appears. Cap 12 links.
- Dir /work/nvme/bivj/jma20/bsgprod_1012 (run/, smoke/, logs/link.<job>.out). STARTED line follows as soon as link 1 runs.

**STARTED**: link 1 job **3367085** RUNNING on gh016 since 11:23:07 CDT (smoke 25 cycles, then fresh production in the same job). Smoke numbers and ic0check follow.

## Smoke (link 1, job 3367085, gh016, 11:23:07-11:23:35 CDT): PASS -> production started fresh at 11:23:35 (`-t 00:54:35`)
- rc 0, FATAL 0, NaN 0, NOT ADMISSIBLE 0, DIVERGED 0; `implicit transport: solves=50 Picard iterations mean=2.68 max=7 NON-CONVERGED=0`.
- dt = 87.5 s on every cycle (26 of 26 lines), 25 cycles to t = 2187.5 s. s/cycle 0.421 (cycles 2-24), init 2.7 s.
- GPU peak 50.3 GB on each of the 4 GH200 (4 blocks/GPU, lowmem on).
- hesdirk2: stage steps 24, BE 1, fallbacks 0, vet clips 0, `vet_gd_rotate_slope=correct solves=1 failed=0`.
- ic0check (`ic0check.py check smoke ic0_ref_1009.json`):
  ```
  dens     BITWISE
  velx     BITWISE
  vely     BITWISE
  velz     BITWISE
  eint     differs: max rel shell-mean 4.729e-08 sum 9.151e-13 min 0.000e+00 max 0.000e+00
  m1_e     differs: max rel shell-mean 4.579e-08 sum 1.002e-11 min 0.000e+00 max 0.000e+00
  m1_f1    BITWISE
  m1_f2    BITWISE
  m1_f3    BITWISE
  ```
  The eint / m1_e shell means differ only in shells 64-112 = **37.62-45.12 Rsun** (cell centres), i.e. exactly the seeded band;
  ic0check reports shell means (angle-averaged), so the ~2e-6 / ~9e-6 per-cell seed shows as ~5e-8 here.

## Links
| link | job | wall (CDT) | t end [s] (d) | cycle end | dt | s/cycle | Picard mean / max | FATAL / NaN / NOT ADM / DIVERGED / NON-CONV | slope solves (failed) | newest rst |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 3367085 | 11:23:35-12:18:17 (-t 54:35) | 5.317375e5 (6.15 d) | 6077 | 87.5 s throughout | 0.540 (cycles 100-6050, incl. dump cycles) | 3.34 / 13 (12456 solves) | 0 / 0 / 0 / 0 / 0 | 303 (0) | rst/bsg3d.00007.rst |

Link 2 = job 3367462, submitted 12:18:47 by the chain (restart from bsg3d.00007.rst + rst_info last_time keys).
Projection at 0.54 s/cycle: ~56.7k cycles = ~8.5 h compute = ~9-10 links to tlim 4.96e6 s.
Page: the BSG True Reproduction artifact now shows this run (t = 4.1 d at publish; refreshed from DeltaAI with the Caltech pagework, 822df54e).

- 12:30 CDT: links shortened to 40 min (code -t 35 min) to backfill sooner; the 1 h link 2 (3367462) was cancelled while pending.
- 12:58 CDT (user): also a 12 h ghx4 batch job **3367703** for the remainder (est. 10-13 afternoon). link.sub now takes an atomic run lock
  (run/LOCKD, owner job id; a job that finds the run held by a running job exits at once); the chain waits while the ghx4 job
  holds the run and cancels it if still pending at DONE.
- **Link 2 STARTED**: job **3367704** (ghx4-interactive, 40 min) on gh081 since 15:41:44 CDT, holds the lock, restarted from
  rst/bsg3d.00007.rst (t 5.317e5 s, cycle 6077); cycle 6100 reached at 15:42.

| link | job | wall (CDT) | t end [s] (d) | cycle end | dt | s/cycle | Picard mean / max | FATAL / NaN / NOT ADM / DIVERGED / NON-CONV | slope solves (failed) | newest rst |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 3367704 | 15:41:49-16:17:05 (-t 35:00) | 8.101625e5 (9.38 d) | 9259 | 87.5 s throughout | 0.657 (cycles 6100-9250) | 4.42 / 12 (6523 solves) | 0 / 0 / 0 / 0 / 0 | 159 (0) | rst/bsg3d.00011.rst |

Link 3 = job 3368862 (40 min) submitted 16:17:47 by the chain; ghx4 3367703 still pending. Note: s/cycle 0.54 -> 0.66 and Picard mean 3.34 -> 4.42 from link 1 to link 2 (the seeded lateral motion growing; still 0 NON-CONV).
