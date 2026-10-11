# NOTE DeltaAI -> raven: AG Car A face/icrad (TASK-2026-10-11-deltaai-agcarA-face-icrad) -- DONE

Picked up 2026-10-10 21:30 CDT; result 21:58 CDT. **All 11 runs rc 0, FATAL 0, NON-CONVERGED 0.**

## Setup
- Job **3357130**, ghx4-interactive gh031, 1 node x 4 GH200, 1 rank per GPU (4 blocks 480x64x64), plain srun, wall 7 min 52 s.
  (First attempt 3357111 failed at startup on every run: my job script lacked MPICH_GPU_SUPPORT_ENABLED=1 -> Cray MPICH
  `process_vm_readv: Bad address`. Environment error, not code; cancelled after 2 min.)
- A: deepatm-1011 e84d2513, `athena_hes_gpu_e84d2513d2aa` md5 1e0134f8f8dd23669954dec55d1e469d.
  B: icrad-1011 7045bc14, `athena_hes_gpu_7045bc14cdba` md5 ddbb4dfce537ac2cd9d89665d5d8368c.
- Inputs: all four md5s match. Every key set in a per-run copy of agcar_A_va.athinput (none on the command line):
  ndiag 1, nlim, implicit_cfl, vet_gd_rotate_every, vet_bc_diff / ic_rad (+ ic_rad_relax_n 40, ic_rad_relax_tol 1e-4) added to <rad_m1>.
- **B (7045bc14) does not contain vet_bc_diff** (e84d2513 is not its ancestor): all B runs use the old cell boundary.
- dt = 1239.996 s (implicit_cfl 1251.3); dt/2 = 619.998 s (625.65, rotate_every 20).

## Every run
Picard = the code's end-of-run line (it counts the ic_rad relax solves as well: 74 solves for 60 cycles + 14 relax steps);
s/cycle = median (mean) of per-cycle wall times from `elapsed=`, excluding rotation cycles and the last.

| run | cycles | rc | FATAL | NON-CONV | Picard mean / max | s/cycle median (mean) |
|---|---|---|---|---|---|---|
| A cell | 50 | 0 | 0 | 0 | 3.90 / 9 | 0.199 (0.271) |
| A face | 50 | 0 | 0 | 0 | 3.90 / 9 | 0.199 (0.270) |
| A face repeat | 50 | 0 | 0 | 0 | 3.90 / 9 | 0.198 (0.270) |
| B default | 60 | 0 | 0 | 0 | 3.78 / 9 | 0.211 (0.272) |
| B fs | 60 | 0 | 0 | 0 | 3.72 / 8 | 0.207 (0.251) |
| B relax | 60 | 0 | 0 | 0 | 3.69 / 9 | 0.203 (0.245) |
| B fs_relax | 60 | 0 | 0 | 0 | 3.66 / 10 | 0.206 (0.251) |
| B default dt/2 | 120 | 0 | 0 | 0 | 3.35 / 8 | 0.198 (0.236) |
| B fs dt/2 | 120 | 0 | 0 | 0 | 3.14 / 8 | 0.194 (0.219) |
| B relax dt/2 | 120 | 0 | 0 | 0 | 3.04 / 9 | 0.192 (0.215) |
| B fs_relax dt/2 | 120 | 0 | 0 | 0 | 3.04 / 8 | 0.192 (0.215) |

## Run 1: face vs cell, cycle 50 (sin theta-weighted shell means; max |face-cell|/|cell| in the band, at r/R_ph)
| quantity | < 0.8 R_ph | 0.95-1.05 R_ph | > 1.05 R_ph |
|---|---|---|---|
| E | 0 (bitwise) | 8.88e-4 (0.977) | 2.9e-7 |
| F1 | 1.4e-11 | 7.1e-7 | 2.9e-7 |
| dens | 0 (bitwise) | 2.42e-4 (0.977) | 1.9e-7 |
| velx | 9.2e-10 | 1.97e-2 (0.984; velx near zero there) | 9.3e-7 |

- L_top: max rel diff 4.97e-6 (cycle 3); cycle 50: cell 2.844094e38, face 2.844093e38.
- **Repeat face run is bitwise**: 12/12 bin files and both hst files identical.
- Cost identical (0.199 s/cycle both).
- Agrees with Caltech's reduced CPU run (NOTE-2026-10-11-caltech-cpu-agcarA-face: E 8.9e-4 at 0.977, dens 2.3e-4, velx 1.8e-2,
  L_top 5.0e-6 at cycle 3) -> the effect does not depend on resolution in theta/phi and is the same on GPU.

## Run 2: IC radiation start (binary B, cell boundary)
L_top (erg/s) at fixed dt, cycles 0-10 (full table: deltaai-agcarA-1011/analysis.txt):

| cycle | default | fs | relax | fs_relax |
|---|---|---|---|---|
| 0 | 2.843095e38 | 2.843095e38 | 2.844903e38 | 2.844868e38 |
| 1 | **4.365807e38 (+53.6 %)** | 3.093375e38 (+8.8 %) | 2.844750e38 | 2.844724e38 |
| 2 | 4.107854e38 | 3.017300e38 | 2.844658e38 | 2.844638e38 |
| 3 | 3.551379e38 | 2.943511e38 | 2.844589e38 | 2.844572e38 |
| 5 | 2.957766e38 | 2.870911e38 | 2.844466e38 | 2.844455e38 |
| 10 | 2.847022e38 | 2.846799e38 | 2.844264e38 | 2.844258e38 |
| 30 | 2.844217e38 | 2.844213e38 | 2.844063e38 | 2.844061e38 |
| 60 | 2.844085e38 | 2.844084e38 | 2.844066e38 | 2.844065e38 |

- default: the +54 % step-1 spike, decays over ~10 cycles. fs: spike cut to +8.8 % (fs ends with max|<E>/<J>-1| 0.15).
  relax / fs_relax: **no spike**, L_top flat to 6e-4 from cycle 0; all four reach the same L_top within 7e-6 by cycle 60.
- **dt vs dt/2**, max (rms) of |L_top(dt) - L_top(dt/2)| / L_top(dt/2) over cycles 1-30 (dt/2 matched at cycle 2c):
  default 8.2e-2 (2.4e-2); fs 9.4e-3 (2.8e-3); relax 2.65e-4 (8.3e-5); **fs_relax 2.63e-4 (8.3e-5)**.
- **Best = fs_relax** (relax is the same to 1 %; fs_relax ends its fs pass with max|<E>/<J>-1| 1.2e-2 vs 0.15 for fs alone).
- Cost of the IC step (code's own timer): fs 3 passes 0.41 s; relax 14 steps 3.8 s (dt) / 20 steps 4.8 s (dt/2);
  fs_relax 4.3 s (dt) / 5.5 s (dt/2). All four relax / fs_relax runs stopped on relax_tol
  (|dL/L| < 1e-4 after 14 steps at dt, 20 at dt/2), well inside relax_n 40. One-off cost = ~19-28 cycles of run time.
- Picard falls with relax/fs_relax (dt/2: 3.04 vs 3.35 default) (likely because the early transient is gone).

## GPU-h
3357130: 7.9 min x 4 GH200 = 0.52 GPU-h; failed 3357111: 0.14 GPU-h. **Total 0.66 GPU-h.**

## Files
In git: docs/handover/deltaai-agcarA-1011/ (analysis.txt, job.sub, analyse.py). On DeltaAI:
/work/nvme/bivj/jma20/agcarA_1011: job.sub, analyse.py, runs/3357130/<run>/{deck.athinput,out.txt,*.hst,bin/},
logs/agcA_face_icrad.3357130.out.
