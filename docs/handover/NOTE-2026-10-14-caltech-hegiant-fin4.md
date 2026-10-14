# NOTE Caltech -> Raven: He giant FIN4 (fin4-1014 @ d505d21f): smoke gate BITWISE, link 1 RELEASED per GO; BG skipped per GO

Reply to TASK-2026-10-14-caltech-hegiant-fin4.md and GO-2026-10-14-caltech-hegiant-fin4.md (ca1d5d65). Work dir
`/resnick/groups/carnegie_poc/jingze/hegiant_1014`.

## 1. Build
- `athena_gpu_he_star_m1_d505d21f_nofma`, d505d21f (on origin/fin4-1014), same toolchain/flags as 9972a56d/dad00933
  (Kokkos CUDA HOPPER90, -ffp-contract=off). **md5 1440c05320a7bef88c4ad960019b9840**.

## 2. Deck
- `hegiant_fin4.athinput.in` md5 1b4d0592 (checked) -> sed @FILES@ = `/resnick/groups/carnegie_poc/jingze/hegiant_1009/files`
  (ic_giant_own.txt ddc72692, rosseland c1ef063c, planck b699c60d, all checked); no key changed; read-only.
  **hegiant_fin4.athinput md5 f266cd9bc230fdc805a48cdb9dcabfe5**.
- Layout as FIN3: 1 node, 2 ranks x 2 MeshBlocks, 1 rank per H200 (CUDA_VISIBLE_DEVICES=$SLURM_LOCALID wrapper), --cpu-bind=cores.

## 3. Validation (job 4489483, hpc-sm-02-16)
- **A smoke C** (fresh, nlim 24, ndiag 1, bin dump every cycle, rst 00001 at t = 240 s = cycle 12): rc 0; FATAL 0, NaN 0,
  DIVERGED 0, NON-CONVERGED 0. Startup line `time2_enth_be=true implicit_cs2_vimp_npass=3` present.
  Picard mean 5.33 / max 10 (88 solves). dt 20 s on all cycles. s/cycle (cycles 3-24) 3.06 (smoke with a bin dump + diag
  every cycle; FIN3 dad00933 same smoke: 2.74). GPU peak memory 85749 MiB (A/B, same as FIN3).
- **B restart R** (new srun, `-r A/rst/hegiant.00001.rst time/nlim=24`, fresh dir): rc 0; FATAL/NaN/DIVERGED/NC 0;
  Picard mean 5.00 / max 6 (25 solves). **BITWISE vs C**: final bin dumps (cycle 24, t 480) hydro_w (dens velx vely velz eint),
  m1, m1_vet: data after the header SAME; user.hst and hydro.hst: all 12 rows after cycle 12 verbatim in A.
- **C background run BG: SKIPPED per GO** (validated on Raven/Viper); the validation job was cancelled at the start of arm C.
  No budget / band numbers from Caltech.

## 4. Production chain (`hegiant_1014/fin4`, fresh t = 0 -> tlim 2.592e6 s)
- link_hg4.sh (-J hgF4, FIN3 link pattern: run.cfg at job start, FRESH if no rst else newest rst + rst_info.py keys,
  STOP on rc != 0 / FATAL / NaN / DIVERGED / NON-CONVERGED > 0, DONE at rst t >= tlim, prune_rst.sh keep 3 EXIT trap);
  run.cfg BIN d505d21f + md5 check, deck hegiant_fin4.athinput.
- 4 x 12 h links, afterany: **4489484 (link 1), 4489485, 4489486, 4489487**. Link 1 was submitted HELD and RELEASED after
  the smoke gate passed (GO ca1d5d65).
- BSG Linhao chain (bsgli_1014) and the stopped FIN3 run (hegiant_1013/fin3) untouched; all data kept.

## 5. STARTED
- STARTED: link 1 job 4489484 RUNNING on hpc-sm-02-16 since 2026-10-13T21:42:06 PDT, FRESH t = 0 (`-i hegiant_fin4.athinput`, -t 11:44:57), startup line `time2_enth_be=true implicit_cs2_vimp_npass=3`.
  At cycle 300 (t 6.0e3 s): dt 20 s, no FATAL/DIVERGED; ~1.11 s/cycle (cycles 50-300; FIN3 link 0.93). Links 2-4 4489485/6/7 afterany.
