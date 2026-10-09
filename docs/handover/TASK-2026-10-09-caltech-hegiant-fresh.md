# TASK for Caltech (10-09): queue a copy of the He giant FRESH N897 run (first start wins vs Raven)

**User GO 10-09.** From: viper. To: Caltech (Resnick `gpu` partition, H200 141 GB, 4 GPUs per node).
Everything is on fork branch **`bsg-files-1009`**: this file + `docs/handover/hegiant-fresh-1009/`.
**Do not touch the BSG chain 4286643-45** (it keeps its 2 H200). This run is queued now; it may start when BSG frees
GPUs or by backfill (short 6 h links help backfill).

## What it is
t = 0 start of the He giant directly on the N897 grid (fine band 60-95 Rsun, 4 MeshBlocks of 897 x 64 x 64) with the
10-09 production radiation (half-range VET blend, vet_gd_twin + noesrc + fuse/det, vet_scatter ma / kappa_e eos, F2/F3,
precond mg), tlim 2.592e6 s (30 d), cfl 0.3. It replaces the remapped N897 continuation
(`/resnick/groups/carnegie_poc/jingze/hegiant_1007/n897`, stopped at 23.50 d; keep that dir, do not touch it).
Setup and viper/Raven smokes: viper `/viper/ptmp2/jinma/he_giant_1006/FRESH_1009.md` (summary in the table below).
The same run is queued on Raven: chain 31032793 -> 94 -> 95 -> 96 -> 97 (12 h links, 4 A100 x 1 block). **First start wins.**

## 1. Binary (no new build)
Use your existing `athena_gpu_he_star_m1_6c5d8fb2_nofma` (md5 e9873137c319b5f29d85680b535a450c; mem-1009 6c5d8fb2 =
98835d99 + `vet_gd_twin_lowmem`, default off; bitwise to 98835d99 with lowmem off). The input does not set lowmem: keep it off.

## 2. Files and input
```
git fetch <fork remote> bsg-files-1009
git archive FETCH_HEAD docs/handover/hegiant-fresh-1009 | tar -x -C <work>
bash <work>/docs/handover/hegiant-fresh-1009/SETUP.sh <abs files dir>   # gunzip, md5 check, writes the input
```
| file | md5 (uncompressed) |
|---|---|
| `ic_giant_own.txt` (`problem/he_ic_file`) | ddc72692d8999ea33730ae97898f4c96 |
| `rosseland_tops_hegiant_blend.txt` (`problem/he_opac_table`) | c1ef063c50fabbda2ecc029d2dd4e695 |
| `planck_tops_hegiant_blend.txt` (`problem/he_planck_table`) | b699c60db8cc55f6da25c61cd62d6926 |
| `hegiant_fresh897.athinput.in` (template) | 0a01113d5cd5515d146c572b9df8181a |

- The IC and tables are the same files (same md5) as `hegiant-files-1007/` that your N445 run used; reuse those if you
  prefer, the bundle is only so this task is self-contained.
- `hegiant_fresh897.athinput` = viper's input (md5 b30d21be0a2cf897d692b37285021939) with only the 3 file paths ->
  `@FILES@` (checked: substituting the viper paths back gives b30d21be). Every key is unchanged. The grid is analytic
  (18 `f_stretch_r_*` coefficients in the input); no other file is read.
- Command-line keys: `XKEYS = time/cfl_number=0.3` exactly as viper/Raven run.cfg (the same value is in the input).
  Nothing else on the command line except restart `last_time` keys and `-d`/`-t`.
- **Before link 1 starts**, check the fork for `docs/handover/NOTE-2026-10-09-viper-caltech-hegiant-INPUT.md`.
  viper is running an inner-boundary check; if it changes the input, that NOTE gives the revised `.athinput.in` + md5.
  Use the revised input (and re-smoke) if present.

## 3. Layout
1 node, **2 MPI ranks x 2 MeshBlocks, 1 rank per H200** (`--ntasks-per-node=2 --gres=gpu:h200:2`, your
CUDA_VISIBLE_DEVICES=$SLURM_LOCALID wrapper, `--cpu-bind=cores`, as in `bsg_1009/link_bsg.sh`). Expected GPU memory
**~72 GiB/GPU** (Raven measured 36645 MiB per A100 for 1 block/GPU, same input class). Fits 141 GB.

## 4. Smoke (separate dir, same binary + input + XKEYS, plus `time/nlim=10 time/ndiag=1`)
Report: job id, rc, FATAL / NaN / NON-CONVERGED counts in the log, s/cycle (cycles 3-10), GPU memory peak (nvidia-smi).
Production only if clean. viper/Raven reference (25 cycles, Raven 4 x A100, 1 block/GPU, precond mg, jobs 31031521 /
31032400): rc 0, FATAL/DIVERGED/NaN 0, Picard mean/max 9.52/12, NON-CONV 0, **1.086 s/cycle**, dt 20.16 s constant.

## 5. Production
- Run dir e.g. `/resnick/groups/carnegie_poc/jingze/hegiant_1009/fresh897`.
- Fresh start (t = 0) from the input; **6 h chained links** on `gpu` (`--dependency=afterany`), queue ~5-6 links now
  (30 d at dt ~20 s is ~1.3e5 cycles). Each later link restarts from the newest `rst/hegiant.*.rst` with the
  `outputN/last_time` keys printed by `rst_info.py` (line 2), plus XKEYS. rst every 0.5 d (input).
- Template: `hegiant-fresh-1009/link_hegiant.sh.txt` (FRESH if no rst, else newest rst + last_time keys; DONE at tlim;
  STOP on rc != 0 / FATAL / NaN / NON-CONVERGED; md5 check of BIN; CANCEL check). Adapt account / modules / wrapper.
- **Stop rule:** rc != 0, FATAL, NaN or NON-CONVERGED -> STOP file, no further links. DONE when newest rst t >= 2.592e6.

## 6. First start wins + notes back (push to `bsg-files-1009`)
- **The moment link 1 STARTS (RUNNING)**: push `docs/handover/NOTE-2026-10-09-caltech-hegiant-fresh.md` with start time,
  chain job ids, smoke job id + s/cycle + GPU peak, binary md5, input md5. viper then cancels the Raven chain.
- **Before each new link** (in the link script, see template): fetch the fork and stop if
  `docs/handover/NOTE-2026-10-09-viper-caltech-hegiant-CANCEL.md` exists (Raven started first). Then `scancel
  --state=PENDING` the remaining links (fresh `squeue -u $USER` in the same command), keep everything written, and
  append the stop to the NOTE.
- **After each link**, append to the same NOTE: t, cycle, dt, s/cycle, FATAL/NaN/NON-CONV counts, newest rst.
