# NOTE DeltaAI -> Raven: bsgli-val GATE SKIPPED (user); backup production STARTED on 2 h ghx4-interactive links

USER (DeltaAI side, 10-13 ~22:50 CDT): "submit 2 h links", then "skip the gate, start production now".
This overrides TASK-2026-10-14-deltaai-bsgli-val item 4 and the update NOTE's conditions (gate BITWISE + batch < 15 h, no interactive links).

## Build
- fin4-linhao-bc @ **8e0182723da23d8c864c8b079bd24e75689f2d28** (per the redirect/update NOTEs, not eb1e032c)
- `athena_hes_gpu_8e0182723da2` md5 **781f4eb2470ede5dd91a72756054f081** (build_inc_deltaai.sh hes_gpu, same flags as 9ceb9cdc)

## Inputs
- bsgli-val-1014 files: MD5SUMS all OK; linhao_profile_full.txt unpacked md5 41edc225; ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805.
- Deck = bsgLI_bc_deltaai.athinput.in with only @FILES@ resolved: md5 **e7c35b832797b8c013f725c6ff150cc8** (he_bc_linhao = true, halo_ipc + twin_lowmem on).

## Gate
- NOT RUN (user). No C60 / L30+30 / Coff comparison, no separate s/cycle + memory measurement. The first link logs s/cycle and
  GPU memory (5 s sampler); I will report them in a follow-up NOTE.

## Queue estimates (10-13 22:47 CDT, sbatch --test-only, 1 node x 4 GH200)
- ghx4 BATCH 24 h: start 2026-10-19 08:09 CDT (~5.4 d out). ghx4 batch 45 min / 2 h: same.
- ghx4-interactive 45 min / 2 h: start ~15 min out (MaxTime 2 h, 1 submitted job per user).

## Backup production (running)
- Dir /work/nvme/bivj/jma20/bsgliprod_1014; same binary/deck/layout (4 ranks x 4 MeshBlocks, 1 rank/GPU); fresh start, TLIM 4.96e6
  (deck outputs as is: hst 1000 s, bins 0.5 d, rst 1 d).
- 2 h ghx4-interactive links, one at a time (chain.sh); each link restarts from rst_pick.py --quarantine's newest valid rst with its
  last_time keys, runs to wall - 5 min, then prunes rst keep 3.
- STOP on rc!=0 / FATAL / NaN / DIVERGED / NON-CONVERGED > 0 / no valid rst; MAXLINKS 12 (~13-16 h compute expected -> ~8-9 links).
- Link 1: job **3377323** (submitted 10-13 22:55 CDT, pending).
- To stop it from Raven: push any file under docs/handover with "bsgli" and one of cancel|withdrawn|otherfirst|stop in its name
  (checked before every link submit), or say so in a NOTE.
