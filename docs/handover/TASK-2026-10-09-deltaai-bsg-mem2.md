# TASK for DeltaAI: BSG memory round 2 -- vet_gd_twin_lowmem at 8 blocks per GPU (mem-1009 6c5d8fb2, user 10-09)

Follow-up of NOTE-2026-10-09-deltaai-bsg-mem (thanks: it found the 2nd band copy). The 2nd copy is the fused
twin's intensity array (vgd_itw). New opt-in key `rad_m1/vet_gd_twin_lowmem = true` (with vet_gd_twin_fuse) sweeps
the twin first in the main array instead: expected bitwise identical, ~half the band memory, some extra s/cycle.
**Budget: 1 node (4 GH200), <= 30 min, <= 2 GPU-h.** Same dir/account/rules as the round-1 task.

## Build
BN = `mem-1009` @ **6c5d8fb2c** (`git fetch fork mem-1009`; `build_inc_deltaai.sh hes_gpu`). BO = the 98835d99
binary of round 1 (md5 c9c6d167...). Record md5s.

## Job (1 node, same SRUN2/SRUN4 binding as round 1, `--time=00:30:00`)
```
git fetch fork bsg-files-1009; git archive fork/bsg-files-1009 docs/handover/mem-1009 | tar -x -C $G/bundle
SRUN2=... SRUN4=... bash $G/bundle/docs/handover/mem-1009/mem_cuda2.sh <BO> <BN> $G/files $G/runs/mem2_$SLURM_JOB_ID $KT
```
Arms (in order): lm_c8 (FULL mesh, 2 ranks x 8 blocks, lowmem = the Caltech layout: does it fit in 96 GB, s/cycle),
lm_c4 / base_c4 / new_c4 (full mesh 4 x 4: bitwise lowmem vs 98835d99, keys-as-input vs 98835d99, s/cycle cost),
lm_r4 / base_r4 (1 block per GPU, bitwise + peak per block). OOM is a result; the script continues.
Smoke rule: if lm_c8 dies with anything but OOM in the first minute, cancel and report.

## NOTE (NOTE-2026-10-09-deltaai-bsg-mem2.md on this branch)
Verdict line first: `lm_c8 fits: yes|no (peak X GiB/GPU), bitwise lm/new vs base: PASS|FAIL`. Then the per-arm lines,
the three `bitwise ... vs ...` lines, steady s/cycle table (cycles 2-4 from elapsed=), md5s, job id, GPU-h, and for
lm_c8 and lm_r4 the top 15 lines of `ana_mem.py <rank-0 .mem_events> 1283800 <blocks per rank>`. Keep run dirs.
