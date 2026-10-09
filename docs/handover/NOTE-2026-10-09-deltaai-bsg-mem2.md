# NOTE DeltaAI: BSG memory round 2 (vet_gd_twin_lowmem) -- RESULTS

**lm_c8 fits: yes (peak 65.0 GiB/GPU = 66567 MiB nvidia-smi; KT 69.2 GB), bitwise lm/new vs base: PASS**

Answer to TASK-2026-10-09-deltaai-bsg-mem2 (bsg-files-1009 9df6f02f). Job **3348593**, 1 node gh054 (4 GH200),
ghx4-interactive, elapsed 1 min 44 s = 0.12 GPU-h. With the failed first try (below) the total is **0.17 GPU-h** (budget 2).

**Headline**
- **Caltech layout (2 ranks x 8 blocks, full mesh) with `vet_gd_twin_lowmem = true` runs:** 6 cycles, rc 0, no FATAL.
  Peak 65.0 GiB per GPU, i.e. **8.65 GB per block** (KT 69.2 GB / 8), under the ~15 GB per block target.
  There is now ONE 40.96 GB band array (the log reports `vet_gd_twin_fuse=0 vet_gd_twin_lowmem=1 vet_gd_twin_det=1`).
- Steady s/cycle at 8 blocks per GPU: **1.37** (16 blocks on 2 GH200).
- **Cost of lowmem:** +8% at 4 blocks per GPU (0.91 vs 0.84-0.85 s/cycle), +35% at 1 block per GPU (0.31 vs
  0.23-0.24).
- **Memory saved:** 4 blocks per GPU 73.6 -> 42.9 GB (nvidia-smi); 1 block 17.7 -> 10.7 GB (KT 10.47 GB).
- **Bitwise: PASS for all three pairs.** All 12 bins + 2 hst; data after `<par_end>` identical 12/12, and bitcmp.py
  BITWISE_PASS. base_c4 vs new_c4 is 14/14 byte-identical (cmp). For base vs lm, `cmp` reports 12 of 14 differing:
  ONLY the bin header, which lists the extra input line `vet_gd_twin_lowmem = true` (+138 bytes); the 2 hst are
  byte-identical.

**Script change (please fix mem_cuda2.sh; 2cad31a2 only adds compares of optional arms, not run here):** the first try (job 3348542) died in 4 s in lm_c8 and lm_c4 with
`### FATAL ERROR in parameter_input.cpp at line 430: Parameter 'vet_gd_twin_lowmem' in block 'rad_m1' on command line
not found`. AthenaK accepts a command-line key only if it already exists in the input file. Per the smoke rule it was
cancelled after 46 s (0.05 GPU-h); base_c4 had already finished fine (73.6 GB, 0.87 s/cycle).
Fix used here: `files/bsg_hr_dc5_lm.athinput` = bsg_hr_dc5.athinput + one line `vet_gd_twin_lowmem = true` after
`vet_gd_twin_fuse` in `<rad_m1>`. The lm arms read that file instead of the command-line key (local copy
`mem_cuda2_local.sh`, a 2-line diff); everything else is unchanged.

## Per-arm lines (mem_cuda2.sh)
```
== lm_c8 eec8bf61c27318a87678b496ce5120f4 n=2 full nlim=6 keys=[(in file)] start 11:23:42
rc=0 fatal=0 oom=0 wall=19s s/cycle(2-4)=1.414 peak: gpu0=66567MiB gpu1=66567MiB gpu2=0MiB gpu3=0MiB  11:24:01
<rad_m1> vet_gd: ragged per-shell band intensity array 4.096431e+01 GB on rank 0 (dense band 7.417627e+01 GB)
== lm_c4 eec8bf61c27318a87678b496ce5120f4 n=4 full nlim=6 keys=[(in file)] start 11:24:01
rc=0 fatal=0 oom=0 wall=17s s/cycle(2-4)=0.942 peak: gpu0=42873MiB gpu1=42873MiB gpu2=42873MiB gpu3=42873MiB  11:24:18
<rad_m1> vet_gd: ragged per-shell band intensity array 2.993509e+01 GB on rank 0 (dense band 3.708813e+01 GB)
== base_c4 c9c6d16705d07dc6e97de6c808d02182 n=4 full nlim=6 keys=[] start 11:24:19
rc=0 fatal=0 oom=0 wall=17s s/cycle(2-4)=0.869 peak: gpu0=73601MiB gpu1=73601MiB gpu2=73601MiB gpu3=73601MiB  11:24:36
<rad_m1> vet_gd: ragged per-shell band intensity array 2.993509e+01 GB per rank (dense band 3.708813e+01 GB)
<rad_m1> vet_gd: ragged per-shell band intensity array 2.993509e+01 GB per rank (dense band 3.708813e+01 GB)
== new_c4 eec8bf61c27318a87678b496ce5120f4 n=4 full nlim=6 keys=[] start 11:24:37
rc=0 fatal=0 oom=0 wall=19s s/cycle(2-4)=0.868 peak: gpu0=72691MiB gpu1=72691MiB gpu2=72691MiB gpu3=72691MiB  11:24:56
<rad_m1> vet_gd: ragged per-shell band intensity array 2.993509e+01 GB on rank 0 (dense band 3.708813e+01 GB)
<rad_m1> vet_gd: ragged per-shell band intensity array 2.993509e+01 GB on rank 0 (dense band 3.708813e+01 GB)
== lm_r4 eec8bf61c27318a87678b496ce5120f4 n=4 red nlim=5 keys=[(in file)] start 11:24:56
rc=0 fatal=0 oom=0 wall=11s s/cycle(2-4)=0.321 peak: gpu0=10681MiB gpu1=10681MiB gpu2=10681MiB gpu3=10681MiB  11:25:07
<rad_m1> vet_gd: ragged per-shell band intensity array 6.800993e+00 GB on rank 0 (dense band 9.272033e+00 GB)
== base_r4 c9c6d16705d07dc6e97de6c808d02182 n=4 red nlim=5 keys=[] start 11:25:07
rc=0 fatal=0 oom=0 wall=9s s/cycle(2-4)=0.245 peak: gpu0=17741MiB gpu1=17741MiB gpu2=17741MiB gpu3=17741MiB  11:25:16
<rad_m1> vet_gd: ragged per-shell band intensity array 6.800993e+00 GB per rank (dense band 9.272033e+00 GB)
<rad_m1> vet_gd: ragged per-shell band intensity array 6.800993e+00 GB per rank (dense band 9.272033e+00 GB)
  DIFF (all 12 bins: header only, see above)
bitwise base_c4 vs lm_c4: 14 files, 12 differ
bitwise base_c4 vs new_c4: 14 files, 0 differ
  DIFF (all 12 bins: header only, see above)
bitwise base_r4 vs lm_r4: 14 files, 12 differ
```
(The lm arms also log `accel-1009 (rank 0): vet_gd_twin_fuse=0 vet_gd_twin_lowmem=1 vet_gd_twin_det=1 ...`.)

## Steady s/cycle (cycles 2-4 from elapsed=)
| arm | binary | layout | key | per-cycle wall (s), cycles 1..n | steady (2-4) | peak nvidia-smi |
|---|---|---|---|---|---|---|
| lm_c8 | BN 6c5d8fb2 | 2 ranks x 8 blocks (full) | lowmem | 2.12 1.51 1.38 1.36 2.36 4.33 | 1.414 | 66567 MiB |
| lm_c4 | BN | 4 x 4 (full) | lowmem | 1.31 0.99 0.92 0.91 1.50 2.35 | 0.942 | 42873 MiB |
| base_c4 | BO 98835d99 | 4 x 4 (full) | -- | 1.26 0.92 0.85 0.84 3.42 2.75 | 0.869 | 73601 MiB |
| new_c4 | BN | 4 x 4 (full) | as input (fused) | 1.20 0.86 0.97 0.77 1.32 2.16 | 0.868 | 72691 MiB |
| lm_r4 | BN | 4 x 1 (reduced) | lowmem | 0.53 0.34 0.32 0.31 1.05 | 0.321 | 10681 MiB |
| base_r4 | BO | 4 x 1 (reduced) | -- | 0.46 0.26 0.24 0.23 0.92 | 0.245 | 17741 MiB |

The script's "steady" averages cycles 2-4, which still includes the slower cycle 2. Cycles 3-4 alone: lm_c8 1.37,
lm_c4 0.915, base_c4 0.845, lm_r4 0.315, base_r4 0.235.

## Binaries
- BO `athena_hes_gpu_98835d99a16f` md5 `c9c6d16705d07dc6e97de6c808d02182`
- BN `athena_hes_gpu_6c5d8fb2c169` md5 `eec8bf61c27318a87678b496ce5120f4` (mem-1009 6c5d8fb2, build_inc_deltaai.sh hes_gpu)
- Binding as round 1 (`srun -n 2|4 -c 72 --cpu-bind=cores` + `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID`). KT on.

## ana_mem.py, top 15 (rank 0, ncells_per_block 1283800)

### lm_c8
## PEAK at t=14.03 s: total 69.202 GB (1082 allocations)
| label | space | n | GB per rank | GB per block | doubles/cell |
|---|---|---|---|---|---|
|  | Cuda | 940 | 62.821 | 7.853 | 764.6 |
| derived-var | Cuda | 7 | 0.657 | 0.082 | 8.0 |
| cons | Cuda | 1 | 0.411 | 0.051 | 5.0 |
| prim | Cuda | 1 | 0.411 | 0.051 | 5.0 |
| cons1 | Cuda | 1 | 0.411 | 0.051 | 5.0 |
| uflx.x1f | Cuda | 1 | 0.411 | 0.051 | 5.0 |
| uflx.x2f | Cuda | 1 | 0.411 | 0.051 | 5.0 |
| uflx.x3f | Cuda | 1 | 0.411 | 0.051 | 5.0 |
| utest | Cuda | 1 | 0.411 | 0.051 | 5.0 |
| rst-m1s | Cuda | 1 | 0.411 | 0.051 | 5.0 |
| m1_u0 | Cuda | 1 | 0.329 | 0.041 | 4.0 |
| dxe.x1e | Cuda | 1 | 0.085 | 0.011 | 1.0 |

### lm_r4
## PEAK at t=3.77 s: total 10.470 GB (1082 allocations)
| label | space | n | GB per rank | GB per block | doubles/cell |
|---|---|---|---|---|---|
|  | Cuda | 940 | 9.634 | 9.634 | 938.0 |
| derived-var | Cuda | 7 | 0.082 | 0.082 | 8.0 |
| cons | Cuda | 1 | 0.051 | 0.051 | 5.0 |
| prim | Cuda | 1 | 0.051 | 0.051 | 5.0 |
| cons1 | Cuda | 1 | 0.051 | 0.051 | 5.0 |
| uflx.x1f | Cuda | 1 | 0.051 | 0.051 | 5.0 |
| uflx.x2f | Cuda | 1 | 0.051 | 0.051 | 5.0 |
| uflx.x3f | Cuda | 1 | 0.051 | 0.051 | 5.0 |
| utest | Cuda | 1 | 0.051 | 0.051 | 5.0 |
| rst-m1s | Cuda | 1 | 0.051 | 0.051 | 5.0 |
| m1_u0 | Cuda | 1 | 0.041 | 0.041 | 4.0 |
| m1_hm_sbuf | Cuda | 1 | 0.023 | 0.023 | 2.2 |

lm_c8 largest device allocations (rank 0): 40.964 GB (the single band array), 4.355 GB, 1.561 GB.

Run dirs kept: `/work/nvme/bivj/jma20/bsg_hrdet_1009/runs/mem2_3348593/` (all 6 arms with .mem_events; mem2_3348542_keyfail
= the failed first try). Nothing else queued on DeltaAI.
