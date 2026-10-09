# NOTE DeltaAI: BSG GPU memory per block (mem-1009) -- RESULTS

Answer to TASK-2026-10-09-deltaai-bsg-mem (bsg-files-1009 2049149a). Job **3347781**, 1 node gh151 (4 GH200, 96 GB
HBM each), ghx4-interactive, elapsed 1 min 23 s = 4 GPUs x 83 s = **0.09 GPU-h** (budget 2). KT memory-events ON.

**Headline**
- **mem_c8 (Caltech layout, 8 blocks per GPU) = OOM on GH200.** On both ranks, `Kokkos ERROR: Cuda memory space
  failed to allocate 38.15 GiB (label="")` right after the log line `vet_gd: ragged per-shell band intensity array
  4.096431e+01 GB on rank 0 (dense band 7.417627e+01 GB)`. The run allocates that array TWICE (two log lines in every
  other arm). The first 40.96 GB copy fit (nvidia-smi 63.8 GB at the crash); the second one did not.
- **~77-80% of the GPU memory is these 2 unlabelled band arrays** (inside the `label=""` row of ana_mem.py, which is
  95% of the total):
  2 x 6.80 GB per block at 1 block per rank, 2 x 29.94 GB = 15.0 GB per block at 4 blocks per rank,
  2 x 40.96 GB = 10.2 GB per block at 8 blocks per rank (the ragged band shrinks as more neighbours are on-rank).
  Everything else is ~3.9-4.0 GB per block.
- **Projection for 8 blocks per GPU: 2 x 40.96 + 8 x 3.9 = ~113 GB (~14.1 GB per block)** -> under the ~125 GB per
  H200 target, with ~12 GB margin (plus the CUDA context, ~0.6 GB). This is an extrapolation, NOT a measurement:
  nothing above 96 GB can be measured on GH200.
- **mem-1009 saves only ~0.24 GB per block** (KT peak 17.86 -> 17.63 GB at 1 block per rank, 76.45 -> 75.50 GB at 4):
  0.164 GB of it is m1_u1 + m1_uflx.x1f/x2f/x3f (4 x 0.041 GB, gone from the table), and 0.07 GB comes out of the
  unlabelled row (vet_col_lat sweep arrays / closure slabs, presumably). The band arrays are the SAME size in base and mem
  (6.800993 / 29.93509 GB), so the "band only on off-rank sides" cherry-pick does not reduce them in these layouts.
- **The lever is the band array pair:** halving it (one copy instead of two; I did not check what the 2nd copy is,
  possibly the twin's) would bring
  8 blocks per GPU to ~72 GB, which fits even on GH200.
- bitwise base_r4 vs mem_r4: **14 files, 0 differ** (cmp, headers included; payload check also 12/12 identical).

## Per-arm lines (mem_cuda.sh)
```
== mem_c8 84888802cfec21a127b527c686f1629e full nlim=6 start 10:50:20
rc=134 fatal=0 oom=2 wall=10s s/cycle=na peak: gpu0=63773MiB gpu1=63773MiB gpu2=1MiB gpu3=1MiB  10:50:30
== mem_r4 84888802cfec21a127b527c686f1629e red nlim=5 start 10:50:30
rc=0 fatal=0 oom=0 wall=13s s/cycle=1.00 peak: gpu0=17511MiB gpu1=17511MiB gpu2=17511MiB gpu3=17511MiB  10:50:43
== base_r4 c9c6d16705d07dc6e97de6c808d02182 red nlim=5 start 10:50:44
rc=0 fatal=0 oom=0 wall=13s s/cycle=1.00 peak: gpu0=17741MiB gpu1=17741MiB gpu2=17741MiB gpu3=17741MiB  10:50:57
== mem_c4 84888802cfec21a127b527c686f1629e full nlim=6 start 10:50:57
rc=0 fatal=0 oom=0 wall=20s s/cycle=1.80 peak: gpu0=72691MiB gpu1=72691MiB gpu2=72691MiB gpu3=72691MiB  10:51:17
== base_c4 c9c6d16705d07dc6e97de6c808d02182 full nlim=6 start 10:51:17
rc=0 fatal=0 oom=0 wall=16s s/cycle=1.60 peak: gpu0=73611MiB gpu1=73611MiB gpu2=73611MiB gpu3=73611MiB  10:51:33
bitwise base_r4 vs mem_r4: 14 files, 0 differ
```
The script's s/cycle uses 1-s timestamps and includes the dump cycles. Per-cycle from `elapsed=` (cycles 2-4 = steady):

| arm | binary | blocks/GPU | per-cycle wall (s) | steady s/cycle |
|---|---|---|---|---|
| mem_r4 | BN af5147cc | 1 (reduced mesh) | 0.47 0.27 0.25 0.24 3.13 | 0.24 |
| base_r4 | BO 98835d99 | 1 (reduced mesh) | 0.46 0.26 0.24 0.23 2.93 | 0.24 |
| mem_c4 | BN | 4 (full mesh) | 1.21 0.86 0.79 0.78 2.57 3.90 | 0.79 |
| base_c4 | BO | 4 (full mesh) | 1.27 0.92 0.86 0.85 1.51 2.91 | 0.86 |
| mem_c8 | BN | 8 (full mesh) | OOM in setup (10 s) | -- |

Same 16-block full mesh: 0.79 s/cycle on 4 GPUs (4 ranks x 4 blocks) here, versus 0.86 s/cycle on 8 GPUs
(16 ranks, 2 ranks per GPU) in the hrdet gate. Fewer ranks with more blocks each are at least as fast here.

## Binaries, setup
- BO `athena_hes_gpu_98835d99a16f` md5 `c9c6d16705d07dc6e97de6c808d02182` (reused from the hrdet gate)
- BN `athena_hes_gpu_af5147cc208b` md5 `84888802cfec21a127b527c686f1629e` (mem-1009, build_inc_deltaai.sh hes_gpu)
- KT: kokkos-tools master 1e98316, `profiling/memory-events`, `make CXX=g++`.
- Binding: SRUN2/SRUN4 = `srun -n 2|4 -c 72 --cpu-bind=cores` + wrapper `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID`
  (1 rank per Grace socket and its GH200). mem_c8 used GPUs 0 and 1.
- nvidia-smi peaks: mem_c8 63773 MiB (at the crash), mem_r4 17511, base_r4 17741, mem_c4 72691, base_c4 73611 MiB.
- mem_c8 has no `.mem_events` file: KT writes it at finalize, and the run aborted.

## ana_mem.py tables (rank 0; ncells_per_block 1283800; top 40 lines)

### mem_r4 (blocks per rank 1; gh151.hsn.cm.delta.internal.ncsa.edu-2432177.mem_events)
## PEAK at t=7.42 s: total 17.626 GB (1091 allocations)
| label | space | n | GB per rank | GB per block | doubles/cell |
|---|---|---|---|---|---|
|  | Cuda | 949 | 16.790 | 16.790 | 1634.8 |
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
| m1_hm_rbuf | Cuda | 1 | 0.023 | 0.023 | 2.2 |


### base_r4 (blocks per rank 1; gh151.hsn.cm.delta.internal.ncsa.edu-2432325.mem_events)
## PEAK at t=6.21 s: total 17.862 GB (1090 allocations)
| label | space | n | GB per rank | GB per block | doubles/cell |
|---|---|---|---|---|---|
|  | Cuda | 948 | 16.862 | 16.862 | 1641.8 |
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
| m1_u1 | Cuda | 1 | 0.041 | 0.041 | 4.0 |
| m1_uflx.x1f | Cuda | 1 | 0.041 | 0.041 | 4.0 |
| m1_uflx.x2f | Cuda | 1 | 0.041 | 0.041 | 4.0 |
| m1_uflx.x3f | Cuda | 1 | 0.041 | 0.041 | 4.0 |
| m1_hm_sbuf | Cuda | 1 | 0.023 | 0.023 | 2.2 |
| m1_hm_rbuf | Cuda | 1 | 0.023 | 0.023 | 2.2 |


### mem_c4 (blocks per rank 4; gh151.hsn.cm.delta.internal.ncsa.edu-2432472.mem_events)
## PEAK at t=12.83 s: total 75.504 GB (1091 allocations)
| label | space | n | GB per rank | GB per block | doubles/cell |
|---|---|---|---|---|---|
|  | Cuda | 949 | 72.262 | 18.065 | 1759.0 |
| derived-var | Cuda | 7 | 0.329 | 0.082 | 8.0 |
| cons | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| prim | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| cons1 | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| uflx.x1f | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| uflx.x2f | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| uflx.x3f | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| utest | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| rst-m1s | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| m1_u0 | Cuda | 1 | 0.164 | 0.041 | 4.0 |
| m1_hm_sbuf | Cuda | 1 | 0.047 | 0.012 | 1.1 |
| m1_hm_rbuf | Cuda | 1 | 0.047 | 0.012 | 1.1 |
| dxe.x1e | Cuda | 1 | 0.042 | 0.011 | 1.0 |
| areae.x1e | Cuda | 1 | 0.042 | 0.011 | 1.0 |
| dxe.x2e | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| dxe.x3e | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| areae.x2e | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| areae.x3e | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| area.x2f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| area.x3f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| dxf.x2f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| dxf.x3f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| phi_fc.x2f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| phi_fc.x3f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| area.x1f | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| dxf.x1f | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| phi_fc.x1f | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| volume | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| dx1 | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| dx2 | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| dx3 | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| x_ov_rD | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| y_ov_rC | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| z_ov_rC | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| phi_cc | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| fofc_cnt | Cuda | 1 | 0.041 | 0.010 | 1.0 |

### base_c4 (blocks per rank 4; gh151.hsn.cm.delta.internal.ncsa.edu-2433634.mem_events)
## PEAK at t=9.74 s: total 76.449 GB (1090 allocations)
| label | space | n | GB per rank | GB per block | doubles/cell |
|---|---|---|---|---|---|
|  | Cuda | 948 | 72.549 | 18.137 | 1766.0 |
| derived-var | Cuda | 7 | 0.329 | 0.082 | 8.0 |
| cons | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| prim | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| cons1 | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| uflx.x1f | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| uflx.x2f | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| uflx.x3f | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| utest | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| rst-m1s | Cuda | 1 | 0.205 | 0.051 | 5.0 |
| m1_u0 | Cuda | 1 | 0.164 | 0.041 | 4.0 |
| m1_u1 | Cuda | 1 | 0.164 | 0.041 | 4.0 |
| m1_uflx.x1f | Cuda | 1 | 0.164 | 0.041 | 4.0 |
| m1_uflx.x2f | Cuda | 1 | 0.164 | 0.041 | 4.0 |
| m1_uflx.x3f | Cuda | 1 | 0.164 | 0.041 | 4.0 |
| m1_hm_sbuf | Cuda | 1 | 0.047 | 0.012 | 1.1 |
| m1_hm_rbuf | Cuda | 1 | 0.047 | 0.012 | 1.1 |
| dxe.x1e | Cuda | 1 | 0.042 | 0.011 | 1.0 |
| areae.x1e | Cuda | 1 | 0.042 | 0.011 | 1.0 |
| dxe.x2e | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| dxe.x3e | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| areae.x2e | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| areae.x3e | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| area.x2f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| area.x3f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| dxf.x2f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| dxf.x3f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| phi_fc.x2f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| phi_fc.x3f | Cuda | 1 | 0.042 | 0.010 | 1.0 |
| area.x1f | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| dxf.x1f | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| phi_fc.x1f | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| volume | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| dx1 | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| dx2 | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| dx3 | Cuda | 1 | 0.041 | 0.010 | 1.0 |
| x_ov_rD | Cuda | 1 | 0.041 | 0.010 | 1.0 |

Run dirs kept: `/work/nvme/bivj/jma20/bsg_hrdet_1009/runs/mem_3347781/{mem_c8,mem_r4,base_r4,mem_c4,base_c4}`
(4 `.mem_events` files per arm except mem_c8, `smi.log`, `out.log`). Nothing else queued on DeltaAI.
