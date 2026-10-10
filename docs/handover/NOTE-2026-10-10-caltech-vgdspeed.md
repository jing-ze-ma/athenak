# NOTE Caltech -> viper: H200 timing of vgdspeed-1009 (TASK-2026-10-10-caltech-vgdspeed)

**Result: NEW is bitwise identical to BASE. The median s/cycle drops 30 % on AG Car B and 31 % on BSG reduced. The
vet_gd halo cost drops 2.4-2.5x.**

## Binaries (build_inc_nofma.sh, incremental from one build dir, same options as 98835d99_nofma: CUDA 12.9,
HOPPER90, MPI, -ffp-contract=off; expansion job 4304113)
- BASE = fork/rt-integration **b2f2e897**: athena_gpu_he_star_m1_b2f2e897_nofma, md5 **0c43b8066c31cbf7064fa2918154fcd5**
- NEW = fork/vgdspeed-1009 tip **ec0757c1** (the revert of 7cffb17a): athena_gpu_he_star_m1_ec0757c1_nofma,
  md5 **a49754374b6e4d75b9f05676aa673841** (40 files recompiled on top of BASE)

## Job
4304128, hpc-sm-01-17, 1 node x 2 H200 (2 ranks, 2 MeshBlocks per GPU), elapsed 8 min 36 s. prof_cuda.sh (unchanged)
via vgd.sbatch, KOKKOS_TOOLS_TIMER_BINARY=1 (real .dat per rank this time, kp_reader -> kt_*.txt), interleaved BASE,
NEW, BASE, NEW. Every arm: rc 0, no FATAL, OOM, NaN or NON-CONVERGED.
- AG Car case = the input prof_cuda.sh names: **agcar_rcxB_ge_accel_st_mg_hr.athinput (AG Car B hr)**, fresh,
  480x128x128 -- NOT the "AG Car A hr" of the Raven numbers.
- BSG = bsg_hr_dc5.athinput on prof_cuda's reduced 256x128x128 mesh, dumps off.

## s/cycle (mean(6-N) / median), peak GPU memory
| arm | BASE rep 1 | NEW rep 1 | BASE rep 2 | NEW rep 2 |
|---|---|---|---|---|
| agc_plain (60) | 0.5585 / 0.4688 | **0.4760 / 0.3292** | 0.5559 / 0.4689 | **0.4385 / 0.3288** |
| agc_kt10 | 0.8892 / 0.5565 | 0.7620 / 0.3813 | - | - |
| agc_kt60 | 0.6306 / 0.5441 | 0.4809 / 0.3624 | - | - |
| bsg_plain (30) | 0.4271 / 0.2989 | **0.3265 / 0.2076** | 0.4004 / 0.2996 | **0.3289 / 0.2070** |
| bsg_kt5 | 0.9042 / 0.4876 | 0.8102 / 0.3671 | - | - |
| bsg_kt30 | 0.4574 / 0.3564 | 0.3804 / 0.2358 | - | - |
| peak/GPU agc | 39705 MiB | 44263-44281 MiB | 39705 | 44263-44275 |
| peak/GPU bsg | 34121 MiB | 38255-38275 MiB | 34121 | 38255-38275 |

- Plain medians: AG Car 0.469 -> 0.329 (-30 %); BSG 0.299 -> 0.207 (-31 %). Repeats agree to 0.1 % on the median.
- kt steady (long minus short, fenced): AG Car 0.5815 -> 0.4172 (-28 %); BSG 0.3892 -> 0.3039 (-22 %).
- Raven expected: AG Car A hr 0.527 -> 0.403 (-24 %), BSG 0.391 -> 0.299 (-24 %). The H200 gain is somewhat larger.
- Memory: NEW is +4.5 GB per GPU, consistent with the 3.6-3.7 GB of lists per rank.

## Grouped (prof_group.py with `hl_` added to the "vet_gd halo pack/unpack" regex; steady cycles, mean of 2 ranks)
| group | AG Car BASE | AG Car NEW | BSG BASE | BSG NEW |
|---|---|---|---|---|
| vet_gd halo pack/unpack | 31.1 % (0.1505 s/cyc) | **19.3 % (0.0635)** | 36.8 % (0.1150) | **22.2 % (0.0459)** |
| vet_gd sweep (m1_vgd_shell, mom) | 35.8 % (0.1733) | 32.1 % (0.1057) | 29.6 % (0.0924) | 27.0 % (0.0558) |
| vet_gd twin | 1.9 % (0.0091) | 2.8 % (0.0091) | 2.8 % (0.0088) | 4.2 % (0.0088) |
| half-range pass | 2.0 % | 2.9 % | 2.4 % | 3.6 % |
| vet_col/lat | 3.0 % | 4.4 % | 1.1 % | 1.6 % |
| rad assembly / mg-line / Krylov / Picard | 3.9 / 3.8 / 0.8 / 2.3 % | 5.7 / 5.6 / 1.2 / 3.4 % | 2.4 / 7.4 / 2.5 / 5.4 % | 3.6 / 11.2 / 3.8 / 8.0 % |
| unmapped (m1_impl_asm, tcell, opacity, f2/f3face, ...) | 6.1 % | 9.0 % | 3.0 % | 4.5 % |
| hydro (fluxes/update/c2p) | 6.3 % | 9.2 % | 2.9 % | 4.3 % |
| kernels, s/cycle | 0.4837 | 0.3294 | 0.3127 | 0.2067 |
| outside kernels, s/cycle (% of wall) | 0.098 (16.8 %) | 0.088 (21.0 %) | 0.076 (19.7 %) | 0.097 (32.0 %) |
| rad/(rad+hydro) | 93.1 % | 89.4 % | 97.0 % | 95.3 % |

- vet_gd halo time per cycle drops 2.4x (AG Car) and 2.5x (BSG). The shell kernel drops 1.6x (shell list).
- Raven expected halo share 42 -> 27 % / 49 -> 31 %; here 31 -> 19 % / 37 -> 22 %.
- NEW halo kernels, AG Car (BSG): m1_vgd_hl_scatter 0.0206 (0.0179) s/cyc, hl_zero 0.0104 (0.0081), hl_pack 0.0079
  (0.0051), hl_scan 0.0052 (0.0031), the remaining hc_scan/unpack/hc_flag about 0.005-0.01. hl_scatter is now the
  biggest halo kernel.
- New label m1_vgd_shl lands in "unmapped" (0.001 s, negligible).

## `<rad_m1> vgdspeed-1009` lines (rank 0) of each NEW log
```
new1/agc_plain: vet_gd halo lists 2.952000e+03 made, 3.686910e+03 MB of 4.096000e+03 (key 4096, GPU: <= free/2); exchanges by list 2.656800e+04, dense 2.952000e+03
new1/agc_kt10: vet_gd halo lists 4.920000e+02 made, 3.684227e+03 MB of 4.096000e+03 (key 4096, GPU: <= free/2); exchanges by list 4.428000e+03, dense 4.920000e+02
new1/agc_kt60: vet_gd halo lists 2.952000e+03 made, 3.686910e+03 MB of 4.096000e+03 (key 4096, GPU: <= free/2); exchanges by list 2.656800e+04, dense 2.952000e+03
new1/bsg_plain: vet_gd halo lists 6.300000e+02 made, 3.586300e+03 MB of 4.096000e+03 (key 4096, GPU: <= free/2); exchanges by list 5.670000e+03, dense 6.300000e+02
new1/bsg_kt5: vet_gd halo lists 2.100000e+02 made, 3.584245e+03 MB of 4.096000e+03 (key 4096, GPU: <= free/2); exchanges by list 8.400000e+02, dense 2.100000e+02
new1/bsg_kt30: vet_gd halo lists 6.300000e+02 made, 3.586300e+03 MB of 4.096000e+03 (key 4096, GPU: <= free/2); exchanges by list 5.670000e+03, dense 6.300000e+02
new2/agc_plain: vet_gd halo lists 2.952000e+03 made, 3.686910e+03 MB of 4.096000e+03 (key 4096, GPU: <= free/2); exchanges by list 2.656800e+04, dense 2.952000e+03
new2/bsg_plain: vet_gd halo lists 6.300000e+02 made, 3.586300e+03 MB of 4.096000e+03 (key 4096, GPU: <= free/2); exchanges by list 5.670000e+03, dense 6.300000e+02
```
About 49 (AG Car) / 21 (BSG) list builds per cycle. 90 % of exchanges go by list; the dense ones are the
first exchange after each build.

## Bitwise (cheap check)
- BASE vs NEW plain arms are **byte-identical**, in both repeats, for: hydro.hst and user.hst (AG Car and BSG); all
  bin dumps (hydro_w, m1, plus m1_face and m1_vet for BSG); and both rst files, including the final-cycle 00001
  outputs (AG Car: 8 files, BSG: 10).
- BASE rep 1 vs rep 2 is identical too.

## Odd
- In BSG, NEW has more host-side time outside kernels: 0.076 -> 0.097 s/cycle, 32 % of wall. That is plausibly the
  list-build bookkeeping (21 builds per cycle). AG Car does not show it (0.098 -> 0.088).
- The 2-GPU H200 layout runs 2 MeshBlocks per GPU (Raven: 1 per A100).

## Files
docs/handover/caltech-vgdspeed-1010/results/4304128/: the four *_prof.txt, the job log j.4304128.out, vgd.sbatch,
prof_group_hl.py. Run dirs are at /resnick/groups/carnegie_poc/jingze/vgdspeed_1010/runs/4304128 (scratch).
