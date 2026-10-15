# xthinfix-1010 Raven GPU point: AG Car A hr, new arms on e50ec7d6 (1 node x 4 A100)

STATUS: complete. Build 31050046; smoke 31050205; arms 31050237 / 31050238 / 31050239; repeats 31050588
(all on ravg1002, gpudev). Raw: `gpu/RESULTS_gpu_point.txt`; runs `gpu/runs/j.<job>/<arm>` (links in `gpu/runs/all_arms`).

## Verdict

1. Gate PASSES: base (rt-integration 89693b59_vgdmerge, production input unchanged) == all (e50ec7d6, mode all,
   key-complete input). Bin data (after `<par_end>`) and hst data rows are bitwise identical. Repeats are bitwise too:
   all_1 == all_2, sq == sq_2, sqsb == sqsb_2.
2. All 11 runs: rc 0, FATAL 0, NON-CONVERGED 0, 30 cycles. plm positivity fallbacks are 0, except 8 in sqs and 6 in aqs.
3. Cost (s/cycle, cycles 5-30; all = 0.403):

   | arm | s/cycle | x all | Picard mean |
   | --- | --- | --- | --- |
   | sq | 1.85 | 4.6 | 46.6 |
   | aq | 1.76 | 4.4 | 45.9 |
   | sqs (step) | 1.57 | 3.9 | 41.3 |
   | aqs (step) | 1.56 | 3.9 | 42.7 |
   | sqb (bound) | 2.00 | 5.0 | 49.0 |
   | sqp (pred) | 2.60 | 6.4 | 61.2 |
   | **sqsb (step + bound)** | **1.05** | **2.6** | **28.7** |
   | sqsbp (step + bound + pred) | 2.30 | 5.7 | 51.7 |

   - step alone saves about 15 %.
   - bound alone costs slightly more than sq.
   - step + bound is the only clear win: -43 % vs sq. It is still 2.6x all.
   - pred makes every arm MORE expensive. In sqp the dc phase takes 198 of the 1836 Picard passes and the plm phase 1638.
4. Accuracy vs all:
   - Inside 0.8 R_ph, every arm is at roundoff level (<= 1e-5).
   - From 0.8 to 1.05 R_ph, every plm arm differs by the same amount as sq (Caltech-style metric: 0.8-0.95 <= 5e-3, 0.95-1.05 velx ~2-4, F2/F3 0.2-0.5).
   - Above 1.05 R_ph, `bound` is what matters. Normalised by max|all|:
     - sq / sqs / sqp / aq / aqs: tangential velocities 1e2-9e2x all's, F2/F3 6-70x.
     - sqb / sqsb / sqsbp: velocities 4-13x, F2/F3 0.5-1.7x.
     - E and eint differences above 1.05 R_ph also drop, from 4e-2..6e-1 to about 2.5e-2.
   - So bound removes most of the outer-wind tangential blow-up of plm. step leaves accuracy essentially unchanged versus its parent arm.
5. Cross-machine check: sq and aq on A100 reproduce Caltech's H200 413b38af point.
   - Picard mean/max and inner BiCGStab mean are identical: 46.567/53/2.947 and 45.9/51/2.946.
   - The sq-vs-all spdiff bands match to 2 digits in nearly every entry.
6. Net: sqsb is the best new arm on the GPU, at 2.6x all with a much tamer outer region. Nothing tested here reaches the cost of all (Picard 4.3).
   Whether sqsb keeps the convergence order and survives ba0 / ba20 is for the CPU battery.

## Per-arm table (cycles 5-30; peak = max nvidia-smi memory.used over the 4 GPUs, 0.5 s sampling)

| arm | rc | FATAL | s/cycle | Picard mean | max | NON-CONV | inner mean | plm pos-fallbacks | bound: cell-axes -> dc | pred dc/plm passes | peak MiB |
|---|---|---|---|---|---|---|---|---|---|---|---|
| base_1 | 0 | 0 | 0.400 | 4.333 | 9 | 0 | 3.331 | - | - | - | 23659 |
| all_1 | 0 | 0 | 0.403 | 4.333 | 9 | 0 | 3.331 | - | - | - | 23643 |
| all_2 | 0 | 0 | 0.400 | 4.333 | 9 | 0 | 3.331 | - | - | - | 23643 |
| sq | 0 | 0 | 1.850 | 46.567 | 53 | 0 | 2.947 | 0 | - | - | 24025 |
| sq_2 | 0 | 0 | 1.839 | 46.567 | 53 | 0 | 2.947 | 0 | - | - | 24025 |
| sqs | 0 | 0 | 1.569 | 41.267 | 52 | 0 | 3.139 | 8 | - | - | 24025 |
| sqb | 0 | 0 | 2.001 | 48.967 | 66 | 0 | 2.834 | 0 | 7.97e7 | - | 24025 |
| sqp | 0 | 0 | 2.599 | 61.200 | 70 | 0 | 3.288 | 0 | - | 198/1638 | 24025 |
| sqsb | 0 | 0 | 1.046 | 28.700 | 48 | 0 | 2.627 | 0 | 8.73e7 | - | 24041 |
| sqsb_2 | 0 | 0 | 1.047 | 28.700 | 48 | 0 | 2.627 | 0 | 8.73e7 | - | 24025 |
| sqsbp | 0 | 0 | 2.295 | 51.733 | 68 | 0 | 2.765 | 0 | 1.05e8 | 195/1357 | 24025 |
| aq | 0 | 0 | 1.763 | 45.900 | 51 | 0 | 2.946 | 0 | - | - | 24025 |
| aqs | 0 | 0 | 1.563 | 42.733 | 53 | 0 | 3.183 | 6 | - | - | 24025 |

- Peak memory is about 23.6 GB per A100 for dc and about 24.0 GB for plm (+0.4 GB). The run.log also reports a 6.68 GB ragged vet_gd band array on rank 0.

## Differences vs all_1 (last dump, t = 3.7577e4 s, cycle 30)

Method:
- Take the last bin dump (cycle 30) of each output. For each radial shell, compute max over (theta, phi) of |arm - all|.
- Divide by max|arm| in that shell (the Caltech spdiff.py convention, which saturates near 1) or by max|all| (shown after `/`).
- Report the maximum over each r/R_ph band, with R_ph = 388.3 Rsun.

Bands 0-0.8 and 0.8-0.95 have the same values under both normalisations, so they are given once.

| arm | var | <0.8 | 0.8-0.95 | 0.95-1.05 | >1.05 (/arm / /all) |
|---|---|---|---|---|---|
| sq | dens | 0 | 1.3e-6 | 1.2e-2 | 6.0e-3 / 6.0e-3 |
| | velx | 9.8e-8 | 5.4e-3 | 3.9 | 1.0e-1 / 1.1e-1 |
| | vely / velz | 4.1e-6 | 2.8e-4 | 3.6e-1 | 1.1 / **1.3e2** |
| | eint / m1_e | 0 | 1.7e-4 | 8.2e-2 | 9.5e-2 / 1.0e-1 |
| | m1_f1 | 7.9e-8 | 8.5e-5 | 7.1e-4 | 5.7e-2 / 6.0e-2 |
| | m1_f2 / f3 | 6.0e-6 | 2.8e-4 | 2.1e-1 | 1.1 / **9.5** |
| sqs | dens | 0 | 1.2e-6 | 1.1e-2 | 3.2e-3 / 3.2e-3 |
| | velx | 9.9e-8 | 5.1e-3 | 3.0 | 1.6e-1 / 1.9e-1 |
| | vely / velz | 3.6e-6 | 2.8e-4 | 3.7e-1 | 1.1 / **8.9e2** |
| | eint / m1_e | 0 | 1.6e-4 | 7.4e-2 | 1.0e-1 / 1.1e-1 |
| | m1_f1 | 1.1e-7 | 8.1e-5 | 7.8e-4 | 2.4e-2 / 2.4e-2 |
| | m1_f2 / f3 | 5.3e-6 | 2.9e-4 | 4.4e-1 | 1.3 / **8.0** |
| sqb | dens | 0 | 1.3e-6 | 1.2e-2 | 1.5e-3 / 1.5e-3 |
| | velx | 9.8e-8 | 5.4e-3 | 3.9 | 2.1e-2 / 2.0e-2 |
| | vely / velz | 4.3e-6 | 2.8e-4 | 3.2e-1 | 1.0 / 1.3e1 |
| | eint / m1_e | 0 | 1.7e-4 | 8.2e-2 | 5.5e-2 / 5.8e-2 |
| | m1_f1 | 8.0e-8 | 8.5e-5 | 7.1e-4 | 4.0e-3 / 4.0e-3 |
| | m1_f2 / f3 | 6.5e-6 | 2.8e-4 | 2.1e-1 | 4.3e-1 / 5.0e-1 |
| sqp | dens | 0 | 1.3e-6 | 1.2e-2 | 6.0e-3 / 6.0e-3 |
| | velx | 9.8e-8 | 5.4e-3 | 3.9 | 3.3e-1 / 5.0e-1 |
| | vely / velz | 3.7e-6 | 2.9e-4 | 3.7e-1 | 1.0 / **8.7e2** |
| | eint / m1_e | 0 | 1.7e-4 | 8.2e-2 | 3.7e-1 / 5.8e-1 |
| | m1_f1 | 8.0e-8 | 8.5e-5 | 6.6e-4 | 2.9e-1 / 4.0e-1 |
| | m1_f2 / f3 | 5.7e-6 | 2.8e-4 | 2.0e-1 | 1.2 / **6.9e1** |
| sqsb | dens | 0 | 1.2e-6 | 1.1e-2 | 1.5e-3 / 1.5e-3 |
| | velx | 9.8e-8 | 5.1e-3 | 3.0 | 1.9e-2 / 1.9e-2 |
| | vely / velz | 3.1e-6 | 2.8e-4 | 3.5e-1 | 1.0 / 5.3 |
| | eint / m1_e | 0 | 1.6e-4 | 7.5e-2 | 5.6e-2 / 5.9e-2 |
| | m1_f1 | 7.9e-8 | 8.1e-5 | 8.2e-4 | 2.8e-3 / 2.8e-3 |
| | m1_f2 / f3 | 4.6e-6 | 2.9e-4 | 5.2e-1 | 8.1e-1 / 1.7 |
| sqsbp | dens | 0 | 1.2e-6 | 1.2e-2 | 1.5e-3 / 1.5e-3 |
| | velx | 9.8e-8 | 5.3e-3 | 4.5 | 2.1e-2 / 2.0e-2 |
| | vely / velz | 4.0e-6 | 2.7e-4 | 3.2e-1 | 1.0 / 1.0e1 |
| | eint / m1_e | 0 | 1.7e-4 | 8.2e-2 | 5.6e-2 / 5.9e-2 |
| | m1_f1 | 8.0e-8 | 8.6e-5 | 7.4e-4 | 4.8e-3 / 4.8e-3 |
| | m1_f2 / f3 | 6.0e-6 | 2.8e-4 | 2.0e-1 | 5.2e-1 / 6.8e-1 |
| aq | dens | 0 | 1.2e-6 | 1.4e-2 | 1.5e-3 / 1.5e-3 |
| | velx | 9.8e-8 | 5.0e-3 | 3.1 | 5.1e-2 / 5.3e-2 |
| | vely / velz | 4.2e-6 | 2.8e-4 | 3.4e-1 | 1.0 / **1.9e2** |
| | eint / m1_e | 0 | 1.6e-4 | 8.9e-2 | 9.5e-2 / 1.0e-1 |
| | m1_f1 | 7.9e-8 | 8.0e-5 | 7.1e-4 | 5.6e-2 / 5.9e-2 |
| | m1_f2 / f3 | 6.3e-6 | 2.8e-4 | 2.3e-1 | 1.1 / **1.1e1** |
| aqs | dens | 0 | 1.2e-6 | 1.2e-2 | 1.9e-3 / 1.9e-3 |
| | velx | 9.8e-8 | 5.3e-3 | 4.2 | 1.3e-1 / 1.5e-1 |
| | vely / velz | 3.2e-6 | 2.7e-4 | 3.5e-1 | 1.1 / **5.2e2** |
| | eint / m1_e | 0 | 1.7e-4 | 8.0e-2 | 7.3e-2 / 7.9e-2 |
| | m1_f1 | 7.9e-8 | 8.6e-5 | 8.6e-4 | 1.3e-2 / 1.4e-2 |
| | m1_f2 / f3 | 4.8e-6 | 2.8e-4 | 4.6e-1 | 1.2 / **6.3** |

- The worst entry from each variable pair is shown. The full per-variable output is in RESULTS_gpu_point.txt.
- The >1.05 /all maxima for sq sit in tangential velocity at 2.6-3.0 R_ph: all has max|vely| ~0.65 there, sq has ~70.
- Caltech's metric (/arm) saturates at 1 there, which hides how large this is.

## Binaries

- New: `bin/athena_he_a100_e50ec7d6`, md5 `1dfa47c43c96a0d03ad42a83b96fb36f`.
  - Source: src/ = fork/xthinfix-1009 e50ec7d6 + kokkos 4.6.2.
  - Build: gcc/13, cuda/12.6, openmpi_gpu/5.0, cmake/3.30, nvcc_wrapper, AMPERE80, MPI, Release, PROBLEM=he_star_m1.
  - Build dir `build_gpu`; logs `logs/{cmake,make}_gpu.log`, `logs/build_gpu.31050046.out`.
- Base: `/raven/ptmp/jinma/vgdspeed_1009/bin/athena_he_a100_89693b59_vgdmerge`, md5 `a10415764b328ff5c86893521c6b0ba6`.

## Jobs / scripts

- `gpu/build_gpu.sh`: job 31050046, interactive -c 64.
- `gpu/gpu_arms_raven.sh` (ARMS / NL env), gpudev, 1 node x 4 A100, 4 ranks x 1 MeshBlock:
  - 31050205: smoke, sqsbp, 3 cycles. rc 0, Picard 21 / max 33, pred 22/41, bound 4.9e5, fallbacks 0.
  - 31050237: base_1, all_1, sq, sqs, all_2.
  - 31050238: sqb, sqp, sqsb.
  - 31050239: sqsbp, aq, aqs.
  - 31050588: sqsb_2, sq_2.
- Analysis:
  - `gpu/ana_gpu.py`: log parsing plus a spdiff.py-equivalent with 4 bands and both normalisations. It runs with `/u/jinma/miniconda3/envs/magritte/bin/python`, because the base miniconda and system python lack numpy.
  - `gpu/bitcmp.sh`: bin data after `<par_end>` plus hst data rows.

## Adaptations

- Input: `gpu/inp/agcar_A_xtf.athinput` is a copy of `/raven/ptmp/jinma/rc_1009/inputs/agcar_rcxA_ge_raven_accel_st_mg_hr.athinput`.
  - After `implicit_blend_xthin = 30` it adds every key used on the command line, set to the production / default value:
    - xthin_mode all, xthin_wmin 0, xthin_r0 1.5
    - hr_recon dc, recon_fresh 2, hr_damp 0
    - recon_qs 0, qs_rel 0
    - hr_pos kill, recon_qs_mode pass, recon_pred false
  - Fresh runs make a command-line key that is missing from the input FATAL, which is why these are added.
  - base runs on the unchanged production input, as Caltech does.
- Arm keys are on top of the production hr keys, as in xtf/arms/gpu_arms.sh / battery_arms.sh:
  - sq = mode steep, xthin 60, recon plm, qs 0.1
  - aq = mode all, recon plm, qs 0.1
  - s = `implicit_hr_recon_qs_mode=step`, b = `implicit_hr_pos=bound`, p = `implicit_hr_recon_pred=true`
- Run setup: `time/nlim=30 time/ndiag=1 output5/dt=1e30`, as Caltech.
  - Raven uses 4 ranks x 1 MeshBlock on 4 A100. Caltech used 2 ranks x 2 MeshBlocks on 2 H200.
  - Each arm takes about 25 s of setup outside the timed cycles.
- Timings are single-shot, one arm after another on one node. The repeats agree to 1 %: sq 1.850 vs 1.839, sqsb 1.046 vs 1.047.
- The CPU battery dirs (cpu/, build_cpu) and the held agcA_rc / agcB_rc jobs were not touched.
