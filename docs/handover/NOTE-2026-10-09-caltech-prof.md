# NOTE Caltech -> viper: GPU cost breakdown of the production radiation scheme on H200 (TASK-2026-10-09-caltech-prof)

Done. No DeltaAI/Delta prof NOTE existed at submit, at job start (checked inside the job) or at push.

## Setup
- Job 4297456, node hpc-sm-02-11, **1 node x 2 H200** (2 ranks, 1 GPU per rank, srun --mpi=pmix -c 8 --cpu-bind=cores
  + gpuwrap CUDA_VISIBLE_DEVICES=SLURM_LOCALID), gcc 13.2 / cuda 12.9 / hpcx 2.17.1. Queue wait 55 min, elapsed 4 min 58 s
  for all 8 arms (0.17 GPU-h).
- Binary **athena_gpu_he_star_m1_98835d99_nofma** (98835d99, NOT 6c5d8fb2), md5 bc6f4e9e11c106a3606246d2ac40d296.
- Layout: both cases have 4 MeshBlocks -> **2 MeshBlocks per H200** (DeltaAI: 1 per GH200).
- (a) AG Car B: agcar_rcxB_ge_accel_st_mg_hr.athinput (md5 84eb5faa), fresh t = 0, 480x128x128 = 4 x (480x64x64),
  production outputs kept, 60 cycles.
- (b) BSG: **bsg_hr_dc5.athinput** (md5 be4d08ae, as prof_cuda.sh; not bsg3d_truerepro2_hr_lm, whose
  vet_gd_twin_lowmem key is 6c5d8fb2-only) on the reduced mesh 256x128x128 = 4 x (256x64x64), bin/rst dumps off
  (prof_cuda.sh KB keys unchanged), 30 cycles.
- Scripts: prof_cuda.sh unchanged, driven by prof_caltech.sbatch (SRUN4 override). kokkos-tools master 1e98316
  (2026-10-04), simple-kernel-timer built with g++ 13.2.
- **Tool deviation:** current kokkos-tools master no longer writes `<host>-<pid>.dat` by default; it prints the
  kp_reader-format table to stdout at finalize (a .dat only with KOKKOS_TOOLS_TIMER_BINARY=1). So "kt: 0 rank files"
  and the 2 ranks' tables ended up interleaved (in whole-line chunks) in each kt arm's out.log. I recovered them with
  split_kt_ascii.py (attached): each numbers line pairs with the name line above it, else with the dangling name.
  Checks pass in all 4 kt arms: every kernel has exactly 2 entries, and the sum of entries equals the sum of the two
  ranks' "Total Time in Kokkos kernels". It writes one kt_mean.txt per arm (times = mean over the 2 ranks), which
  prof_group.py reads unchanged. So "mean over 1 rank files" in the outputs = the 2-rank mean, and the "rank spread"
  block is meaningless (max = mean). prof_group.py was not changed. Next time: export KOKKOS_TOOLS_TIMER_BINARY=1 in
  prof_cuda.sh.

## s/cycle (prof_cuda `s/cycle: mean(6-N) median` line), peak GPU memory, health
| arm | s/cycle mean | median | peak/GPU | rc / FATAL / OOM / NaN / NON-CONV |
|---|---|---|---|---|
| agc_plain (60) | **0.5637** | **0.4504** | 40549 MiB | 0/0/0/0/0 |
| agc_kt10 | 1.1143 (2-10) | 0.5472 | 40549 | 0/0/0/0/0 |
| agc_kt60 | 0.6745 | 0.5275 | 40549 | 0/0/0/0/0 |
| agc_tmr | 0.5891 | 0.4511 | 40549 | 0/0/0/0/0 |
| bsg_plain (30) | **0.4615** | **0.2970** | 34591 MiB | 0/0/0/0/0 |
| bsg_kt5 | 1.3509 (2-5) | 0.4871 | 34581 | 0/0/0/0/0 |
| bsg_kt30 | 0.5322 | 0.3548 | 34581 | 0/0/0/0/0 |
| bsg_tmr | 0.4981 | 0.2972 | 34581 | 0/0/0/0/0 |

kt steady (long minus short, fenced): AG Car 0.5647 s/cycle (kernels 0.4672 = 82.7 %); BSG 0.3910 s/cycle (kernels
0.3113 = 79.6 %). The tool adds about +17-20 % to the median. In both cases mean >> median: a few slow cycles, which in
BSG are not bin/rst dumps (those were off). implicit_timers costs about 0 on the median.

## Code counters (plain arms)
- AG Car B: solves 60, **Picard mean 4.12, max 7, NON-CONVERGED 0**; bicgstab outer passes 247, **inner its mean
  3.76, max 17** (total 928), breakdowns 0, 3.53 global reductions per inner it; final linear residual mean 1.3e-12,
  max 7.2e-11 (tol 1e-10).
- BSG reduced: solves 30, **Picard mean 2.70, max 6, NON-CONVERGED 0**; outer passes 81, **inner its mean 26.4,
  max 52** (total 2142), **breakdowns 1**, 3.08 reductions/it; **final linear residual mean 1.8e-9, max 1.1e-8 >
  tol 1e-10** (see Odd).

## <rad_m1> timers (tmr arms, from cycle 3)
- AG Car B: stages 114, solves 57, passes 228, kry_it 835. s: closure 0.0070, opacity 0.145, solve_pre 0.180,
  **tensor 19.81**, pred 0.0097, pass_setup 1.979, krylov 2.247, pass_post 0.248, solve_end 0.279, hyd_c2p 0.160,
  m1_bvals 0.047 (sum 25.1 s of about 33.6 s wall for cycles 3-60).
- BSG: stages 54, solves 27, passes 65, kry_it 1705. s: closure 0.0022, opacity 0.0099, solve_pre 0.0289,
  **tensor 6.567**, pred 0.0034, pass_setup 0.281, krylov 1.500, pass_post 0.039, solve_end 0.070, hyd_c2p 0.021,
  m1_bvals 0.013.
- The "tensor" category (the vet_gd sweep) is 79 % (AG Car) and 77 % (BSG) of the fenced implicit time, matching the
  kernel groups below.

## Grouped (steady cycles, % of kernel time; s/cycle per rank)
| group | AG Car B | BSG reduced |
|---|---|---|
| vet_gd halo pack/unpack (m1_vgd_hc_*, unpack) | 34.1 % (0.160) | 41.9 % (0.130) |
| vet_gd sweep (m1_vgd_shell, mom) | 31.4 % (0.147) | 24.5 % (0.076) |
| vet_gd twin | 1.9 % | 2.8 % |
| **vet_gd total** | **67.5 %** | **69.1 %** |
| half-range pass | 2.2 % | 2.5 % |
| vet_col/lat | 3.1 % | 1.1 % |
| rad matrix/operator assembly | 4.0 % | 2.4 % |
| rad mg/line preconditioner (pcrx, mg_pcr) | 4.0 % | 7.4 % |
| rad Krylov/BiCGStab vector ops | 0.8 % | 2.5 % |
| rad Picard update/residual/accel (impl_stop ...) | 2.4 % | 5.3 % |
| rad implicit halo pack/unpack | 0.5 % | 1.0 % |
| unmapped (m1_impl_asm, m1_impl_tcell, m1_opacity, f2/f3face, RecvBuff, ...) | 6.3 % | 3.0 % |
| rad other | 0.2 % | 0.8 % |
| hydro fluxes/recon | 2.6 % | 2.1 % |
| hydro c2p (eos) | 2.2 % | 0.3 % |
| hydro update/other | 1.7 % | 0.5 % |
| halo/MPI pack-unpack (bvals) | 0.02 % | 0.02 % |
| outputs | 0.02 % | 0.01 % |
| srcterms+pgen / mesh+utils / kokkos internal | 2.4 % | 2.2 % |
| outside kernels (host + MPI waits + launch gaps), % of wall | 17.3 % | 20.4 % |

**Radiation vs hydro:** AG Car B: radiation 84.7 % of kernel time, hydro 6.6 %, rad/(rad+hydro) = **92.8 %**
(about 99 % if the unmapped m1_* labels are counted as radiation). BSG: radiation 92.1 %, hydro 2.8 %,
rad/(rad+hydro) = **97.1 %**.

## Top 25 kernels: AG Car B (steady 50 cycles, 2-rank mean)
```
# STEADY = agc_kt60 minus agc_kt10 (50 cycles), mean over 1 rank files
# wall (Total Execution Time) 28.237 s, in kernels 23.362 s (82.7 %), outside kernels (host, MPI waits, launch gaps) 4.874 s
# s/cycle (tool run, fenced): 0.5647, kernels 0.4672

## top 25 kernels (time, % of kernel time, % of wall, calls)
m1_vgd_shell                           7.066 s  30.24 %  25.02 %    24600  vet_gd sweep (other)   0.1413 s/cyc
m1_vgd_unpack                          2.649 s  11.34 %   9.38 %    24650  vet_gd halo pack/unpack   0.0530 s/cyc
m1_vgd_hc_scan                         2.560 s  10.96 %   9.07 %    49200  vet_gd halo pack/unpack   0.0512 s/cyc
m1_vgd_hc_flag                         1.445 s   6.19 %   5.12 %    24600  vet_gd halo pack/unpack   0.0289 s/cyc
m1_vgd_hc_pack                         0.671 s   2.87 %   2.38 %    24600  vet_gd halo pack/unpack   0.0134 s/cyc
m1_impl_asm                            0.526 s   2.25 %   1.86 %      200  unmapped (label not a literal in src/)   0.0105 s/cyc
hyd_c2p_gen                            0.518 s   2.22 %   1.83 %      350  hydro c2p (eos)   0.0104 s/cyc
m1_impl_pcrx                           0.478 s   2.05 %   1.69 %     2860  rad mg/line preconditioner   0.0096 s/cyc
m1_vgd_tw_sumd                         0.433 s   1.85 %   1.53 %      150  vet_gd twin   0.0087 s/cyc
m1_vgd_hc_expand                       0.415 s   1.78 %   1.47 %    24600  vet_gd halo pack/unpack   0.0083 s/cyc
m1_impl_stop                           0.378 s   1.62 %   1.34 %     3654  rad Picard update/residual/accel   0.0076 s/cyc
m1_impl_tcell                          0.363 s   1.55 %   1.29 %      250  unmapped (label not a literal in src/)   0.0073 s/cyc
m1_vlat_src                            0.294 s   1.26 %   1.04 %       50  vet_col/lat   0.0059 s/cyc
wbcache                                0.286 s   1.22 %   1.01 %      100  hydro update/other   0.0057 s/cyc
m1_vgd_hr_sm                           0.273 s   1.17 %   0.97 %      300  half-range pass   0.0055 s/cyc
m1_vgd_mom                             0.266 s   1.14 %   0.94 %      100  vet_gd sweep (other)   0.0053 s/cyc
m1_vlat_stencil                        0.237 s   1.02 %   0.84 %      200  vet_col/lat   0.0047 s/cyc
m1_vgd_hr                              0.228 s   0.97 %   0.81 %      100  half-range pass   0.0046 s/cyc
hflux_x2                               0.225 s   0.97 %   0.80 %      100  hydro fluxes/recon   0.0045 s/cyc
hflux_x3                               0.222 s   0.95 %   0.79 %      100  hydro fluxes/recon   0.0044 s/cyc
m1_vgd_hc_bnd                          0.222 s   0.95 %   0.79 %    24600  vet_gd halo pack/unpack   0.0044 s/cyc
m1_impl_hdir                           0.207 s   0.88 %   0.73 %     4057  rad matrix/operator assembly   0.0041 s/cyc
m1_mg_pcr                              0.196 s   0.84 %   0.70 %     5720  rad mg/line preconditioner   0.0039 s/cyc
m1_vlat_op                             0.168 s   0.72 %   0.60 %      200  vet_col/lat   0.0034 s/cyc
hflux_x1                               0.165 s   0.71 %   0.59 %      100  hydro fluxes/recon   0.0033 s/cyc

```

## Top 25 kernels: BSG reduced (steady 25 cycles, 2-rank mean)
```
# STEADY = bsg_kt30 minus bsg_kt5 (25 cycles), mean over 1 rank files
# wall (Total Execution Time) 9.774 s, in kernels 7.781 s (79.6 %), outside kernels (host, MPI waits, launch gaps) 1.993 s
# s/cycle (tool run, fenced): 0.3910, kernels 0.3113

## top 25 kernels (time, % of kernel time, % of wall, calls)
m1_vgd_shell                           1.839 s  23.63 %  18.81 %     5250  vet_gd sweep (other)   0.0735 s/cyc
m1_vgd_unpack                          1.285 s  16.52 %  13.15 %     5275  vet_gd halo pack/unpack   0.0514 s/cyc
m1_vgd_hc_scan                         0.965 s  12.40 %   9.87 %    10500  vet_gd halo pack/unpack   0.0386 s/cyc
m1_vgd_hc_flag                         0.583 s   7.49 %   5.96 %     5250  vet_gd halo pack/unpack   0.0233 s/cyc
m1_impl_pcrx                           0.578 s   7.43 %   5.91 %     5952  rad mg/line preconditioner   0.0231 s/cyc
m1_impl_stop                           0.375 s   4.81 %   3.83 %     6184  rad Picard update/residual/accel   0.0150 s/cyc
m1_vgd_hc_pack                         0.227 s   2.91 %   2.32 %     5250  vet_gd halo pack/unpack   0.0091 s/cyc
m1_vgd_tw_sumd                         0.213 s   2.73 %   2.18 %       75  vet_gd twin   0.0085 s/cyc
m1_vgd_hc_expand                       0.146 s   1.87 %   1.49 %     5250  vet_gd halo pack/unpack   0.0058 s/cyc
m1_vgd_hr_sm                           0.134 s   1.72 %   1.37 %      150  half-range pass   0.0054 s/cyc
m1_impl_bcgf_upd                       0.106 s   1.37 %   1.09 %     1488  rad Krylov/BiCGStab vector ops   0.0043 s/cyc
m1_impl_hdir                           0.095 s   1.22 %   0.97 %     3406  rad matrix/operator assembly   0.0038 s/cyc
m1_impl_asm                            0.080 s   1.02 %   0.81 %       57  unmapped (label not a literal in src/)   0.0032 s/cyc
Kokkos::View::initialization [_mir     0.076 s   0.98 %   0.78 %       12  kokkos internal (view init, copies)   0.0030 s/cyc
m1_impl_tcell                          0.067 s   0.87 %   0.69 %       82  unmapped (label not a literal in src/)   0.0027 s/cyc
hflux_x2                               0.065 s   0.84 %   0.67 %       50  hydro fluxes/recon   0.0026 s/cyc
hflux_x3                               0.065 s   0.83 %   0.66 %       50  hydro fluxes/recon   0.0026 s/cyc
m1_vgd_mom                             0.060 s   0.77 %   0.61 %       50  vet_gd sweep (other)   0.0024 s/cyc
m1_det_red2                            0.054 s   0.70 %   0.55 %     4580  rad other (opacity, explicit flux, time2, coupling, dt)   0.0022 s/cyc
m1_vgd_hr                              0.053 s   0.68 %   0.54 %       50  half-range pass   0.0021 s/cyc
m1_vgd_hc_bnd                          0.048 s   0.62 %   0.49 %     5250  vet_gd halo pack/unpack   0.0019 s/cyc
m1_impl_bcgf_ts                        0.042 s   0.54 %   0.43 %     1488  rad Krylov/BiCGStab vector ops   0.0017 s/cyc
m1_impl_hm_pack                        0.039 s   0.50 %   0.40 %     3406  rad implicit halo pack/unpack   0.0016 s/cyc
m1_impl_hm_unpack                      0.039 s   0.50 %   0.39 %     3406  rad implicit halo pack/unpack   0.0015 s/cyc
m1_vlat_stencil                        0.037 s   0.47 %   0.38 %       57  vet_col/lat   0.0015 s/cyc

```

## Odd / worth a look
1. **The vet_gd halo exchange is the biggest single cost, not the sweep.** m1_vgd_unpack + hc_scan + hc_flag +
   hc_pack + hc_expand + hc_bnd are 34-42 % of kernel time, with about 490 (AG Car) / 210 (BSG) launches per cycle of
   each, about 50-110 us per call. That looks launch/latency bound and is the obvious lever: fuse hc_flag/scan/pack,
   or exchange less often. m1_vgd_shell alone is 24-30 %.
2. BSG reduced: final 7-point linear residual mean 1.8e-9 / max 1.1e-8 above lin_tol 1e-10, plus 1 bicgstab
   breakdown, with NON-CONVERGED 0 (Picard converged). AG Car meets tol. inner its 26 (BSG) vs 3.8 (AG Car).
3. Labels missing from labels_98835d99.tsv (m1_impl_asm, m1_impl_tcell, m1_opacity, m1_impl_f2face/f3face,
   m1_impl_i1, m1_impl_tsolve) fall into "unmapped"; they are all rad_m1 implicit assembly/opacity work (6.3 % AG Car).
4. 2 blocks per GPU here vs 1 on DeltaAI, so per-launch overheads differ. Not directly comparable 1:1 with GH200.

## Files
Attached under docs/handover/caltech-prof-1009/results/4297456/: agc_prof.txt, bsg_prof.txt, kt_mean_<arm>.txt (the
recovered 2-rank-mean tables), out.log tails (tail -60) of the plain and tmr arms, prof_caltech.sbatch,
split_kt_ascii.py. Run dirs are kept at /resnick/groups/carnegie_poc/jingze/prof_1009/runs/4297456 (scratch, purged
after 14 d).
