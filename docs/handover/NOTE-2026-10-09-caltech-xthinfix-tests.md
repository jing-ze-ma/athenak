# NOTE Caltech -> viper: xthinfix-1009 test battery (f1d355e8d), 10-09

Answers TASK-2026-10-09-caltech-xthinfix-tests.md. Run trees stay on Caltech:
/resnick/groups/carnegie_poc/jingze/xthinfix_1009 (run/, logs/; REPO = git archive of f1d355e8d in repo/).

STATUS: complete (GPU 4299983 + spdiff 4300415; CPU order 4299601-04, gates+beams 4300406, ana 4299605 + re-ana 4302028).

## Verdict (hrb = mode beam)

1. Smooth problems: hrb == cen BITWISE on all pulses (arm_minus_cen 0.00, every kappa, be/hesdirk2, X/T/C); rw_t10 E converges (hesdirk2 X 1.91 1.97 1.69, = hrx0; hr -1.67 -1.18); rw_t1000 hrb-cen <= 4e-7.
2. atm (Hopf L1 32..512): hrb 5.8e-3 1.7e-3 8.3e-4 7.0e-4 6.1e-4 (monotone, hr 5.5e-3 ... 4.2e-3), but flattens like hrx0 (3.7e-4 at 512); 512 is 1.6x hrx0. Gates G1 (1.69e-2), G5 (1280: 1.017, 12800: 1.000), G3 (L1 0.0138) PASS for hr and hrb (G5 kappa 12.8/128 ratio ~0.55 for both arms; not in the PASS list).
3. Beams: no FATAL, NON-CONVERGED 0, but hrb is NOT comparable to hr in cost: BiCGStab inner mean 172-178 on xb20/ba0/ba20 (135 shd3b, 51 cyl; max = 200) vs hr 17-42, Picard 8.6-12 vs 4-5.4 on xb20/ba0/ba20, wall 2-8x; and 20-28 % of ALL faces have |F|/cE > 1 (max 580-970) in xb20/ba0/ba20 (hr 0-0.2 %). The beams themselves match hr (cyl median 4.6 deg, xb20 direction 4.8 deg, shd3b umbra 0 / edges 0.232 / 0.346 identical).
4. AG Car A on GPU: rc 0, NON-CONVERGED 0, Picard identical; s/cycle beam 0.627/0.628 vs all 0.655/0.624 (noise), 1 GPU 0.816 vs 0.810; repeats bitwise (6/6 files); all-vs-beam <= 4e-6 below 0.8 R_ph (not exactly 0), 1e-5..8e-4 at 0.8-0.95, up to 0.33 (vely/velz) at 0.95-1.3 R_ph.
5. Net: the fix restores 2nd order / central behaviour on smooth problems, but in the beam tests the dark-region faces become superluminal and the linear solve stalls at lin_maxit; mode beam not ready as is.

## Caltech adaptations (bundle issues found)

- Gate inputs (atm2d, generated pulse2d, marsh2d) lack `vet_tensor`, so every gate run was FATAL ("Parameter 'vet_tensor' in block 'rad_m1' on command line not found"). Fixed in a copy (battery_cpu_caltech.sh + fixkeys.py): any command-line key missing from the input is added to it (only vet_tensor was missing; beam inputs needed nothing). Bundle fix: add `vet_tensor = full` to those inputs.
- beams() uses /usr/bin/time (absent on Caltech compute nodes): replaced by a shell stand-in with the same "%e s wall" output.
- cylref.py / collref.py read `<run>/cmd.txt`, which beams() never writes (Traceback in the first ana pass). Wrote cmd.txt (= ATHENA_CPU + the beam's command-line keys, c_light=1000) afterwards and re-ran ana (4302028). Bundle fix: write cmd.txt in beams().
- python = spack python 3.11 + numpy/scipy/matplotlib + py-h5py (bin_convert needs h5py). order was split round-robin into 4 jobs x 32 cores (od.py run on 4 sub-job-files); no other changes.
- No FATAL in any order/gate/beam run log; all 860 order runs rc=0; NON-CONVERGED=0 in all beam runs.

## Binaries (f1d355e8d70490c79fd736df2af1f9d3a4bb9dbb)

- ATHENA_GPU (H200, he_star_m1, CUDA 12.9, HOPPER90, MPI hpcx, host -ffp-contract=off, incremental in the nofma tree of
  98835d99, 40 objects recompiled): `athena_gpu_he_star_m1_f1d355e8_nofma` md5 `4f9c64f98000fd443626b1550612465b`
- ATHENA_CPU (no PROBLEM = built_in_pgens, MPI openmpi 5.0.1, Release, gcc 13.2): `athena_cpu_built_in_pgens_f1d355e8`
  md5 `fed1609255749fc2fb1a497cbd60d038`

## GPU (job 4299983, 2 H200, hpc-sm-02-11; smoke 4299546 nlim 3 beam: rc 0, NON-CONVERGED 0)

Caltech adaptations only: #SBATCH header, srun `--mpi=pmix -n N -c 8 --cpu-bind=cores gpuwrap.sh` (CUDA_VISIBLE_DEVICES =
SLURM_LOCALID). spdiff.py needs h5py (via vis/python/bin_convert); the in-job python lacked it, so the last section was
re-run afterwards in CPU job 4300415 with py-h5py (same files). All 6 runs rc=0 fatal=0, 30 cycles.

RESULTS_gpu.txt (raw):
```
all_1    0.655 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  FATAL 0
beam_1   0.627 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.353846e+00  FATAL 0
all_2    0.624 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  FATAL 0
beam_2   0.628 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.353846e+00  FATAL 0
beam_1g  0.816 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.353846e+00  FATAL 0
all_1g   0.810 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  FATAL 0
-- bitwise repeats (hst+bin): all_1 vs all_2, beam_1 vs beam_2
all_1 vs all_2: 6 files, 0 differ
beam_1 vs beam_2: 6 files, 0 differ
-- all vs beam by radius (max over angles of |a-b|/max|b|), last dump
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 2.9e-07 1.3e-03 7.6e-04 2.0e-07 5.5e-07
agcar3d.hydro_w.00001.bin    velx       9.8e-08 5.0e-08 8.3e-04 1.3e-02 5.9e-03 5.7e-05 4.0e-05
agcar3d.hydro_w.00001.bin    vely       2.9e-06 3.8e-08 1.8e-05 1.8e-01 3.3e-01 1.1e-03 9.7e-04
agcar3d.hydro_w.00001.bin    velz       2.8e-06 6.5e-09 1.6e-05 1.6e-01 3.1e-01 1.5e-03 1.4e-03
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 6.3e-06 1.3e-03 7.6e-04 5.0e-05 4.4e-05
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 2.4e-05 9.7e-04 6.5e-05 3.0e-05 1.4e-05
agcar3d.m1.00001.bin         m1_f1      8.0e-08 5.0e-09 1.2e-05 1.5e-05 8.0e-06 8.2e-06 8.1e-06
agcar3d.m1.00001.bin         m1_f2      4.4e-06 1.8e-08 1.1e-05 1.2e-03 7.5e-04 7.6e-04 7.1e-04
agcar3d.m1.00001.bin         m1_f3      4.1e-06 2.8e-08 1.2e-05 1.0e-03 6.9e-04 7.7e-04 7.6e-04
```

## RESULTS/order_eval.txt (raw)

```
atm          cen  be       S e      lev 32-64-128-256-512  L1 1.82e-03 4.59e-04 1.18e-04 3.86e-05 | ord  1.99  1.96  1.61
atm          cen  be       S f1     lev 32-64-128-256-512  L1 9.17e-01 7.37e-01 1.12e+00 5.44e+00 | ord  0.32 -0.60 -2.28
atm          hr   be       S e      lev 32-64-128-256-512  L1 2.11e-03 9.83e-04 8.98e-04 1.10e-03 | ord  1.10  0.13 -0.29
atm          hr   be       S f1     lev 32-64-128-256-512  L1 1.27e+00 1.91e+01 5.60e+01 9.75e-01 | ord -3.91 -1.55  5.85
atm          hrb  be       S e      lev 32-64-128-256-512  L1 1.81e-03 4.67e-04 1.65e-04 1.30e-04 | ord  1.96  1.50  0.34
atm          hrb  be       S f1     lev 32-64-128-256-512  L1 1.27e+00 1.35e+01 3.12e+01 1.22e+00 | ord -3.41 -1.20  4.67
atm          hrx0 be       S e      lev 32-64-128-256-512  L1 1.77e-03 4.41e-04 1.11e-04 3.91e-05 | ord  2.00  1.99  1.51
atm          hrx0 be       S f1     lev 32-64-128-256-512  L1 1.62e+00 1.02e+01 3.29e+01 9.05e-01 | ord -2.65 -1.69  5.18
pulse_k0.128 cen  be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 cen  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 cen  be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 cen  be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 cen  be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 cen  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 cen  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 cen  hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 cen  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 hr   be       X e      lev 32-64-128-256-512  L1 7.25e-02 2.14e-02 6.33e-03 1.11e-02 | ord  1.76  1.76 -0.82
pulse_k0.128 hr   be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.91e-02 7.96e-03 1.71e-02 | ord  1.83  1.87 -1.10
pulse_k0.128 hr   be       T e      lev 1-2-4-8-16  L1 5.35e-02 2.81e-02 1.44e-02 7.33e-03 | ord  0.93  0.96  0.98
pulse_k0.128 hr   be       T f1     lev 1-2-4-8-16  L1 6.54e-02 3.51e-02 1.83e-02 9.31e-03 | ord  0.90  0.94  0.97
pulse_k0.128 hr   be       C e      lev 32-64-128-256-512  L1 8.16e-02 2.65e-02 1.07e-02 4.61e-03 | ord  1.62  1.32  1.21
pulse_k0.128 hr   be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.61e-02 1.90e-02 8.41e-03 | ord  1.33  1.28  1.17
pulse_k0.128 hr   hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.97e-03 3.70e-03 | ord  1.74  1.90  0.69
pulse_k0.128 hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.88e-03 5.25e-03 | ord  1.83  1.96  0.59
pulse_k0.128 hr   hesdirk2 T e      lev 1-2-4-8-16  L1 5.83e-03 2.47e-03 1.13e-03 5.42e-04 | ord  1.24  1.13  1.06
pulse_k0.128 hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 7.47e-03 3.40e-03 1.64e-03 8.11e-04 | ord  1.14  1.05  1.02
pulse_k0.128 hr   hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.93e-03 1.80e-03 | ord  1.75  1.93  1.72
pulse_k0.128 hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.95e-03 2.23e-03 | ord  1.83  1.96  1.83
pulse_k0.128 hrb  be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 hrb  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 hrb  be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 hrb  be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 hrb  be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 hrb  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 hrb  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 hrb  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 hrb  hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 hrb  hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 hrb  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 hrb  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 hrx0 be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 hrx0 be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 hrx0 be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 hrx0 be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 hrx0 be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 hrx0 be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 hrx0 hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 hrx0 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 hrx0 hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 hrx0 hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 hrx0 hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 hrx0 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k12.8  cen  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  cen  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  cen  be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  cen  be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  cen  be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  cen  be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  cen  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  cen  hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  cen  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  hr   be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.36e-03 1.41e-03 2.05e-03 | ord  1.89  2.18 -0.55
pulse_k12.8  hr   be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.62e-02 3.49e-03 4.78e-03 | ord  1.78  2.21 -0.45
pulse_k12.8  hr   be       T e      lev 1-2-4-8-16  L1 1.54e-02 8.30e-03 4.32e-03 2.21e-03 | ord  0.89  0.94  0.97
pulse_k12.8  hr   be       T f1     lev 1-2-4-8-16  L1 3.43e-02 1.91e-02 1.01e-02 5.19e-03 | ord  0.85  0.92  0.96
pulse_k12.8  hr   be       C e      lev 32-64-128-256-512  L1 2.80e-02 1.15e-02 4.24e-03 1.56e-03 | ord  1.28  1.44  1.44
pulse_k12.8  hr   be       C f1     lev 32-64-128-256-512  L1 7.53e-02 3.17e-02 1.23e-02 5.06e-03 | ord  1.25  1.37  1.28
pulse_k12.8  hr   hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.38e-03 6.57e-04 | ord  1.90  2.24  1.07
pulse_k12.8  hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.67e-02 3.65e-03 1.37e-03 | ord  1.76  2.20  1.41
pulse_k12.8  hr   hesdirk2 T e      lev 1-2-4-8-16  L1 4.39e-03 2.16e-03 1.07e-03 5.30e-04 | ord  1.03  1.02  1.01
pulse_k12.8  hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 9.26e-03 4.58e-03 2.27e-03 1.13e-03 | ord  1.02  1.01  1.01
pulse_k12.8  hr   hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.34e-03 1.86e-03 4.83e-04 | ord  1.69  1.98  1.94
pulse_k12.8  hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.26e-03 1.56e-03 | ord  1.63  1.91  1.76
pulse_k12.8  hrb  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  hrb  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  hrb  be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  hrb  be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  hrb  be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  hrb  be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  hrb  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  hrb  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  hrb  hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  hrb  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  hrb  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  hrb  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  hrx0 be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  hrx0 be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  hrx0 be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  hrx0 be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  hrx0 be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  hrx0 be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  hrx0 hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  hrx0 hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  hrx0 hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  hrx0 hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  hrx0 hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  hrx0 hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k128   cen  be       X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   cen  be       X f1     lev 32-64-128-256-512  L1 2.55e-02 6.91e-03 1.68e-03 4.02e-04 | ord  1.89  2.04  2.06
pulse_k128   cen  be       T e      lev 1-2-4-8-16  L1 1.46e-03 7.37e-04 3.70e-04 1.85e-04 | ord  0.99  0.99  1.00
pulse_k128   cen  be       T f1     lev 1-2-4-8-16  L1 5.40e-03 2.73e-03 1.37e-03 6.86e-04 | ord  0.99  0.99  1.00
pulse_k128   cen  be       C e      lev 32-64-128-256-512  L1 1.78e-02 5.18e-03 1.42e-03 4.37e-04 | ord  1.79  1.86  1.70
pulse_k128   cen  be       C f1     lev 32-64-128-256-512  L1 2.51e-02 7.99e-03 2.49e-03 9.03e-04 | ord  1.65  1.68  1.46
pulse_k128   cen  hesdirk2 X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.55e-02 6.92e-03 1.68e-03 4.02e-04 | ord  1.88  2.04  2.06
pulse_k128   cen  hesdirk2 T e      lev 1-2-4-8-16  L1 9.44e-05 3.66e-05 1.58e-05 7.27e-06 | ord  1.37  1.21  1.12
pulse_k128   cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 2.78e-04 1.05e-04 4.48e-05 2.06e-05 | ord  1.40  1.23  1.12
pulse_k128   cen  hesdirk2 C e      lev 32-64-128-256-512  L1 1.66e-02 4.49e-03 1.07e-03 2.60e-04 | ord  1.88  2.07  2.04
pulse_k128   cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.56e-02 6.96e-03 1.71e-03 4.15e-04 | ord  1.88  2.03  2.04
pulse_k128   hr   be       X e      lev 32-64-128-256-512  L1 1.64e-02 4.46e-03 1.08e-03 3.00e-04 | ord  1.88  2.05  1.85
pulse_k128   hr   be       X f1     lev 32-64-128-256-512  L1 2.55e-02 6.93e-03 1.71e-03 4.66e-04 | ord  1.88  2.02  1.88
pulse_k128   hr   be       T e      lev 1-2-4-8-16  L1 1.42e-03 6.98e-04 3.43e-04 1.69e-04 | ord  1.02  1.03  1.02
pulse_k128   hr   be       T f1     lev 1-2-4-8-16  L1 5.32e-03 2.66e-03 1.33e-03 6.62e-04 | ord  1.00  1.00  1.00
pulse_k128   hr   be       C e      lev 32-64-128-256-512  L1 1.79e-02 5.19e-03 1.43e-03 4.37e-04 | ord  1.78  1.86  1.71
pulse_k128   hr   be       C f1     lev 32-64-128-256-512  L1 2.51e-02 8.00e-03 2.49e-03 9.00e-04 | ord  1.65  1.68  1.47
pulse_k128   hr   hesdirk2 X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.07e-03 2.71e-04 | ord  1.88  2.06  1.98
pulse_k128   hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 2.55e-02 6.92e-03 1.69e-03 4.19e-04 | ord  1.88  2.03  2.01
pulse_k128   hr   hesdirk2 T e      lev 1-2-4-8-16  L1 7.37e-05 3.12e-05 1.72e-05 9.18e-06 | ord  1.24  0.86  0.90
pulse_k128   hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 2.62e-04 1.19e-04 6.23e-05 3.24e-05 | ord  1.14  0.94  0.95
pulse_k128   hr   hesdirk2 C e      lev 32-64-128-256-512  L1 1.66e-02 4.49e-03 1.07e-03 2.59e-04 | ord  1.88  2.07  2.05
pulse_k128   hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 2.56e-02 6.97e-03 1.71e-03 4.11e-04 | ord  1.88  2.03  2.05
pulse_k128   hrb  be       X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrb  be       X f1     lev 32-64-128-256-512  L1 2.55e-02 6.91e-03 1.68e-03 4.02e-04 | ord  1.89  2.04  2.06
pulse_k128   hrb  be       T e      lev 1-2-4-8-16  L1 1.46e-03 7.37e-04 3.70e-04 1.85e-04 | ord  0.99  0.99  1.00
pulse_k128   hrb  be       T f1     lev 1-2-4-8-16  L1 5.40e-03 2.73e-03 1.37e-03 6.86e-04 | ord  0.99  0.99  1.00
pulse_k128   hrb  be       C e      lev 32-64-128-256-512  L1 1.78e-02 5.18e-03 1.42e-03 4.37e-04 | ord  1.79  1.86  1.70
pulse_k128   hrb  be       C f1     lev 32-64-128-256-512  L1 2.51e-02 7.99e-03 2.49e-03 9.03e-04 | ord  1.65  1.68  1.46
pulse_k128   hrb  hesdirk2 X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrb  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.55e-02 6.92e-03 1.68e-03 4.02e-04 | ord  1.88  2.04  2.06
pulse_k128   hrb  hesdirk2 T e      lev 1-2-4-8-16  L1 9.44e-05 3.66e-05 1.58e-05 7.27e-06 | ord  1.37  1.21  1.12
pulse_k128   hrb  hesdirk2 T f1     lev 1-2-4-8-16  L1 2.78e-04 1.05e-04 4.48e-05 2.06e-05 | ord  1.40  1.23  1.12
pulse_k128   hrb  hesdirk2 C e      lev 32-64-128-256-512  L1 1.66e-02 4.49e-03 1.07e-03 2.60e-04 | ord  1.88  2.07  2.04
pulse_k128   hrb  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.56e-02 6.96e-03 1.71e-03 4.15e-04 | ord  1.88  2.03  2.04
pulse_k128   hrx0 be       X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrx0 be       X f1     lev 32-64-128-256-512  L1 2.55e-02 6.91e-03 1.68e-03 4.02e-04 | ord  1.89  2.04  2.06
pulse_k128   hrx0 be       T e      lev 1-2-4-8-16  L1 1.46e-03 7.37e-04 3.70e-04 1.85e-04 | ord  0.99  0.99  1.00
pulse_k128   hrx0 be       T f1     lev 1-2-4-8-16  L1 5.40e-03 2.73e-03 1.37e-03 6.86e-04 | ord  0.99  0.99  1.00
pulse_k128   hrx0 be       C e      lev 32-64-128-256-512  L1 1.78e-02 5.18e-03 1.42e-03 4.37e-04 | ord  1.79  1.86  1.70
pulse_k128   hrx0 be       C f1     lev 32-64-128-256-512  L1 2.51e-02 7.99e-03 2.49e-03 9.03e-04 | ord  1.65  1.68  1.46
pulse_k128   hrx0 hesdirk2 X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrx0 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.55e-02 6.92e-03 1.68e-03 4.02e-04 | ord  1.88  2.04  2.06
pulse_k128   hrx0 hesdirk2 T e      lev 1-2-4-8-16  L1 9.44e-05 3.66e-05 1.58e-05 7.27e-06 | ord  1.37  1.21  1.12
pulse_k128   hrx0 hesdirk2 T f1     lev 1-2-4-8-16  L1 2.78e-04 1.05e-04 4.48e-05 2.06e-05 | ord  1.40  1.23  1.12
pulse_k128   hrx0 hesdirk2 C e      lev 32-64-128-256-512  L1 1.66e-02 4.49e-03 1.07e-03 2.60e-04 | ord  1.88  2.07  2.04
pulse_k128   hrx0 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.56e-02 6.96e-03 1.71e-03 4.15e-04 | ord  1.88  2.03  2.04
pulse_k1280  cen  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  cen  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  cen  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  cen  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  cen  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  cen  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  cen  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  cen  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  cen  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  hr   be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 1.01e-03 2.92e-04 | ord  1.89  1.95  1.80
pulse_k1280  hr   be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.08e-03 1.60e-03 4.69e-04 | ord  1.93  1.93  1.77
pulse_k1280  hr   be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.50e-04 1.75e-04 | ord  0.99  1.00  1.00
pulse_k1280  hr   be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.30e-03 6.52e-04 | ord  0.99  0.99  1.00
pulse_k1280  hr   be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.36e-03 4.60e-04 | ord  1.76  1.76  1.56
pulse_k1280  hr   be       C f1     lev 32-64-128-256-512  L1 2.22e-02 6.84e-03 2.31e-03 9.23e-04 | ord  1.70  1.57  1.32
pulse_k1280  hr   hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 1.01e-03 2.89e-04 | ord  1.89  1.95  1.81
pulse_k1280  hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.08e-03 1.60e-03 4.64e-04 | ord  1.93  1.93  1.78
pulse_k1280  hr   hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.10e-05 1.87e-06 3.67e-06 | ord  2.08  2.55 -0.97
pulse_k1280  hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.14e-05 1.23e-05 1.36e-05 | ord  2.01  1.75 -0.14
pulse_k1280  hr   hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.94e-03 1.01e-03 2.82e-04 | ord  1.89  1.96  1.84
pulse_k1280  hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.09e-03 1.60e-03 4.50e-04 | ord  1.92  1.93  1.83
pulse_k1280  hrb  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrb  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrb  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrb  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrb  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  hrb  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  hrb  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrb  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrb  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  hrb  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  hrb  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrb  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  hrx0 be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrx0 be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrx0 be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrx0 be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrx0 be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  hrx0 be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  hrx0 hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrx0 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrx0 hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  hrx0 hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  hrx0 hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrx0 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k12800 cen  be       X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 cen  be       X f1     lev 32-64-128-256-512  L1 2.30e-02 5.88e-03 1.48e-03 3.72e-04 | ord  1.97  1.99  1.99
pulse_k12800 cen  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k12800 cen  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k12800 cen  be       C e      lev 32-64-128-256-512  L1 1.54e-02 4.48e-03 1.28e-03 4.09e-04 | ord  1.79  1.80  1.65
pulse_k12800 cen  be       C f1     lev 32-64-128-256-512  L1 2.19e-02 6.59e-03 2.16e-03 8.28e-04 | ord  1.73  1.61  1.38
pulse_k12800 cen  hesdirk2 X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.29e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 cen  hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.17e-05 2.95e-06 7.45e-07 | ord  1.98  1.99  1.99
pulse_k12800 cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.21e-05 1.06e-05 2.68e-06 | ord  1.98  1.99  1.99
pulse_k12800 cen  hesdirk2 C e      lev 32-64-128-256-512  L1 1.44e-02 3.80e-03 9.41e-04 2.37e-04 | ord  1.92  2.01  1.99
pulse_k12800 cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.28e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hr   be       X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.39e-04 2.39e-04 | ord  1.92  2.01  1.97
pulse_k12800 hr   be       X f1     lev 32-64-128-256-512  L1 2.30e-02 5.88e-03 1.48e-03 3.77e-04 | ord  1.97  1.99  1.97
pulse_k12800 hr   be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k12800 hr   be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k12800 hr   be       C e      lev 32-64-128-256-512  L1 1.54e-02 4.48e-03 1.29e-03 4.12e-04 | ord  1.78  1.80  1.64
pulse_k12800 hr   be       C f1     lev 32-64-128-256-512  L1 2.19e-02 6.59e-03 2.16e-03 8.32e-04 | ord  1.73  1.61  1.38
pulse_k12800 hr   hesdirk2 X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.39e-04 2.39e-04 | ord  1.92  2.01  1.97
pulse_k12800 hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 2.29e-02 5.89e-03 1.48e-03 3.78e-04 | ord  1.96  1.99  1.97
pulse_k12800 hr   hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.17e-05 2.95e-06 7.40e-07 | ord  1.98  1.99  1.99
pulse_k12800 hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.21e-05 1.06e-05 2.67e-06 | ord  1.98  1.99  1.99
pulse_k12800 hr   hesdirk2 C e      lev 32-64-128-256-512  L1 1.44e-02 3.80e-03 9.42e-04 2.40e-04 | ord  1.92  2.01  1.97
pulse_k12800 hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 2.28e-02 5.88e-03 1.48e-03 3.78e-04 | ord  1.96  1.99  1.97
pulse_k12800 hrb  be       X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrb  be       X f1     lev 32-64-128-256-512  L1 2.30e-02 5.88e-03 1.48e-03 3.72e-04 | ord  1.97  1.99  1.99
pulse_k12800 hrb  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrb  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrb  be       C e      lev 32-64-128-256-512  L1 1.54e-02 4.48e-03 1.28e-03 4.09e-04 | ord  1.79  1.80  1.65
pulse_k12800 hrb  be       C f1     lev 32-64-128-256-512  L1 2.19e-02 6.59e-03 2.16e-03 8.28e-04 | ord  1.73  1.61  1.38
pulse_k12800 hrb  hesdirk2 X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrb  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.29e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrb  hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.17e-05 2.95e-06 7.45e-07 | ord  1.98  1.99  1.99
pulse_k12800 hrb  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.21e-05 1.06e-05 2.68e-06 | ord  1.98  1.99  1.99
pulse_k12800 hrb  hesdirk2 C e      lev 32-64-128-256-512  L1 1.44e-02 3.80e-03 9.41e-04 2.37e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrb  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.28e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrx0 be       X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrx0 be       X f1     lev 32-64-128-256-512  L1 2.30e-02 5.88e-03 1.48e-03 3.72e-04 | ord  1.97  1.99  1.99
pulse_k12800 hrx0 be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrx0 be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrx0 be       C e      lev 32-64-128-256-512  L1 1.54e-02 4.48e-03 1.28e-03 4.09e-04 | ord  1.79  1.80  1.65
pulse_k12800 hrx0 be       C f1     lev 32-64-128-256-512  L1 2.19e-02 6.59e-03 2.16e-03 8.28e-04 | ord  1.73  1.61  1.38
pulse_k12800 hrx0 hesdirk2 X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrx0 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.29e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrx0 hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.17e-05 2.95e-06 7.45e-07 | ord  1.98  1.99  1.99
pulse_k12800 hrx0 hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.21e-05 1.06e-05 2.68e-06 | ord  1.98  1.99  1.99
pulse_k12800 hrx0 hesdirk2 C e      lev 32-64-128-256-512  L1 1.44e-02 3.80e-03 9.41e-04 2.37e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrx0 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.28e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
rw_t10       cen  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.23e-03 3.17e-04 | ord  2.17  2.04  1.96
rw_t10       cen  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 1.10e-03 | ord  0.44  0.12  6.95
rw_t10       cen  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       cen  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 3.13e-04 | ord  2.17  2.15  2.13
rw_t10       cen  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       cen  be       T e      lev 1-2-4-8-16  L1 7.83e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       cen  be       T f1     lev 1-2-4-8-16  L1 6.62e-02 3.55e-02 1.85e-02 9.38e-03 | ord  0.90  0.94  0.98
rw_t10       cen  be       T dens   lev 1-2-4-8-16  L1 6.75e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       cen  be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.00e-04 | ord  1.00  1.00  1.00
rw_t10       cen  be       T eint   lev 1-2-4-8-16  L1 6.78e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       cen  be       C e      lev 32-64-128-256-512  L1 2.00e-01 1.33e-01 7.95e-02 4.36e-02 | ord  0.59  0.74  0.87
rw_t10       cen  be       C f1     lev 32-64-128-256-512  L1 2.11e-01 1.22e-01 6.67e-02 3.63e-02 | ord  0.79  0.87  0.88
rw_t10       cen  be       C dens   lev 32-64-128-256-512  L1 4.40e-02 1.51e-02 6.85e-03 3.40e-03 | ord  1.54  1.14  1.01
rw_t10       cen  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.48e-03 3.22e-03 | ord  1.55  1.08  1.01
rw_t10       cen  be       C eint   lev 32-64-128-256-512  L1 4.45e-02 1.54e-02 6.91e-03 3.42e-03 | ord  1.54  1.15  1.01
rw_t10       cen  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 2.81e-04 | ord  1.91  1.97  1.90
rw_t10       cen  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.67e-01 4.76e-02 4.50e-02 1.10e-03 | ord  1.81  0.08  5.36
rw_t10       cen  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.97e-04 | ord  2.18  2.13  2.19
rw_t10       cen  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.16e-04 | ord  2.16  2.14  2.14
rw_t10       cen  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.97e-04 | ord  2.18  2.13  2.19
rw_t10       cen  hesdirk2 T e      lev 1-2-4-8-16  L1 8.57e-04 4.25e-04 2.12e-04 1.06e-04 | ord  1.01  1.00  1.00
rw_t10       cen  hesdirk2 T f1     lev 1-2-4-8-16  L1 7.30e-04 2.67e-04 6.35e-04 8.44e-04 | ord  1.45 -1.25 -0.41
rw_t10       cen  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.33e-04 3.16e-04 1.58e-04 | ord  1.00  1.00  1.00
rw_t10       cen  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.38e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.00
rw_t10       cen  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.00  1.00  1.00
rw_t10       cen  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.26e-04 2.42e-04 | ord  2.35  2.42  1.12
rw_t10       cen  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 3.81e-03 1.32e-03 | ord  2.49  2.88  1.52
rw_t10       cen  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.19e-04 | ord  2.12  1.96  1.46
rw_t10       cen  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.17e-04 | ord  2.16  2.07  1.31
rw_t10       cen  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.19e-04 | ord  2.12  1.96  1.46
rw_t10       hr   be       X e      lev 32-64-128-256-512  L1 2.28e-02 3.93e-02 1.14e-01 1.96e-01 | ord -0.79 -1.53 -0.79
rw_t10       hr   be       X f1     lev 32-64-128-256-512  L1 2.06e-01 2.39e-01 5.75e-01 1.02e+00 | ord -0.22 -1.27 -0.83
rw_t10       hr   be       X dens   lev 32-64-128-256-512  L1 2.66e-02 6.24e-03 2.68e-03 2.73e-03 | ord  2.09  1.22 -0.03
rw_t10       hr   be       X velx   lev 32-64-128-256-512  L1 2.73e-02 6.44e-03 3.36e-03 4.84e-03 | ord  2.08  0.94 -0.53
rw_t10       hr   be       X eint   lev 32-64-128-256-512  L1 2.67e-02 6.30e-03 2.86e-03 2.99e-03 | ord  2.08  1.14 -0.06
rw_t10       hr   be       T e      lev 1-2-4-8-16  L1 8.43e-02 3.58e-02 1.34e-02 2.90e-03 | ord  1.23  1.42  2.21
rw_t10       hr   be       T f1     lev 1-2-4-8-16  L1 6.77e-02 5.03e-02 6.64e-02 7.34e-02 | ord  0.43 -0.40 -0.14
rw_t10       hr   be       T dens   lev 1-2-4-8-16  L1 5.03e-03 2.54e-03 1.30e-03 6.84e-04 | ord  0.99  0.96  0.93
rw_t10       hr   be       T velx   lev 1-2-4-8-16  L1 5.22e-03 2.64e-03 1.36e-03 7.36e-04 | ord  0.98  0.96  0.88
rw_t10       hr   be       T eint   lev 1-2-4-8-16  L1 5.05e-03 2.55e-03 1.30e-03 6.84e-04 | ord  0.99  0.97  0.93
rw_t10       hr   be       C e      lev 32-64-128-256-512  L1 2.63e-01 2.63e-01 2.12e-01 3.34e-02 | ord  0.00  0.31  2.67
rw_t10       hr   be       C f1     lev 32-64-128-256-512  L1 2.58e-01 2.85e-01 7.16e-01 1.06e+00 | ord -0.14 -1.33 -0.57
rw_t10       hr   be       C dens   lev 32-64-128-256-512  L1 4.52e-02 1.71e-02 8.08e-03 2.60e-03 | ord  1.40  1.08  1.64
rw_t10       hr   be       C velx   lev 32-64-128-256-512  L1 4.12e-02 1.52e-02 7.71e-03 4.32e-03 | ord  1.44  0.98  0.84
rw_t10       hr   be       C eint   lev 32-64-128-256-512  L1 4.58e-02 1.75e-02 8.23e-03 2.60e-03 | ord  1.39  1.09  1.66
rw_t10       hr   hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 5.30e-03 1.69e-02 3.83e-02 | ord  1.54 -1.67 -1.18
rw_t10       hr   hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 1.53e-01 4.42e-01 8.66e-01 | ord  0.14 -1.53 -0.97
rw_t10       hr   hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.86e-03 1.22e-03 4.96e-04 | ord  2.19  2.27  1.30
rw_t10       hr   hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.11e-03 2.39e-03 4.60e-03 | ord  2.16  1.36 -0.95
rw_t10       hr   hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.86e-03 1.20e-03 5.21e-04 | ord  2.19  2.29  1.20
rw_t10       hr   hesdirk2 T e      lev 1-2-4-8-16  L1 6.87e-04 1.25e-03 1.31e-03 1.04e-03 | ord -0.86 -0.07  0.34
rw_t10       hr   hesdirk2 T f1     lev 1-2-4-8-16  L1 6.15e-02 7.27e-02 6.83e-02 5.18e-02 | ord -0.24  0.09  0.40
rw_t10       hr   hesdirk2 T dens   lev 1-2-4-8-16  L1 1.28e-03 6.48e-04 3.31e-04 1.69e-04 | ord  0.99  0.97  0.97
rw_t10       hr   hesdirk2 T velx   lev 1-2-4-8-16  L1 1.31e-03 6.93e-04 4.10e-04 2.71e-04 | ord  0.92  0.76  0.60
rw_t10       hr   hesdirk2 T eint   lev 1-2-4-8-16  L1 1.28e-03 6.48e-04 3.31e-04 1.69e-04 | ord  0.99  0.97  0.97
rw_t10       hr   hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 8.81e-03 2.06e-02 3.22e-02 | ord  0.73 -1.22 -0.65
rw_t10       hr   hesdirk2 C f1     lev 32-64-128-256-512  L1 1.79e-01 2.28e-01 5.01e-01 6.77e-01 | ord -0.35 -1.14 -0.43
rw_t10       hr   hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.53e-03 1.43e-03 2.77e-04 | ord  2.13  2.19  2.36
rw_t10       hr   hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.40e-03 2.78e-03 3.61e-03 | ord  2.16  1.21 -0.38
rw_t10       hr   hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.52e-03 1.40e-03 3.08e-04 | ord  2.13  2.22  2.18
rw_t10       hrb  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 3.83e-04 | ord  2.17  2.02  1.70
rw_t10       hrb  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 3.33e-03 | ord  0.44  0.12  5.36
rw_t10       hrb  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       hrb  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 3.13e-04 | ord  2.17  2.15  2.13
rw_t10       hrb  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       hrb  be       T e      lev 1-2-4-8-16  L1 7.84e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       hrb  be       T f1     lev 1-2-4-8-16  L1 6.67e-02 3.57e-02 1.86e-02 9.29e-03 | ord  0.90  0.94  1.00
rw_t10       hrb  be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       hrb  be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       hrb  be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       hrb  be       C e      lev 32-64-128-256-512  L1 2.04e-01 1.38e-01 7.95e-02 4.36e-02 | ord  0.56  0.79  0.87
rw_t10       hrb  be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.74e-02 3.64e-02 | ord  0.76  0.88  0.89
rw_t10       hrb  be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.86e-03 3.40e-03 | ord  1.53  1.15  1.01
rw_t10       hrb  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.47e-03 3.22e-03 | ord  1.55  1.08  1.01
rw_t10       hrb  be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.92e-03 3.42e-03 | ord  1.52  1.16  1.02
rw_t10       hrb  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 3.24e-04 | ord  1.91  1.97  1.69
rw_t10       hrb  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 1.75e-03 | ord  1.82  0.08  4.69
rw_t10       hrb  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.20
rw_t10       hrb  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.13e-04 | ord  2.16  2.14  2.15
rw_t10       hrb  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.94e-04 | ord  2.18  2.13  2.20
rw_t10       hrb  hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.29e-04 2.16e-04 1.06e-04 | ord  1.00  0.99  1.03
rw_t10       hrb  hesdirk2 T f1     lev 1-2-4-8-16  L1 4.62e-04 5.83e-04 2.01e-04 1.23e-04 | ord -0.33  1.54  0.70
rw_t10       hrb  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrb  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       hrb  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrb  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.27e-04 2.49e-04 | ord  2.35  2.41  1.08
rw_t10       hrb  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 3.77e-03 1.90e-03 | ord  2.48  2.90  0.99
rw_t10       hrb  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       hrb  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.17e-04 | ord  2.16  2.07  1.31
rw_t10       hrb  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       hrx0 be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.25e-03 3.83e-04 | ord  2.17  2.02  1.70
rw_t10       hrx0 be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 3.33e-03 | ord  0.44  0.12  5.36
rw_t10       hrx0 be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       hrx0 be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 3.13e-04 | ord  2.17  2.15  2.13
rw_t10       hrx0 be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 2.98e-04 | ord  2.18  2.14  2.16
rw_t10       hrx0 be       T e      lev 1-2-4-8-16  L1 7.84e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       hrx0 be       T f1     lev 1-2-4-8-16  L1 6.67e-02 3.57e-02 1.86e-02 9.29e-03 | ord  0.90  0.94  1.00
rw_t10       hrx0 be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       hrx0 be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       hrx0 be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       hrx0 be       C e      lev 32-64-128-256-512  L1 2.04e-01 1.38e-01 7.95e-02 4.36e-02 | ord  0.56  0.79  0.87
rw_t10       hrx0 be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.74e-02 3.64e-02 | ord  0.76  0.88  0.89
rw_t10       hrx0 be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.86e-03 3.40e-03 | ord  1.53  1.15  1.01
rw_t10       hrx0 be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.47e-03 3.22e-03 | ord  1.55  1.08  1.01
rw_t10       hrx0 be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.92e-03 3.42e-03 | ord  1.52  1.16  1.02
rw_t10       hrx0 hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 3.24e-04 | ord  1.91  1.97  1.69
rw_t10       hrx0 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 1.75e-03 | ord  1.82  0.08  4.69
rw_t10       hrx0 hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.20
rw_t10       hrx0 hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.13e-04 | ord  2.16  2.14  2.15
rw_t10       hrx0 hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.94e-04 | ord  2.18  2.13  2.20
rw_t10       hrx0 hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.29e-04 2.16e-04 1.06e-04 | ord  1.00  0.99  1.03
rw_t10       hrx0 hesdirk2 T f1     lev 1-2-4-8-16  L1 4.62e-04 5.83e-04 2.01e-04 1.23e-04 | ord -0.33  1.54  0.70
rw_t10       hrx0 hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrx0 hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       hrx0 hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrx0 hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.27e-04 2.49e-04 | ord  2.35  2.41  1.08
rw_t10       hrx0 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 3.77e-03 1.90e-03 | ord  2.48  2.90  0.99
rw_t10       hrx0 hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t10       hrx0 hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.17e-04 | ord  2.16  2.07  1.31
rw_t10       hrx0 hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.14e-04 | ord  2.12  1.96  1.47
rw_t1000     cen  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     cen  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     cen  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     cen  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     cen  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     cen  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     cen  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     cen  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     cen  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     cen  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     cen  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     cen  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     cen  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     cen  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     cen  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     cen  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     cen  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     cen  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     cen  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     cen  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     cen  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     cen  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     cen  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     cen  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     cen  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     cen  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     cen  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     cen  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     cen  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     cen  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     hr   be       X e      lev 32-64-128-256  L1 8.56e-03 2.20e-03 8.67e-04 | ord  1.96  1.34
rw_t1000     hr   be       X f1     lev 32-64-128-256  L1 4.89e-03 1.49e-02 9.42e-04 | ord -1.61  3.98
rw_t1000     hr   be       X dens   lev 32-64-128-256  L1 8.81e-03 2.20e-03 8.60e-04 | ord  2.00  1.36
rw_t1000     hr   be       X velx   lev 32-64-128-256  L1 8.48e-03 2.20e-03 8.68e-04 | ord  1.95  1.34
rw_t1000     hr   be       X eint   lev 32-64-128-256  L1 8.78e-03 2.19e-03 8.59e-04 | ord  2.00  1.35
rw_t1000     hr   be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.46e-04 | ord  0.97  0.99
rw_t1000     hr   be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.48e-04 | ord  4.03  0.99
rw_t1000     hr   be       T dens   lev 2-4-8-16  L1 3.27e-03 1.66e-03 8.37e-04 | ord  0.98  0.99
rw_t1000     hr   be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.47e-04 | ord  0.97  0.99
rw_t1000     hr   be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.40e-04 | ord  0.97  0.99
rw_t1000     hr   be       C e      lev 32-64-128-256  L1 2.49e-02 1.16e-02 5.91e-03 | ord  1.10  0.98
rw_t1000     hr   be       C f1     lev 32-64-128-256  L1 2.63e-02 1.33e-02 3.03e-02 | ord  0.98 -1.19
rw_t1000     hr   be       C dens   lev 32-64-128-256  L1 2.63e-02 1.19e-02 5.95e-03 | ord  1.14  1.00
rw_t1000     hr   be       C velx   lev 32-64-128-256  L1 2.52e-02 1.19e-02 6.00e-03 | ord  1.08  0.99
rw_t1000     hr   be       C eint   lev 32-64-128-256  L1 2.57e-02 1.18e-02 5.94e-03 | ord  1.12  0.99
rw_t1000     hr   hesdirk2 X e      lev 32-64-128-256  L1 8.56e-03 2.20e-03 8.62e-04 | ord  1.96  1.35
rw_t1000     hr   hesdirk2 X f1     lev 32-64-128-256  L1 5.84e-03 6.04e-03 9.37e-04 | ord -0.05  2.69
rw_t1000     hr   hesdirk2 X dens   lev 32-64-128-256  L1 8.84e-03 2.20e-03 8.51e-04 | ord  2.01  1.37
rw_t1000     hr   hesdirk2 X velx   lev 32-64-128-256  L1 8.52e-03 2.20e-03 8.64e-04 | ord  1.95  1.35
rw_t1000     hr   hesdirk2 X eint   lev 32-64-128-256  L1 8.80e-03 2.17e-03 8.52e-04 | ord  2.02  1.35
rw_t1000     hr   hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     hr   hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.78e-03 | ord  1.12  0.99
rw_t1000     hr   hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     hr   hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     hr   hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.01
rw_t1000     hr   hesdirk2 C e      lev 32-64-128-256  L1 8.00e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     hr   hesdirk2 C f1     lev 32-64-128-256  L1 6.33e-02 3.35e-02 1.79e-02 | ord  0.92  0.90
rw_t1000     hr   hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.86e-02 1.89e-02 | ord  1.04  1.03
rw_t1000     hr   hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.73e-02 | ord  0.94  0.96
rw_t1000     hr   hesdirk2 C eint   lev 32-64-128-256  L1 7.93e-02 3.82e-02 1.86e-02 | ord  1.05  1.04
rw_t1000     hrb  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrb  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     hrb  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     hrb  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     hrb  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     hrb  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     hrb  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     hrb  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     hrb  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     hrb  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     hrb  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     hrb  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     hrb  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     hrb  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     hrb  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     hrb  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrb  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     hrb  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     hrb  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     hrb  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     hrb  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     hrb  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     hrb  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     hrb  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     hrb  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     hrb  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     hrb  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     hrb  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     hrb  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     hrb  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     hrx0 be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrx0 be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     hrx0 be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     hrx0 be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     hrx0 be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     hrx0 be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     hrx0 be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     hrx0 be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     hrx0 be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     hrx0 be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     hrx0 be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     hrx0 be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     hrx0 be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     hrx0 be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     hrx0 be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     hrx0 hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrx0 hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     hrx0 hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     hrx0 hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     hrx0 hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     hrx0 hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     hrx0 hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     hrx0 hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     hrx0 hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     hrx0 hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     hrx0 hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     hrx0 hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     hrx0 hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     hrx0 hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     hrx0 hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
```

## RESULTS/atm_hopf.txt (raw)

```
cen n=  32  L1 1.398e-02  max|dE/E| 1.633e-01  max(tau<1) 1.633e-01  max(tau>1) 5.547e-02  top 1.484e-01  tau_cell top 4.0e-03 bot 2.8e+01
cen n=  64  L1 9.676e-03  max|dE/E| 1.675e-01  max(tau<1) 1.675e-01  max(tau>1) 5.224e-02  top 1.516e-01  tau_cell top 2.0e-03 bot 1.4e+01
cen n= 128  L1 8.664e-03  max|dE/E| 1.697e-01  max(tau<1) 1.697e-01  max(tau>1) 5.578e-02  top 1.531e-01  tau_cell top 1.0e-03 bot 7.0e+00
cen n= 256  L1 8.443e-03  max|dE/E| 1.708e-01  max(tau<1) 1.708e-01  max(tau>1) 5.544e-02  top 1.539e-01  tau_cell top 5.0e-04 bot 3.5e+00
cen n= 512  L1 8.396e-03  max|dE/E| 1.713e-01  max(tau<1) 1.713e-01  max(tau>1) 5.590e-02  top 1.543e-01  tau_cell top 2.5e-04 bot 1.7e+00
cen Hopf L1 orders 0.53 0.16 0.04 0.01
cen self m1_e 1.82e-03 4.59e-04 1.18e-04 3.86e-05 | ord 1.99 1.96 1.61
cen self m1_f1 9.17e-01 7.37e-01 1.12e+00 5.44e+00 | ord 0.32 -0.60 -2.28
hr n=  32  L1 5.472e-03  max|dE/E| 1.788e-02  max(tau<1) 1.788e-02  max(tau>1) 1.206e-02  top 8.650e-04  tau_cell top 4.0e-03 bot 2.8e+01
hr n=  64  L1 1.498e-03  max|dE/E| 1.722e-02  max(tau<1) 1.722e-02  max(tau>1) 1.325e-02  top 6.658e-04  tau_cell top 2.0e-03 bot 1.4e+01
hr n= 128  L1 1.946e-03  max|dE/E| 1.774e-02  max(tau<1) 1.774e-02  max(tau>1) 1.468e-02  top 8.608e-04  tau_cell top 1.0e-03 bot 7.0e+00
hr n= 256  L1 3.054e-03  max|dE/E| 1.956e-02  max(tau<1) 1.956e-02  max(tau>1) 1.303e-02  top 9.828e-04  tau_cell top 5.0e-04 bot 3.5e+00
hr n= 512  L1 4.236e-03  max|dE/E| 2.015e-02  max(tau<1) 2.015e-02  max(tau>1) 9.668e-03  top 8.396e-04  tau_cell top 2.5e-04 bot 1.7e+00
hr Hopf L1 orders 1.87 -0.38 -0.65 -0.47
hr self m1_e 2.11e-03 9.83e-04 8.98e-04 1.10e-03 | ord 1.10 0.13 -0.29
hr self m1_f1 1.27e+00 1.91e+01 5.60e+01 9.75e-01 | ord -3.91 -1.55 5.85
hrb n=  32  L1 5.775e-03  max|dE/E| 1.824e-02  max(tau<1) 1.824e-02  max(tau>1) 1.075e-02  top 6.644e-04  tau_cell top 4.0e-03 bot 2.8e+01
hrb n=  64  L1 1.654e-03  max|dE/E| 1.740e-02  max(tau<1) 1.740e-02  max(tau>1) 1.095e-02  top 3.152e-04  tau_cell top 2.0e-03 bot 1.4e+01
hrb n= 128  L1 8.282e-04  max|dE/E| 1.689e-02  max(tau<1) 1.689e-02  max(tau>1) 1.152e-02  top 2.698e-04  tau_cell top 1.0e-03 bot 7.0e+00
hrb n= 256  L1 6.970e-04  max|dE/E| 1.843e-02  max(tau<1) 1.843e-02  max(tau>1) 9.568e-03  top 2.108e-04  tau_cell top 5.0e-04 bot 3.5e+00
hrb n= 512  L1 6.082e-04  max|dE/E| 1.916e-02  max(tau<1) 1.916e-02  max(tau>1) 7.854e-03  top 1.412e-04  tau_cell top 2.5e-04 bot 1.7e+00
hrb Hopf L1 orders 1.80 1.00 0.25 0.20
hrb self m1_e 1.81e-03 4.67e-04 1.65e-04 1.30e-04 | ord 1.96 1.50 0.34
hrb self m1_f1 1.27e+00 1.35e+01 3.12e+01 1.22e+00 | ord -3.41 -1.20 4.67
hrx0 n=  32  L1 5.936e-03  max|dE/E| 1.536e-02  max(tau<1) 1.536e-02  max(tau>1) 6.389e-03  top 8.197e-04  tau_cell top 4.0e-03 bot 2.8e+01
hrx0 n=  64  L1 1.684e-03  max|dE/E| 1.732e-02  max(tau<1) 1.732e-02  max(tau>1) 5.168e-03  top 2.821e-04  tau_cell top 2.0e-03 bot 1.4e+01
hrx0 n= 128  L1 6.658e-04  max|dE/E| 1.840e-02  max(tau<1) 1.840e-02  max(tau>1) 5.471e-03  top 8.283e-05  tau_cell top 1.0e-03 bot 7.0e+00
hrx0 n= 256  L1 4.304e-04  max|dE/E| 1.893e-02  max(tau<1) 1.893e-02  max(tau>1) 4.958e-03  top 1.520e-05  tau_cell top 5.0e-04 bot 3.5e+00
hrx0 n= 512  L1 3.714e-04  max|dE/E| 1.928e-02  max(tau<1) 1.928e-02  max(tau>1) 4.742e-03  top 3.278e-05  tau_cell top 2.5e-04 bot 1.7e+00
hrx0 Hopf L1 orders 1.82 1.34 0.63 0.21
hrx0 self m1_e 1.77e-03 4.41e-04 1.11e-04 3.91e-05 | ord 2.00 1.99 1.51
hrx0 self m1_f1 1.62e+00 1.02e+01 3.29e+01 9.05e-01 | ord -2.65 -1.69 5.18
```

## RESULTS/arm_minus_cen.txt (raw)

```
pulse_k0.128  hr   hesdirk2 X |hr-cen| 1.43e-05 6.13e-05 2.56e-04 1.05e-03 4.23e-03
pulse_k0.128  hr   hesdirk2 C |hr-cen| 1.05e-04 2.29e-04 4.80e-04 9.82e-04 1.99e-03
pulse_k0.128  hr   hesdirk2 T |hr-cen| 7.83e-03 3.92e-03 1.96e-03 9.82e-04 4.91e-04
pulse_k0.128  hr   be       X |hr-cen| 4.92e-05 2.10e-04 8.77e-04 3.58e-03 1.44e-02
pulse_k0.128  hr   be       C |hr-cen| 3.73e-04 7.89e-04 1.65e-03 3.35e-03 6.76e-03
pulse_k0.128  hr   be       T |hr-cen| 2.67e-02 1.34e-02 6.70e-03 3.35e-03 1.68e-03
pulse_k0.128  hrb  hesdirk2 X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrb  hesdirk2 C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrb  hesdirk2 T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrb  be       X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrb  be       C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrb  be       T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hr   hesdirk2 X |hr-cen| 3.38e-06 1.21e-05 4.96e-05 2.02e-04 8.20e-04
pulse_k12.8   hr   hesdirk2 C |hr-cen| 2.59e-05 4.55e-05 9.23e-05 1.90e-04 3.87e-04
pulse_k12.8   hr   hesdirk2 T |hr-cen| 1.42e-03 7.37e-04 3.76e-04 1.90e-04 9.55e-05
pulse_k12.8   hr   be       X |hr-cen| 1.12e-05 4.03e-05 1.66e-04 6.78e-04 2.74e-03
pulse_k12.8   hr   be       C |hr-cen| 7.07e-05 1.40e-04 3.02e-04 6.36e-04 1.31e-03
pulse_k12.8   hr   be       T |hr-cen| 4.17e-03 2.34e-03 1.24e-03 6.36e-04 3.23e-04
pulse_k12.8   hrb  hesdirk2 X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrb  hesdirk2 C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrb  hesdirk2 T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrb  be       X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrb  be       C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrb  be       T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hr   hesdirk2 X |hr-cen| 2.66e-06 5.35e-06 1.07e-05 2.13e-05 4.21e-05
pulse_k128    hr   hesdirk2 C |hr-cen| 1.51e-05 1.78e-05 1.94e-05 2.00e-05 1.97e-05
pulse_k128    hr   hesdirk2 T |hr-cen| 1.27e-04 7.26e-05 3.89e-05 2.00e-05 1.01e-05
pulse_k128    hr   be       X |hr-cen| 7.88e-06 1.58e-05 3.15e-05 6.20e-05 1.19e-04
pulse_k128    hr   be       C |hr-cen| 2.80e-05 4.04e-05 5.14e-05 5.87e-05 6.11e-05
pulse_k128    hr   be       T |hr-cen| 2.28e-04 1.61e-04 1.02e-04 5.87e-05 3.18e-05
pulse_k128    hrb  hesdirk2 X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrb  hesdirk2 C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrb  hesdirk2 T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrb  be       X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrb  be       C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrb  be       T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hr   hesdirk2 X |hr-cen| 4.12e-06 8.16e-06 1.62e-05 3.24e-05 6.49e-05
pulse_k1280   hr   hesdirk2 C |hr-cen| 4.71e-06 9.11e-06 1.74e-05 3.21e-05 5.58e-05
pulse_k1280   hr   hesdirk2 T |hr-cen| 3.72e-05 3.63e-05 3.48e-05 3.21e-05 2.79e-05
pulse_k1280   hr   be       X |hr-cen| 4.58e-06 9.06e-06 1.80e-05 3.60e-05 7.21e-05
pulse_k1280   hr   be       C |hr-cen| 4.74e-06 9.35e-06 1.84e-05 3.59e-05 6.85e-05
pulse_k1280   hr   be       T |hr-cen| 3.74e-05 3.73e-05 3.68e-05 3.59e-05 3.42e-05
pulse_k1280   hrb  hesdirk2 X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrb  hesdirk2 C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrb  hesdirk2 T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrb  be       X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrb  be       C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrb  be       T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hr   hesdirk2 X |hr-cen| 4.82e-07 9.52e-07 1.89e-06 3.78e-06 7.56e-06
pulse_k12800  hr   hesdirk2 C |hr-cen| 4.83e-07 9.53e-07 1.89e-06 3.78e-06 7.55e-06
pulse_k12800  hr   hesdirk2 T |hr-cen| 3.79e-06 3.79e-06 3.78e-06 3.78e-06 3.77e-06
pulse_k12800  hr   be       X |hr-cen| 4.82e-07 9.52e-07 1.89e-06 3.78e-06 7.56e-06
pulse_k12800  hr   be       C |hr-cen| 4.79e-07 9.49e-07 1.89e-06 3.78e-06 7.56e-06
pulse_k12800  hr   be       T |hr-cen| 3.76e-06 3.77e-06 3.78e-06 3.78e-06 3.78e-06
pulse_k12800  hrb  hesdirk2 X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrb  hesdirk2 C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrb  hesdirk2 T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrb  be       X |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrb  be       C |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrb  be       T |hrb-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
rw_t10        hr   hesdirk2 X |hr-cen| 8.09e-06 1.05e-06 5.90e-05 3.03e-04 8.68e-04
rw_t10        hr   hesdirk2 C |hr-cen| 2.53e-05 1.36e-05 9.20e-05 3.83e-04 8.61e-04
rw_t10        hr   hesdirk2 T |hr-cen| 9.18e-05 7.76e-05 5.76e-05 3.73e-05 2.19e-05
rw_t10        hr   be       X |hr-cen| 8.36e-05 2.74e-04 9.04e-04 2.59e-03 5.17e-03
rw_t10        hr   be       C |hr-cen| 1.51e-03 2.68e-03 4.67e-03 6.04e-03 4.98e-03
rw_t10        hr   be       T |hr-cen| 4.58e-03 2.06e-03 8.48e-04 2.98e-04 7.85e-05
rw_t10        hrb  hesdirk2 X |hrb-cen| 7.33e-07 6.13e-07 4.32e-07 3.22e-07 4.91e-06
rw_t10        hrb  hesdirk2 C |hrb-cen| 2.73e-07 8.34e-08 7.10e-08 9.16e-08 4.65e-06
rw_t10        hrb  hesdirk2 T |hrb-cen| 4.03e-08 2.70e-08 3.72e-07 8.19e-08 1.39e-08
rw_t10        hrb  be       X |hrb-cen| 3.23e-08 2.88e-08 3.19e-08 1.08e-06 6.51e-07
rw_t10        hrb  be       C |hrb-cen| 2.76e-08 1.39e-04 1.67e-05 9.26e-07 4.48e-07
rw_t10        hrb  be       T |hrb-cen| 1.63e-05 4.53e-06 2.90e-08 3.53e-07 1.34e-07
rw_t1000      hr   hesdirk2 X |hr-cen| 8.19e-05 1.64e-04 3.28e-04 6.56e-04
rw_t1000      hr   hesdirk2 C |hr-cen| 8.09e-05 1.65e-04 3.32e-04 6.60e-04
rw_t1000      hr   hesdirk2 T |hr-cen| 3.30e-04 3.28e-04 3.22e-04 3.10e-04
rw_t1000      hr   be       X |hr-cen| 8.33e-05 1.67e-04 3.33e-04 6.66e-04
rw_t1000      hr   be       C |hr-cen| 8.41e-05 1.68e-04 3.34e-04 6.68e-04
rw_t1000      hr   be       T |hr-cen| 3.34e-04 3.33e-04 3.31e-04 3.27e-04
rw_t1000      hrb  hesdirk2 X |hrb-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      hrb  hesdirk2 C |hrb-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      hrb  hesdirk2 T |hrb-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      hrb  be       X |hrb-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      hrb  be       C |hrb-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      hrb  be       T |hrb-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
```

## RESULTS/gates.txt (raw)

```
g1_hrb_1e4 rc=0 Picard iterations mean=2.164706e+01 | HOPF: max|dE/E| 1.689e-02 top5 4.783e-03 L1 8.282e-04 (E_top 1.7349 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hr_1e4 rc=0 Picard iterations mean=2.150980e+01 | HOPF: max|dE/E| 1.774e-02 top5 5.409e-03 L1 1.946e-03 (E_top 1.7359 exact 1.7344, tau_top-cell 5.09e-04) 
g5_hrb_k12.8 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.940470e-02 2D=5.208333e-02 ratio=0.564570 (target 1.000000 +- 0.02) 
g5_hrb_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.780463e-03 2D=5.208333e-03 ratio=0.533849 (target 1.000000 +- 0.02) 
g5_hr_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.756269e-03 2D=5.208333e-03 ratio=0.529204 (target 1.000000 +- 0.02) 
g5_hrb_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_hr_k12.8 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.949960e-02 2D=5.208333e-02 ratio=0.566392 (target 1.000000 +- 0.02) 
g5_hr_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208163e-05 2D=5.208333e-05 ratio=0.999967 (target 1.000000 +- 0.02) 
g5_hrb_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295065e-04 2D=5.208333e-04 ratio=1.016652 (target 1.000000 +- 0.02) 
g5_hr_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295020e-04 2D=5.208333e-04 ratio=1.016644 (target 1.000000 +- 0.02) 
g1_hrb_1e2 rc=0 Picard iterations mean=2.058801e+00 | HOPF: max|dE/E| 1.690e-02 top5 4.667e-03 L1 7.360e-04 (E_top 1.7347 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hr_1e2 rc=0 Picard iterations mean=2.060168e+00 | HOPF: max|dE/E| 1.751e-02 top5 5.118e-03 L1 1.557e-03 (E_top 1.7354 exact 1.7344, tau_top-cell 5.09e-04) 
g3_hrb_c10 rc=0 Picard iterations mean=2.523110e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/xthinfix_1009/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hr_c10 rc=0 Picard iterations mean=2.526977e+00 | PASS T6 marshak: L1(E)=0.0140 L1(material)=0.0140 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/xthinfix_1009/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hr_c1 rc=0 Picard iterations mean=2.006111e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/xthinfix_1009/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hrb_c1 rc=0 Picard iterations mean=2.006053e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/xthinfix_1009/repo/tests_m1/runs_3a/t6_ref_sn.txt 
```

## RESULTS/beams.txt (raw)

```
shd3b hr rc=0 fatal=0 318.78 s wall | Picard iterations mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=2.065333e+01 max=5.100000e+01 | breakdowns=0.000000e+00
shd3b hrb rc=0 fatal=0 883.46 s wall | Picard iterations mean=2.990000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.348094e+02 max=2.000000e+02 | breakdowns=0.000000e+00
shd3b cen rc=0 fatal=0 566.94 s wall | Picard iterations mean=3.030000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=7.533993e+01 max=2.000000e+02 | breakdowns=0.000000e+00
cyl hr rc=0 fatal=0 128.44 s wall | Picard iterations mean=5.000000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.678400e+01 max=1.620000e+02 | breakdowns=3.000000e+00
cyl hrb rc=0 fatal=0 224.74 s wall | Picard iterations mean=4.980000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=5.128715e+01 max=2.000000e+02 | breakdowns=6.000000e+00
cyl cen rc=0 fatal=0 129.01 s wall | Picard iterations mean=2.180000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=4.687156e+01 max=2.000000e+02 | breakdowns=2.000000e+00
xb20 hr rc=0 fatal=0 331.75 s wall | Picard iterations mean=4.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=3.731500e+01 max=8.000000e+01 | breakdowns=0.000000e+00
xb20 hrb rc=0 fatal=0 2594.32 s wall | Picard iterations mean=1.194000e+01 max=2.300000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.719179e+02 max=2.000000e+02 | breakdowns=0.000000e+00
xb20 cen rc=0 fatal=0 501.12 s wall | Picard iterations mean=3.700000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=9.399189e+01 max=2.000000e+02 | breakdowns=0.000000e+00
ba0 hr rc=0 fatal=0 388.89 s wall | Picard iterations mean=4.750000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=3.984421e+01 max=1.000000e+02 | breakdowns=0.000000e+00
ba0 hrb rc=0 fatal=0 1948.34 s wall | Picard iterations mean=8.590000e+00 max=1.800000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.749232e+02 max=2.000000e+02 | breakdowns=0.000000e+00
ba0 cen rc=0 fatal=0 684.99 s wall | Picard iterations mean=4.750000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.061305e+02 max=2.000000e+02 | breakdowns=0.000000e+00
ba20 hr rc=0 fatal=0 443.60 s wall | Picard iterations mean=5.440000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=4.164154e+01 max=1.400000e+02 | breakdowns=0.000000e+00
ba20 hrb rc=0 fatal=0 2307.02 s wall | Picard iterations mean=9.950000e+00 max=1.900000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.781678e+02 max=2.000000e+02 | breakdowns=2.000000e+00
ba20 cen rc=0 fatal=0 492.63 s wall | Picard iterations mean=4.200000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=7.870238e+01 max=2.000000e+02 | breakdowns=0.000000e+00
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/cyl/hr t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.377  max|E/J-1| 0.571  <|F|/|cH|> 0.673  median dangle 4.6 deg  p90 11.0  <f_M1> 0.502 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.673  median dangle 4.6  face max |F|/cE 0.65  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/cyl/hrb t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.383  max|E/J-1| 0.575  <|F|/|cH|> 0.667  median dangle 4.6 deg  p90 10.9  <f_M1> 0.502 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.667  median dangle 4.6  face max |F|/cE 0.65  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/cyl/cen t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.414  max|E/J-1| 0.635  <|F|/|cH|> 22.300  median dangle 87.8 deg  p90 143.2  <f_M1> 18.064 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 1.262  median dangle 87.8  face max |F|/cE 130.06  frac(|F_face| > cE) 0.994
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/xb20/hr t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 2.930  sum E/sum J 2.821  |F|/|cH| median 2.508  flux-weighted <|F|>/<|cH|> 2.426
  beam: direction error J-weighted mean 4.8 deg  median 3.8  p90 12.6   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.94
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.1008  max E/Jmax 1.0909
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.811
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.158  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.396
    peaks y: E [0.965 1.035]  J [1.004]   E max/J max 2.047  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.763
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.177  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.365
    peaks y: E [0.699 0.879 1.004 1.121 1.301]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.698  width(>50% max) E 0.391 J 0.281
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/xb20/hrb t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 2.848  sum E/sum J 2.736  |F|/|cH| median 2.454  flux-weighted <|F|>/<|cH|> 2.354
  beam: direction error J-weighted mean 4.8 deg  median 3.8  p90 12.6   face |F|/cE > 1: frac 0.0000 (beam) 0.2775 (all)  max 971.70
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.1714  max E/Jmax 0.9299
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.004
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 5.037  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.505
    peaks y: E [0.965 1.035]  J [1.004]   E max/J max 2.002  width(>50% max) E 0.172 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.740
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 3.883  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.348
    peaks y: E [0.707 0.887 1.004 1.113 1.293]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.618  width(>50% max) E 0.375 J 0.281
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/xb20/cen t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 1.788  sum E/sum J 1.802  |F|/|cH| median 884.853  flux-weighted <|F|>/<|cH|> 1178.540
  beam: direction error J-weighted mean 64.5 deg  median 64.1  p90 160.4   face |F|/cE > 1: frac 1.0000 (beam) 0.9999 (all)  max 3656.39
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.2839  max E/Jmax 0.7855
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.378
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 5.514  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.693
    peaks y: E [0.074 0.965 1.035 1.926]  J [1.004]   E max/J max 1.305  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.966
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 3.355  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.449
    peaks y: E [1.004]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.652  width(>50% max) E 0.281 J 0.281
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/ba0/hr t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.073  sum E/sum J 5.599  |F|/|cH| median 7.574  flux-weighted <|F|>/<|cH|> 5.862
  beam: direction error J-weighted mean 5.1 deg  median 1.0  p90 4.9   face |F|/cE > 1: frac 0.0096 (beam) 0.0018 (all)  max 78.88
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.1141  max E/Jmax 6.6123
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.942
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 11.140  E-centroid 0.9806 (exact 0.9732)  E outside the J>5% band / total 0.452
    peaks y: E [0.965]  J [0.980 0.996]   E max/J max 5.894  width(>50% max) E 0.070 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.040
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 12.732  E-centroid 0.9883 (exact 0.9633)  E outside the J>5% band / total 0.516
    peaks y: E [0.957]  J [0.988]   E max/J max 5.715  width(>50% max) E 0.078 J 0.062
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/ba0/hrb t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.056  sum E/sum J 5.539  |F|/|cH| median 7.518  flux-weighted <|F|>/<|cH|> 5.804
  beam: direction error J-weighted mean 5.2 deg  median 1.0  p90 4.9   face |F|/cE > 1: frac 0.0096 (beam) 0.2202 (all)  max 792.64
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.1945  max E/Jmax 25.5545
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.062
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 12.872  E-centroid 0.9768 (exact 0.9732)  E outside the J>5% band / total 0.525
    peaks y: E [0.965]  J [0.980 0.996]   E max/J max 5.882  width(>50% max) E 0.062 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.185
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 14.980  E-centroid 0.9448 (exact 0.9633)  E outside the J>5% band / total 0.594
    peaks y: E [0.957]  J [0.988]   E max/J max 5.669  width(>50% max) E 0.078 J 0.062
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/ba0/cen t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.967  sum E/sum J 6.898  |F|/|cH| median 2825.624  flux-weighted <|F|>/<|cH|> 3257.489
  beam: direction error J-weighted mean 55.2 deg  median 85.9  p90 169.3   face |F|/cE > 1: frac 0.9521 (beam) 0.9890 (all)  max 9542.21
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.3118  max E/Jmax 5.9405
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.542
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 35.426  E-centroid 0.9595 (exact 0.9732)  E outside the J>5% band / total 0.778
    peaks y: E [0.973]  J [0.980 0.996]   E max/J max 7.487  width(>50% max) E 0.062 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.309
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 23.945  E-centroid 0.9979 (exact 0.9633)  E outside the J>5% band / total 0.662
    peaks y: E [0.965]  J [0.988]   E max/J max 7.892  width(>50% max) E 0.070 J 0.062
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/ba20/hr t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 2.425  sum E/sum J 2.335  |F|/|cH| median 2.145  flux-weighted <|F|>/<|cH|> 1.999
  beam: direction error J-weighted mean 12.2 deg  median 7.0  p90 31.4   face |F|/cE > 1: frac 0.0113 (beam) 0.0015 (all)  max 78.23
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0784  max E/Jmax 3.0565
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.024
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 6.069  E-centroid 0.8618 (exact 0.8959)  E outside the J>5% band / total 0.489
    peaks y: E [0.895]  J [0.918]   E max/J max 2.855  width(>50% max) E 0.078 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.256
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 7.469  E-centroid 0.9191 (exact 0.9590)  E outside the J>5% band / total 0.592
    peaks y: E [0.918]  J [0.980]   E max/J max 2.468  width(>50% max) E 0.156 J 0.070
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/ba20/hrb t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 1.902  sum E/sum J 1.830  |F|/|cH| median 1.670  flux-weighted <|F|>/<|cH|> 1.569
  beam: direction error J-weighted mean 12.1 deg  median 6.6  p90 32.1   face |F|/cE > 1: frac 0.0113 (beam) 0.2022 (all)  max 579.37
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.2151  max E/Jmax 10.1671
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.667
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 15.658  E-centroid 1.0702 (exact 0.8959)  E outside the J>5% band / total 0.834
    peaks y: E [0.441 0.895 1.551 1.637]  J [0.918]   E max/J max 2.195  width(>50% max) E 0.078 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.827
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 29.656  E-centroid 0.5314 (exact 0.9590)  E outside the J>5% band / total 0.916
    peaks y: E [0.301 0.387 0.410]  J [0.980]   E max/J max 10.629  width(>50% max) E 0.102 J 0.070
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/ba20/cen t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 1.523  sum E/sum J 1.527  |F|/|cH| median 1129.363  flux-weighted <|F|>/<|cH|> 1382.489
  beam: direction error J-weighted mean 95.7 deg  median 74.2  p90 164.3   face |F|/cE > 1: frac 0.9435 (beam) 0.9883 (all)  max 16003.76
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.3362  max E/Jmax 12.1155
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.875
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 41.456  E-centroid 1.2620 (exact 0.8959)  E outside the J>5% band / total 0.937
    peaks y: E [0.863 1.504 1.551 1.637]  J [0.918]   E max/J max 6.334  width(>50% max) E 0.203 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.715
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 19.569  E-centroid 0.8809 (exact 0.9590)  E outside the J>5% band / total 0.857
    peaks y: E [0.121 0.309 0.434 0.895]  J [0.980]   E max/J max 2.421  width(>50% max) E 0.266 J 0.070
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/shd3b/hr shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.232  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.593
  x=0.85  E depth 0.000 edge 0.346  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.608
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/shd3b/hrb shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.232  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.593
  x=0.85  E depth 0.000 edge 0.346  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.608
/resnick/groups/carnegie_poc/jingze/xthinfix_1009/run/beams/shd3b/cen shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.284 edge 0.017  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.144
  x=0.85  E depth 0.208 edge 0.011  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.288
```
