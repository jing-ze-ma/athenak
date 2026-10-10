# NOTE Caltech -> viper: xthinfix-1009 STEEP battery (21795e2a), 10-10

Answers TASK-2026-10-10-caltech-steep-battery.md. Run trees stay on Caltech: /resnick/groups/carnegie_poc/jingze/steep_1010
(run/, logs/; REPO = git archive of 21795e2a in repo/).

STATUS: complete. GPU 4306468 (smoke 4306394); CPU order 8 x 32 cores 4306383-90 (1300 runs, all rc 0), gates 4306391,
beams 4306392, ana 4306393.

## Verdict

1. Pulses: hrs / hrs15 / hrs60 / hrsp equal central to <= 3e-10 at every kappa (X/T/C, be/hesdirk2); hr breaks (k0.128 X 0.69, k12.8 1.07). But in rw_t10 E (hesdirk2 X) hrs, hrsp and hrs15 break at 512 (orders -3.91 / -4.13 / -2.87); only hrs60 converges (1.91 1.97 1.91 = central).
2. atm Hopf L1 at 512: hrs 1.8e-3 (orders 1.76 0.62 -0.40 -0.30, grows), hrs15 3.0e-1 (blows up, -6.58), hrs60 9.3e-4 (grows, as viper), hrsp 3.89e-4 / hrp 3.89e-4 (1.77 1.33 0.72 0.30), about the same as hrx0 3.7e-4. G5 / G3 / G1 PASS for every dc steep arm (G1 hrs60 L1 7.2e-4).
3. hrsp CRASHES: G1 (atm2d, cfl 1e2 / 1e4) and G3 (marsh2d, cfl 1 / 10) segfault at cycle 0-1 (rc 139, null address, ucx handler, no symbols); beams xb20 / ba0 / ba20 DIVERGE (Picard resid 1.9e27 / 2.7e26 / 7.9e25, FATAL at cycle 3 / 10 / 15), confirming viper. hrp (all + plm) does not crash, but Picard is 13-38 vs 3-5 (plm positivity fallbacks 1.2e3-3.3e3).
4. dc steep beams: no FATAL, NON-CONVERGED 0, Picard 4.7-5.8 (hr 4.0-5.0), BiCGStab inner 11-22, LOWER than hr (17-45); superluminal-face fraction as hr (hrs: xb20 0 %, ba0 0.16 %, ba20 0.12 % vs hr 0 / 0.18 / 0.15 %); shd3b hrs umbra 0, edges 0.232 / 0.346 = hr.
5. GPU AG Car A: base (5b304cf9) == all bin data (and repeats bitwise). steep costs the same as all (0.492-0.502 vs 0.495-0.503 s/cycle, Picard 4.30 vs 4.33), and differs from all only at >= 0.8 R_ph (<= 6e-6 inside). steep+plm costs +67 % (0.837 s/cycle, Picard 19.9, max 52) and gives O(1) differences above 1.05 R_ph. Summary: dc steep is cheap and stable but not convergent on rw_t10 / atm (except X0 60 on rw_t10); plm fixes atm but crashes / diverges.

## Caltech adaptations (same bundle issues as on 10-09; scripts copied, not edited in place)

- Gate inputs (atm2d, generated pulse2d, marsh2d) have no `vet_tensor`. A command-line key missing from the input is FATAL, so `fixkeys.py` adds every missing key to a copy of the input. Only vet_tensor was missing; the beam inputs were complete.
- /usr/bin/time is absent on the compute nodes, so a shell stand-in writes the same "%e s wall" line.
- beams() now writes `<run>/cmd.txt` (ATHENA_CPU + the command-line keys) for cylref/collref.
- python: spack 3.11 + numpy/scipy/matplotlib + py-h5py. order ran as 8 round-robin job files (od.py run, 32 cores each).
- Copy: steep/battery_steep_caltech.sh. GPU: only the #SBATCH header + srun (gpuwrap.sh) changed.
- FATAL list: beams xb20/ba0/ba20 hrsp (Picard DIVERGED); gates g1_hrsp_1e2, g1_hrsp_1e4, g3_hrsp_c1, g3_hrsp_c10 (SIGSEGV, rc 139; the 4 'g?|...' lines at the top of gates.txt are bash echoing the crashed commands). No other FATAL; NON-CONVERGED > 0 only in the 3 diverged hrsp beams.

## Binaries

- ATHENA_CPU xthinfix-1009 21795e2a, built-in pgens, MPI (openmpi 5.0.1), Release: `athena_cpu_built_in_pgens_21795e2a`
  md5 `bc3e6892f6d8aa96a97519de6d51fdcf`
- ATHENA_GPU xthinfix-1009 21795e2a, he_star_m1 H200 (CUDA 12.9, HOPPER90, hpcx MPI, host -ffp-contract=off):
  `athena_gpu_he_star_m1_21795e2a_nofma` md5 `32a3ef746051f842fa06736497896fe1`
- ATHENA_BASE rt-integration 5b304cf9, he_star_m1 H200, same flags: `athena_gpu_he_star_m1_5b304cf9_nofma`
  md5 `a9724eea21f363dc55ef0bc664ab3cf4`

## GPU point (job 4306468, 2 H200, hpc-sm-02-10; smoke 4306394 steep+plm nlim 3: rc 0, NON-CONVERGED 0, Picard mean 34)

Adapted only: #SBATCH header, srun `--mpi=pmix -n 2 -c 8 --cpu-bind=cores gpuwrap.sh` (CUDA_VISIBLE_DEVICES =
SLURM_LOCALID); python with h5py for spdiff. All 6 arms rc=0 fatal=0, 30 cycles.

RESULTS_gpu_steep.txt (raw):
```
base_1    0.495 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
all_1     0.503 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
steep_1   0.502 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.300000e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.271318e+00  plm pos-fallbacks -  FATAL 0
steepp_1  0.837 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('1.993333e+01', '5.200000e+01', '0.000000e+00')  inner mean 3.075251e+00  plm pos-fallbacks 1.090000e+02  FATAL 0
base_2    0.496 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
steep_2   0.492 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.300000e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.271318e+00  plm pos-fallbacks -  FATAL 0
-- bin data base_1 vs all_1 (after the parameter header <par_end>), and repeats
base_1 vs all_1: 4 bin files, 0 differ in data
base_1 vs base_2: 4 bin files, 0 differ in data
steep_1 vs steep_2: 4 bin files, 0 differ in data
-- all vs steep / steep+plm by radius (max over angles |a-b|/max|b|), last dump
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 2.9e-07 1.3e-03 1.9e-03 2.0e-07 6.1e-07
agcar3d.hydro_w.00001.bin    velx       9.8e-08 5.0e-08 8.7e-04 7.6e-02 1.3e-02 7.4e-05 5.0e-05
agcar3d.hydro_w.00001.bin    vely       3.0e-06 6.0e-09 2.3e-05 1.8e-01 6.2e-01 1.2e-02 9.8e-03
agcar3d.hydro_w.00001.bin    velz       4.0e-06 5.9e-09 2.1e-05 1.6e-01 4.6e-01 1.5e-02 1.5e-02
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 6.6e-06 4.7e-03 1.9e-03 1.3e-04 1.1e-04
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 2.5e-05 3.6e-03 1.6e-04 7.5e-05 3.7e-05
agcar3d.m1.00001.bin         m1_f1      8.0e-08 2.0e-08 1.2e-05 2.3e-05 3.0e-05 3.1e-05 2.5e-05
agcar3d.m1.00001.bin         m1_f2      4.2e-06 1.4e-08 2.0e-05 1.8e-02 6.0e-03 7.0e-03 5.4e-03
agcar3d.m1.00001.bin         m1_f3      5.8e-06 2.8e-08 2.1e-05 1.8e-02 5.8e-03 6.2e-03 5.5e-03
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 1.4e-06 1.4e-02 1.5e-02 2.6e-02 1.3e-02
agcar3d.hydro_w.00001.bin    velx       9.9e-08 6.2e-08 5.9e-03 2.4e+00 8.6e-02 1.8e-01 5.2e-01
agcar3d.hydro_w.00001.bin    vely       3.4e-06 5.5e-08 3.0e-04 2.6e-01 1.0e+00 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    velz       3.4e-06 5.6e-08 3.0e-04 2.5e-01 1.0e+00 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 4.5e-05 8.4e-02 1.6e-02 5.8e-01 8.4e-01
agcar3d.m1.00001.bin         m1_e       1.0e-07 0.0e+00 1.8e-04 8.2e-02 1.0e-02 2.4e-01 5.2e-01
agcar3d.m1.00001.bin         m1_f1      1.1e-07 4.4e-08 8.9e-05 6.2e-04 3.6e-03 2.3e-01 5.1e-01
agcar3d.m1.00001.bin         m1_f2      4.8e-06 5.1e-08 2.9e-04 2.6e-01 7.2e-01 1.1e+00 1.0e+00
agcar3d.m1.00001.bin         m1_f3      4.7e-06 9.0e-08 2.9e-04 2.6e-01 4.9e-01 1.0e+00 1.0e+00
```

## RESULTS/order_eval.txt (raw)

```
atm          cen  be       S e      lev 32-64-128-256-512  L1 1.82e-03 4.59e-04 1.18e-04 3.86e-05 | ord  1.99  1.96  1.61
atm          cen  be       S f1     lev 32-64-128-256-512  L1 9.17e-01 7.37e-01 1.12e+00 5.44e+00 | ord  0.32 -0.60 -2.28
atm          hr   be       S e      lev 32-64-128-256-512  L1 2.11e-03 9.83e-04 8.98e-04 1.10e-03 | ord  1.10  0.13 -0.29
atm          hr   be       S f1     lev 32-64-128-256-512  L1 1.27e+00 1.91e+01 5.60e+01 9.75e-01 | ord -3.91 -1.55  5.85
atm          hrx0 be       S e      lev 32-64-128-256-512  L1 1.77e-03 4.41e-04 1.11e-04 3.91e-05 | ord  2.00  1.99  1.51
atm          hrx0 be       S f1     lev 32-64-128-256-512  L1 1.62e+00 1.02e+01 3.29e+01 9.05e-01 | ord -2.65 -1.69  5.18
atm          hrs  be       S e      lev 32-64-128-256-512  L1 1.88e-03 6.65e-04 4.44e-04 4.97e-04 | ord  1.50  0.58 -0.16
atm          hrs  be       S f1     lev 32-64-128-256-512  L1 1.44e+00 1.68e+01 5.47e+01 8.50e-01 | ord -3.54 -1.70  6.01
atm          hrs15 be       S e      lev 32-64-128-256-512  L1 2.25e-03 1.13e-03 1.02e-03 2.50e-01 | ord  0.99  0.15 -7.94
atm          hrs15 be       S f1     lev 32-64-128-256-512  L1 1.55e+00 1.79e+01 5.67e+01 2.46e+04 | ord -3.53 -1.66 -8.76
atm          hrs60 be       S e      lev 32-64-128-256-512  L1 1.79e-03 4.83e-04 2.10e-04 2.05e-04 | ord  1.89  1.20  0.03
atm          hrs60 be       S f1     lev 32-64-128-256-512  L1 1.25e+00 1.65e+01 5.32e+01 7.66e-01 | ord -3.72 -1.69  6.12
atm          hrsp be       S e      lev 32-64-128-256-512  L1 2.08e-03 5.78e-04 1.76e-04 6.65e-05 | ord  1.85  1.72  1.40
atm          hrsp be       S f1     lev 32-64-128-256-512  L1 1.72e+00 3.74e+00 1.57e+01 1.55e+00 | ord -1.12 -2.07  3.34
atm          hrp  be       S e      lev 32-64-128-256-512  L1 2.08e-03 5.67e-04 1.69e-04 6.41e-05 | ord  1.87  1.74  1.40
atm          hrp  be       S f1     lev 32-64-128-256-512  L1 1.52e+00 6.04e+00 1.49e+01 1.29e+00 | ord -1.99 -1.30  3.53
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
pulse_k0.128 hrs  be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 hrs  be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 hrs  be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 hrs  be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 hrs  be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 hrs  be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 hrs  hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 hrs  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 hrs  hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 hrs  hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 hrs  hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 hrs  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 hrs15 be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 hrs15 be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 hrs15 be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 hrs15 be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 hrs15 be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 hrs15 be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 hrs15 hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 hrs15 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 hrs15 hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 hrs15 hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 hrs15 hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 hrs15 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 hrs60 be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 hrs60 be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 hrs60 be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 hrs60 be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 hrs60 be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 hrs60 be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 hrs60 hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 hrs60 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 hrs60 hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 hrs60 hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 hrs60 hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 hrs60 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
pulse_k0.128 hrsp be       X e      lev 32-64-128-256-512  L1 7.26e-02 2.14e-02 5.49e-03 1.24e-03 | ord  1.76  1.97  2.15
pulse_k0.128 hrsp be       X f1     lev 32-64-128-256-512  L1 1.04e-01 2.90e-02 7.32e-03 1.67e-03 | ord  1.83  1.99  2.13
pulse_k0.128 hrsp be       T e      lev 1-2-4-8-16  L1 4.31e-02 2.29e-02 1.18e-02 6.02e-03 | ord  0.91  0.95  0.97
pulse_k0.128 hrsp be       T f1     lev 1-2-4-8-16  L1 5.40e-02 2.95e-02 1.55e-02 7.95e-03 | ord  0.87  0.93  0.96
pulse_k0.128 hrsp be       C e      lev 32-64-128-256-512  L1 8.18e-02 2.68e-02 1.16e-02 5.81e-03 | ord  1.61  1.21  1.00
pulse_k0.128 hrsp be       C f1     lev 32-64-128-256-512  L1 1.16e-01 4.63e-02 1.95e-02 8.79e-03 | ord  1.32  1.25  1.15
pulse_k0.128 hrsp hesdirk2 X e      lev 32-64-128-256-512  L1 7.49e-02 2.23e-02 5.80e-03 1.41e-03 | ord  1.75  1.94  2.04
pulse_k0.128 hrsp hesdirk2 X f1     lev 32-64-128-256-512  L1 1.09e-01 3.06e-02 7.81e-03 1.87e-03 | ord  1.83  1.97  2.06
pulse_k0.128 hrsp hesdirk2 T e      lev 1-2-4-8-16  L1 2.38e-03 6.88e-04 2.11e-04 7.09e-05 | ord  1.79  1.71  1.57
pulse_k0.128 hrsp hesdirk2 T f1     lev 1-2-4-8-16  L1 3.03e-03 8.44e-04 2.51e-04 8.28e-05 | ord  1.84  1.75  1.60
pulse_k0.128 hrsp hesdirk2 C e      lev 32-64-128-256-512  L1 7.62e-02 2.26e-02 5.83e-03 1.41e-03 | ord  1.75  1.95  2.04
pulse_k0.128 hrsp hesdirk2 C f1     lev 32-64-128-256-512  L1 1.10e-01 3.09e-02 7.93e-03 1.92e-03 | ord  1.83  1.96  2.05
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
pulse_k12.8  hrs  be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  hrs  be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  hrs  be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  hrs  be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  hrs  be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  hrs  be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  hrs  hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  hrs  hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  hrs  hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  hrs  hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  hrs  hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  hrs  hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  hrs15 be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  hrs15 be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  hrs15 be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  hrs15 be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  hrs15 be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  hrs15 be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  hrs15 hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  hrs15 hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  hrs15 hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  hrs15 hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  hrs15 hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  hrs15 hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  hrs60 be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  hrs60 be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  hrs60 be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  hrs60 be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  hrs60 be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  hrs60 be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  hrs60 hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  hrs60 hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  hrs60 hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  hrs60 hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  hrs60 hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  hrs60 hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
pulse_k12.8  hrsp be       X e      lev 32-64-128-256-512  L1 2.36e-02 6.37e-03 1.33e-03 2.80e-04 | ord  1.89  2.25  2.25
pulse_k12.8  hrsp be       X f1     lev 32-64-128-256-512  L1 5.56e-02 1.63e-02 3.58e-03 7.68e-04 | ord  1.77  2.18  2.22
pulse_k12.8  hrsp be       T e      lev 1-2-4-8-16  L1 1.37e-02 7.28e-03 3.77e-03 1.92e-03 | ord  0.91  0.95  0.97
pulse_k12.8  hrsp be       T f1     lev 1-2-4-8-16  L1 3.21e-02 1.77e-02 9.29e-03 4.77e-03 | ord  0.86  0.93  0.96
pulse_k12.8  hrsp be       C e      lev 32-64-128-256-512  L1 2.81e-02 1.17e-02 4.51e-03 2.04e-03 | ord  1.27  1.37  1.15
pulse_k12.8  hrsp be       C f1     lev 32-64-128-256-512  L1 7.55e-02 3.18e-02 1.25e-02 5.48e-03 | ord  1.24  1.35  1.19
pulse_k12.8  hrsp hesdirk2 X e      lev 32-64-128-256-512  L1 2.44e-02 6.51e-03 1.37e-03 2.89e-04 | ord  1.90  2.25  2.24
pulse_k12.8  hrsp hesdirk2 X f1     lev 32-64-128-256-512  L1 5.69e-02 1.68e-02 3.71e-03 7.99e-04 | ord  1.76  2.18  2.22
pulse_k12.8  hrsp hesdirk2 T e      lev 1-2-4-8-16  L1 3.83e-03 1.85e-03 9.03e-04 4.47e-04 | ord  1.05  1.03  1.02
pulse_k12.8  hrsp hesdirk2 T f1     lev 1-2-4-8-16  L1 8.30e-03 4.06e-03 2.00e-03 9.93e-04 | ord  1.03  1.02  1.01
pulse_k12.8  hrsp hesdirk2 C e      lev 32-64-128-256-512  L1 2.36e-02 7.37e-03 1.91e-03 6.01e-04 | ord  1.68  1.95  1.67
pulse_k12.8  hrsp hesdirk2 C f1     lev 32-64-128-256-512  L1 6.12e-02 1.98e-02 5.35e-03 1.69e-03 | ord  1.63  1.89  1.67
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
pulse_k128   hrs  be       X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrs  be       X f1     lev 32-64-128-256-512  L1 2.55e-02 6.91e-03 1.68e-03 4.02e-04 | ord  1.89  2.04  2.06
pulse_k128   hrs  be       T e      lev 1-2-4-8-16  L1 1.46e-03 7.37e-04 3.70e-04 1.85e-04 | ord  0.99  0.99  1.00
pulse_k128   hrs  be       T f1     lev 1-2-4-8-16  L1 5.40e-03 2.73e-03 1.37e-03 6.86e-04 | ord  0.99  0.99  1.00
pulse_k128   hrs  be       C e      lev 32-64-128-256-512  L1 1.78e-02 5.18e-03 1.42e-03 4.37e-04 | ord  1.79  1.86  1.70
pulse_k128   hrs  be       C f1     lev 32-64-128-256-512  L1 2.51e-02 7.99e-03 2.49e-03 9.03e-04 | ord  1.65  1.68  1.46
pulse_k128   hrs  hesdirk2 X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrs  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.55e-02 6.92e-03 1.68e-03 4.02e-04 | ord  1.88  2.04  2.06
pulse_k128   hrs  hesdirk2 T e      lev 1-2-4-8-16  L1 9.44e-05 3.66e-05 1.58e-05 7.27e-06 | ord  1.37  1.21  1.12
pulse_k128   hrs  hesdirk2 T f1     lev 1-2-4-8-16  L1 2.78e-04 1.05e-04 4.48e-05 2.06e-05 | ord  1.40  1.23  1.12
pulse_k128   hrs  hesdirk2 C e      lev 32-64-128-256-512  L1 1.66e-02 4.49e-03 1.07e-03 2.60e-04 | ord  1.88  2.07  2.04
pulse_k128   hrs  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.56e-02 6.96e-03 1.71e-03 4.15e-04 | ord  1.88  2.03  2.04
pulse_k128   hrs15 be       X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrs15 be       X f1     lev 32-64-128-256-512  L1 2.55e-02 6.91e-03 1.68e-03 4.02e-04 | ord  1.89  2.04  2.06
pulse_k128   hrs15 be       T e      lev 1-2-4-8-16  L1 1.46e-03 7.37e-04 3.70e-04 1.85e-04 | ord  0.99  0.99  1.00
pulse_k128   hrs15 be       T f1     lev 1-2-4-8-16  L1 5.40e-03 2.73e-03 1.37e-03 6.86e-04 | ord  0.99  0.99  1.00
pulse_k128   hrs15 be       C e      lev 32-64-128-256-512  L1 1.78e-02 5.18e-03 1.42e-03 4.37e-04 | ord  1.79  1.86  1.70
pulse_k128   hrs15 be       C f1     lev 32-64-128-256-512  L1 2.51e-02 7.99e-03 2.49e-03 9.03e-04 | ord  1.65  1.68  1.46
pulse_k128   hrs15 hesdirk2 X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrs15 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.55e-02 6.92e-03 1.68e-03 4.02e-04 | ord  1.88  2.04  2.06
pulse_k128   hrs15 hesdirk2 T e      lev 1-2-4-8-16  L1 9.44e-05 3.66e-05 1.58e-05 7.27e-06 | ord  1.37  1.21  1.12
pulse_k128   hrs15 hesdirk2 T f1     lev 1-2-4-8-16  L1 2.78e-04 1.05e-04 4.48e-05 2.06e-05 | ord  1.40  1.23  1.12
pulse_k128   hrs15 hesdirk2 C e      lev 32-64-128-256-512  L1 1.66e-02 4.49e-03 1.07e-03 2.60e-04 | ord  1.88  2.07  2.04
pulse_k128   hrs15 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.56e-02 6.96e-03 1.71e-03 4.15e-04 | ord  1.88  2.03  2.04
pulse_k128   hrs60 be       X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrs60 be       X f1     lev 32-64-128-256-512  L1 2.55e-02 6.91e-03 1.68e-03 4.02e-04 | ord  1.89  2.04  2.06
pulse_k128   hrs60 be       T e      lev 1-2-4-8-16  L1 1.46e-03 7.37e-04 3.70e-04 1.85e-04 | ord  0.99  0.99  1.00
pulse_k128   hrs60 be       T f1     lev 1-2-4-8-16  L1 5.40e-03 2.73e-03 1.37e-03 6.86e-04 | ord  0.99  0.99  1.00
pulse_k128   hrs60 be       C e      lev 32-64-128-256-512  L1 1.78e-02 5.18e-03 1.42e-03 4.37e-04 | ord  1.79  1.86  1.70
pulse_k128   hrs60 be       C f1     lev 32-64-128-256-512  L1 2.51e-02 7.99e-03 2.49e-03 9.03e-04 | ord  1.65  1.68  1.46
pulse_k128   hrs60 hesdirk2 X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrs60 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.55e-02 6.92e-03 1.68e-03 4.02e-04 | ord  1.88  2.04  2.06
pulse_k128   hrs60 hesdirk2 T e      lev 1-2-4-8-16  L1 9.44e-05 3.66e-05 1.58e-05 7.27e-06 | ord  1.37  1.21  1.12
pulse_k128   hrs60 hesdirk2 T f1     lev 1-2-4-8-16  L1 2.78e-04 1.05e-04 4.48e-05 2.06e-05 | ord  1.40  1.23  1.12
pulse_k128   hrs60 hesdirk2 C e      lev 32-64-128-256-512  L1 1.66e-02 4.49e-03 1.07e-03 2.60e-04 | ord  1.88  2.07  2.04
pulse_k128   hrs60 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.56e-02 6.96e-03 1.71e-03 4.15e-04 | ord  1.88  2.03  2.04
pulse_k128   hrsp be       X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrsp be       X f1     lev 32-64-128-256-512  L1 2.55e-02 6.91e-03 1.68e-03 4.02e-04 | ord  1.89  2.04  2.06
pulse_k128   hrsp be       T e      lev 1-2-4-8-16  L1 1.46e-03 7.37e-04 3.70e-04 1.85e-04 | ord  0.99  0.99  1.00
pulse_k128   hrsp be       T f1     lev 1-2-4-8-16  L1 5.40e-03 2.73e-03 1.37e-03 6.86e-04 | ord  0.99  0.99  1.00
pulse_k128   hrsp be       C e      lev 32-64-128-256-512  L1 1.78e-02 5.18e-03 1.42e-03 4.37e-04 | ord  1.79  1.86  1.70
pulse_k128   hrsp be       C f1     lev 32-64-128-256-512  L1 2.51e-02 7.99e-03 2.49e-03 9.03e-04 | ord  1.65  1.68  1.46
pulse_k128   hrsp hesdirk2 X e      lev 32-64-128-256-512  L1 1.64e-02 4.45e-03 1.06e-03 2.55e-04 | ord  1.88  2.07  2.05
pulse_k128   hrsp hesdirk2 X f1     lev 32-64-128-256-512  L1 2.55e-02 6.92e-03 1.68e-03 4.02e-04 | ord  1.88  2.04  2.06
pulse_k128   hrsp hesdirk2 T e      lev 1-2-4-8-16  L1 9.44e-05 3.66e-05 1.58e-05 7.27e-06 | ord  1.37  1.21  1.12
pulse_k128   hrsp hesdirk2 T f1     lev 1-2-4-8-16  L1 2.78e-04 1.05e-04 4.48e-05 2.06e-05 | ord  1.40  1.23  1.12
pulse_k128   hrsp hesdirk2 C e      lev 32-64-128-256-512  L1 1.66e-02 4.49e-03 1.07e-03 2.60e-04 | ord  1.88  2.07  2.04
pulse_k128   hrsp hesdirk2 C f1     lev 32-64-128-256-512  L1 2.56e-02 6.96e-03 1.71e-03 4.15e-04 | ord  1.88  2.03  2.04
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
pulse_k1280  hrs  be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs  be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrs  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrs  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrs  be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  hrs  be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  hrs  hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrs  hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  hrs  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  hrs  hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  hrs15 be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs15 be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrs15 be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrs15 be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrs15 be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  hrs15 be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  hrs15 hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs15 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrs15 hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  hrs15 hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  hrs15 hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs15 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  hrs60 be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs60 be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrs60 be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrs60 be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrs60 be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  hrs60 be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  hrs60 hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs60 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrs60 hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  hrs60 hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  hrs60 hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrs60 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
pulse_k1280  hrsp be       X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrsp be       X f1     lev 32-64-128-256-512  L1 2.32e-02 6.07e-03 1.57e-03 4.15e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrsp be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrsp be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k1280  hrsp be       C e      lev 32-64-128-256-512  L1 1.57e-02 4.61e-03 1.35e-03 4.36e-04 | ord  1.77  1.78  1.63
pulse_k1280  hrsp be       C f1     lev 32-64-128-256-512  L1 2.21e-02 6.83e-03 2.29e-03 8.82e-04 | ord  1.70  1.58  1.38
pulse_k1280  hrsp hesdirk2 X e      lev 32-64-128-256-512  L1 1.45e-02 3.92e-03 9.99e-04 2.62e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrsp hesdirk2 X f1     lev 32-64-128-256-512  L1 2.31e-02 6.07e-03 1.58e-03 4.16e-04 | ord  1.93  1.95  1.92
pulse_k1280  hrsp hesdirk2 T e      lev 1-2-4-8-16  L1 4.69e-05 1.21e-05 3.13e-06 8.37e-07 | ord  1.96  1.94  1.90
pulse_k1280  hrsp hesdirk2 T f1     lev 1-2-4-8-16  L1 1.67e-04 4.28e-05 1.10e-05 2.86e-06 | ord  1.97  1.96  1.94
pulse_k1280  hrsp hesdirk2 C e      lev 32-64-128-256-512  L1 1.46e-02 3.93e-03 1.00e-03 2.63e-04 | ord  1.89  1.97  1.93
pulse_k1280  hrsp hesdirk2 C f1     lev 32-64-128-256-512  L1 2.30e-02 6.07e-03 1.58e-03 4.17e-04 | ord  1.92  1.94  1.92
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
pulse_k12800 hrs  be       X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs  be       X f1     lev 32-64-128-256-512  L1 2.30e-02 5.88e-03 1.48e-03 3.72e-04 | ord  1.97  1.99  1.99
pulse_k12800 hrs  be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrs  be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrs  be       C e      lev 32-64-128-256-512  L1 1.54e-02 4.48e-03 1.28e-03 4.09e-04 | ord  1.79  1.80  1.65
pulse_k12800 hrs  be       C f1     lev 32-64-128-256-512  L1 2.19e-02 6.59e-03 2.16e-03 8.28e-04 | ord  1.73  1.61  1.38
pulse_k12800 hrs  hesdirk2 X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs  hesdirk2 X f1     lev 32-64-128-256-512  L1 2.29e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrs  hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.17e-05 2.95e-06 7.45e-07 | ord  1.98  1.99  1.99
pulse_k12800 hrs  hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.21e-05 1.06e-05 2.68e-06 | ord  1.98  1.99  1.99
pulse_k12800 hrs  hesdirk2 C e      lev 32-64-128-256-512  L1 1.44e-02 3.80e-03 9.41e-04 2.37e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs  hesdirk2 C f1     lev 32-64-128-256-512  L1 2.28e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrs15 be       X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs15 be       X f1     lev 32-64-128-256-512  L1 2.30e-02 5.88e-03 1.48e-03 3.72e-04 | ord  1.97  1.99  1.99
pulse_k12800 hrs15 be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrs15 be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrs15 be       C e      lev 32-64-128-256-512  L1 1.54e-02 4.48e-03 1.28e-03 4.09e-04 | ord  1.79  1.80  1.65
pulse_k12800 hrs15 be       C f1     lev 32-64-128-256-512  L1 2.19e-02 6.59e-03 2.16e-03 8.28e-04 | ord  1.73  1.61  1.38
pulse_k12800 hrs15 hesdirk2 X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs15 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.29e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrs15 hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.17e-05 2.95e-06 7.45e-07 | ord  1.98  1.99  1.99
pulse_k12800 hrs15 hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.21e-05 1.06e-05 2.68e-06 | ord  1.98  1.99  1.99
pulse_k12800 hrs15 hesdirk2 C e      lev 32-64-128-256-512  L1 1.44e-02 3.80e-03 9.41e-04 2.37e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs15 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.28e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrs60 be       X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs60 be       X f1     lev 32-64-128-256-512  L1 2.30e-02 5.88e-03 1.48e-03 3.72e-04 | ord  1.97  1.99  1.99
pulse_k12800 hrs60 be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrs60 be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrs60 be       C e      lev 32-64-128-256-512  L1 1.54e-02 4.48e-03 1.28e-03 4.09e-04 | ord  1.79  1.80  1.65
pulse_k12800 hrs60 be       C f1     lev 32-64-128-256-512  L1 2.19e-02 6.59e-03 2.16e-03 8.28e-04 | ord  1.73  1.61  1.38
pulse_k12800 hrs60 hesdirk2 X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs60 hesdirk2 X f1     lev 32-64-128-256-512  L1 2.29e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrs60 hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.17e-05 2.95e-06 7.45e-07 | ord  1.98  1.99  1.99
pulse_k12800 hrs60 hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.21e-05 1.06e-05 2.68e-06 | ord  1.98  1.99  1.99
pulse_k12800 hrs60 hesdirk2 C e      lev 32-64-128-256-512  L1 1.44e-02 3.80e-03 9.41e-04 2.37e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrs60 hesdirk2 C f1     lev 32-64-128-256-512  L1 2.28e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrsp be       X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrsp be       X f1     lev 32-64-128-256-512  L1 2.30e-02 5.88e-03 1.48e-03 3.72e-04 | ord  1.97  1.99  1.99
pulse_k12800 hrsp be       T e      lev 1-2-4-8-16  L1 1.39e-03 6.99e-04 3.51e-04 1.76e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrsp be       T f1     lev 1-2-4-8-16  L1 5.15e-03 2.60e-03 1.31e-03 6.54e-04 | ord  0.99  0.99  1.00
pulse_k12800 hrsp be       C e      lev 32-64-128-256-512  L1 1.54e-02 4.48e-03 1.28e-03 4.09e-04 | ord  1.79  1.80  1.65
pulse_k12800 hrsp be       C f1     lev 32-64-128-256-512  L1 2.19e-02 6.59e-03 2.16e-03 8.28e-04 | ord  1.73  1.61  1.38
pulse_k12800 hrsp hesdirk2 X e      lev 32-64-128-256-512  L1 1.43e-02 3.79e-03 9.38e-04 2.36e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrsp hesdirk2 X f1     lev 32-64-128-256-512  L1 2.29e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
pulse_k12800 hrsp hesdirk2 T e      lev 1-2-4-8-16  L1 4.62e-05 1.17e-05 2.95e-06 7.45e-07 | ord  1.98  1.99  1.99
pulse_k12800 hrsp hesdirk2 T f1     lev 1-2-4-8-16  L1 1.66e-04 4.21e-05 1.06e-05 2.68e-06 | ord  1.98  1.99  1.99
pulse_k12800 hrsp hesdirk2 C e      lev 32-64-128-256-512  L1 1.44e-02 3.80e-03 9.41e-04 2.37e-04 | ord  1.92  2.01  1.99
pulse_k12800 hrsp hesdirk2 C f1     lev 32-64-128-256-512  L1 2.28e-02 5.88e-03 1.48e-03 3.73e-04 | ord  1.96  1.99  1.99
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
rw_t10       hrs  be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.08e-03 1.21e-02 4.27e-01 | ord  2.17 -1.26 -5.14
rw_t10       hrs  be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.46e-01 2.29e+00 | ord  0.44  0.02 -3.97
rw_t10       hrs  be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.46e-03 6.78e-03 | ord  2.18  2.01 -2.22
rw_t10       hrs  be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.48e-03 1.23e-02 | ord  2.17  2.04 -3.06
rw_t10       hrs  be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.48e-03 7.41e-03 | ord  2.18  1.99 -2.33
rw_t10       hrs  be       T e      lev 1-2-4-8-16  L1 7.77e-02 4.26e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       hrs  be       T f1     lev 1-2-4-8-16  L1 6.62e-02 3.55e-02 1.86e-02 9.24e-03 | ord  0.90  0.94  1.01
rw_t10       hrs  be       T dens   lev 1-2-4-8-16  L1 6.75e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       hrs  be       T velx   lev 1-2-4-8-16  L1 6.40e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       hrs  be       T eint   lev 1-2-4-8-16  L1 6.78e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       hrs  be       C e      lev 32-64-128-256-512  L1 2.03e-01 1.38e-01 1.40e-01 4.01e-01 | ord  0.56 -0.02 -1.52
rw_t10       hrs  be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.25e-01 1.77e-01 2.17e+00 | ord  0.75 -0.50 -3.62
rw_t10       hrs  be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.53e-02 7.84e-03 8.48e-03 | ord  1.52  0.96 -0.11
rw_t10       hrs  be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 7.23e-03 1.16e-02 | ord  1.55  0.92 -0.68
rw_t10       hrs  be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 7.96e-03 8.93e-03 | ord  1.52  0.96 -0.16
rw_t10       hrs  hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.04e-03 1.56e-02 | ord  1.91  1.99 -3.91
rw_t10       hrs  hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.52e-02 4.15e-01 | ord  1.81  0.08 -3.20
rw_t10       hrs  hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.64e-04 | ord  2.18  2.13  2.35
rw_t10       hrs  hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 2.17e-03 | ord  2.16  2.14 -0.64
rw_t10       hrs  hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.58e-04 | ord  2.18  2.13  2.39
rw_t10       hrs  hesdirk2 T e      lev 1-2-4-8-16  L1 8.52e-04 4.30e-04 2.17e-04 1.50e-04 | ord  0.99  0.99  0.54
rw_t10       hrs  hesdirk2 T f1     lev 1-2-4-8-16  L1 6.77e-04 3.86e-04 1.82e-04 2.18e-04 | ord  0.81  1.08 -0.26
rw_t10       hrs  hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrs  hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       hrs  hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrs  hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 1.09e-03 1.39e-02 | ord  2.35  1.37 -3.67
rw_t10       hrs  hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.82e-02 1.96e-02 3.59e-01 | ord  2.48  0.52 -4.19
rw_t10       hrs  hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.69e-03 4.39e-04 | ord  2.12  1.97  1.95
rw_t10       hrs  hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 1.88e-03 | ord  2.16  2.07 -0.30
rw_t10       hrs  hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.69e-03 4.33e-04 | ord  2.12  1.97  1.97
rw_t10       hrs15 be       X e      lev 32-64-128-256-512  L1 2.28e-02 1.32e-02 4.34e-01 1.12e-01 | ord  0.79 -5.04  1.95
rw_t10       hrs15 be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.52e-01 2.08e+00 4.98e-01 | ord  0.40 -3.77  2.06
rw_t10       hrs15 be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.97e-03 7.60e-03 1.17e-03 | ord  2.15 -0.35  2.70
rw_t10       hrs15 be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.17e-03 1.23e-02 1.81e-03 | ord  2.14 -1.00  2.77
rw_t10       hrs15 be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.99e-03 8.28e-03 1.27e-03 | ord  2.15 -0.47  2.70
rw_t10       hrs15 be       T e      lev 1-2-4-8-16  L1 1.07e-01 1.39e-02 1.22e-02 1.07e-02 | ord  2.95  0.18  0.19
rw_t10       hrs15 be       T f1     lev 1-2-4-8-16  L1 1.71e-01 1.28e-01 5.44e-02 1.13e-02 | ord  0.41  1.24  2.27
rw_t10       hrs15 be       T dens   lev 1-2-4-8-16  L1 4.84e-03 2.72e-03 1.56e-03 8.31e-04 | ord  0.83  0.80  0.91
rw_t10       hrs15 be       T velx   lev 1-2-4-8-16  L1 5.12e-03 2.78e-03 1.52e-03 7.94e-04 | ord  0.88  0.87  0.94
rw_t10       hrs15 be       T eint   lev 1-2-4-8-16  L1 4.87e-03 2.72e-03 1.56e-03 8.35e-04 | ord  0.84  0.80  0.90
rw_t10       hrs15 be       C e      lev 32-64-128-256-512  L1 2.05e-01 3.30e-01 9.80e-01 2.72e+01 | ord -0.69 -1.57 -4.79
rw_t10       hrs15 be       C f1     lev 32-64-128-256-512  L1 2.10e-01 3.50e-01 2.43e+00 1.08e+00 | ord -0.73 -2.80  1.18
rw_t10       hrs15 be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.93e-02 1.84e-02 5.84e-03 | ord  1.19  0.07  1.66
rw_t10       hrs15 be       C velx   lev 32-64-128-256-512  L1 4.02e-02 1.69e-02 1.74e-02 4.91e-03 | ord  1.25 -0.05  1.83
rw_t10       hrs15 be       C eint   lev 32-64-128-256-512  L1 4.45e-02 1.99e-02 1.94e-02 6.67e-03 | ord  1.16  0.03  1.54
rw_t10       hrs15 hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.09e-03 1.33e-02 9.70e-02 | ord  1.91 -1.70 -2.87
rw_t10       hrs15 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.03e-01 2.28e+00 | ord  1.82 -3.08 -2.50
rw_t10       hrs15 hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.25e-03 1.31e-03 | ord  2.18  2.24 -0.07
rw_t10       hrs15 hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 2.19e-03 1.20e-02 | ord  2.16  1.49 -2.46
rw_t10       hrs15 hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.23e-03 1.42e-03 | ord  2.18  2.26 -0.21
rw_t10       hrs15 hesdirk2 T e      lev 1-2-4-8-16  L1 1.04e-03 3.17e-04 1.95e-04 1.05e-04 | ord  1.72  0.70  0.90
rw_t10       hrs15 hesdirk2 T f1     lev 1-2-4-8-16  L1 7.17e-02 1.64e-02 1.54e-03 1.28e-04 | ord  2.12  3.42  3.59
rw_t10       hrs15 hesdirk2 T dens   lev 1-2-4-8-16  L1 1.28e-03 6.36e-04 3.16e-04 1.58e-04 | ord  1.02  1.01  1.00
rw_t10       hrs15 hesdirk2 T velx   lev 1-2-4-8-16  L1 1.32e-03 6.41e-04 3.17e-04 1.58e-04 | ord  1.04  1.01  1.01
rw_t10       hrs15 hesdirk2 T eint   lev 1-2-4-8-16  L1 1.28e-03 6.35e-04 3.16e-04 1.58e-04 | ord  1.02  1.01  1.00
rw_t10       hrs15 hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 3.77e-03 5.77e-02 4.96e-02 | ord  1.92 -3.94  0.22
rw_t10       hrs15 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 9.75e-02 1.68e+00 9.03e-01 | ord  0.69 -4.10  0.89
rw_t10       hrs15 hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.59e-03 1.01e-03 3.92e-04 | ord  2.12  2.70  1.37
rw_t10       hrs15 hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.42e-03 8.78e-03 4.92e-03 | ord  2.16 -0.45  0.83
rw_t10       hrs15 hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.59e-03 9.32e-04 4.82e-04 | ord  2.12  2.82  0.95
rw_t10       hrs60 be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.07e-03 1.26e-03 1.18e-02 | ord  2.17  2.01 -3.23
rw_t10       hrs60 be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.37e-01 7.05e-02 | ord  0.44  0.12  0.96
rw_t10       hrs60 be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.37e-04 | ord  2.18  2.14  1.61
rw_t10       hrs60 be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.37e-03 4.68e-04 | ord  2.17  2.15  1.55
rw_t10       hrs60 be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.33e-03 4.57e-04 | ord  2.18  2.14  1.54
rw_t10       hrs60 be       T e      lev 1-2-4-8-16  L1 7.83e-02 4.28e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       hrs60 be       T f1     lev 1-2-4-8-16  L1 6.67e-02 3.59e-02 1.85e-02 9.36e-03 | ord  0.89  0.95  0.99
rw_t10       hrs60 be       T dens   lev 1-2-4-8-16  L1 6.76e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       hrs60 be       T velx   lev 1-2-4-8-16  L1 6.41e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       hrs60 be       T eint   lev 1-2-4-8-16  L1 6.79e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       hrs60 be       C e      lev 32-64-128-256-512  L1 2.04e-01 1.38e-01 7.98e-02 5.37e-02 | ord  0.56  0.79  0.57
rw_t10       hrs60 be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.24e-01 6.76e-02 6.51e-02 | ord  0.75  0.87  0.06
rw_t10       hrs60 be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.52e-02 6.87e-03 3.54e-03 | ord  1.53  1.15  0.96
rw_t10       hrs60 be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 6.48e-03 3.33e-03 | ord  1.55  1.08  0.96
rw_t10       hrs60 be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 6.93e-03 3.57e-03 | ord  1.52  1.16  0.96
rw_t10       hrs60 hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.05e-03 2.79e-04 | ord  1.91  1.97  1.91
rw_t10       hrs60 hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 2.23e-03 | ord  1.82  0.08  4.34
rw_t10       hrs60 hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       hrs60 hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 3.12e-04 | ord  2.16  2.14  2.16
rw_t10       hrs60 hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.93e-04 | ord  2.18  2.13  2.21
rw_t10       hrs60 hesdirk2 T e      lev 1-2-4-8-16  L1 8.55e-04 4.28e-04 4.68e-04 5.21e-04 | ord  1.00 -0.13 -0.15
rw_t10       hrs60 hesdirk2 T f1     lev 1-2-4-8-16  L1 3.99e-04 6.24e-04 2.08e-04 8.66e-05 | ord -0.65  1.59  1.26
rw_t10       hrs60 hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrs60 hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.37e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       hrs60 hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrs60 hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 5.23e-04 3.03e-04 | ord  2.35  2.43  0.79
rw_t10       hrs60 hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.81e-02 4.03e-03 2.36e-03 | ord  2.49  2.80  0.77
rw_t10       hrs60 hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       hrs60 hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.53e-03 6.16e-04 | ord  2.16  2.07  1.31
rw_t10       hrs60 hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.70e-03 6.13e-04 | ord  2.12  1.96  1.47
rw_t10       hrsp be       X e      lev 32-64-128-256-512  L1 2.28e-02 5.08e-03 1.13e-02 4.18e-01 | ord  2.17 -1.16 -5.21
rw_t10       hrsp be       X f1     lev 32-64-128-256-512  L1 2.01e-01 1.48e-01 1.46e-01 2.29e+00 | ord  0.44  0.03 -3.97
rw_t10       hrsp be       X dens   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.45e-03 6.66e-03 | ord  2.18  2.01 -2.20
rw_t10       hrsp be       X velx   lev 32-64-128-256-512  L1 2.72e-02 6.07e-03 1.47e-03 1.23e-02 | ord  2.17  2.04 -3.07
rw_t10       hrsp be       X eint   lev 32-64-128-256-512  L1 2.65e-02 5.85e-03 1.47e-03 7.27e-03 | ord  2.18  2.00 -2.31
rw_t10       hrsp be       T e      lev 1-2-4-8-16  L1 7.77e-02 4.26e-02 2.24e-02 1.15e-02 | ord  0.87  0.93  0.96
rw_t10       hrsp be       T f1     lev 1-2-4-8-16  L1 6.62e-02 3.57e-02 1.85e-02 9.24e-03 | ord  0.89  0.95  1.00
rw_t10       hrsp be       T dens   lev 1-2-4-8-16  L1 6.75e-03 3.37e-03 1.68e-03 8.41e-04 | ord  1.00  1.00  1.00
rw_t10       hrsp be       T velx   lev 1-2-4-8-16  L1 6.40e-03 3.20e-03 1.60e-03 8.01e-04 | ord  1.00  1.00  1.00
rw_t10       hrsp be       T eint   lev 1-2-4-8-16  L1 6.78e-03 3.39e-03 1.69e-03 8.45e-04 | ord  1.00  1.00  1.00
rw_t10       hrsp be       C e      lev 32-64-128-256-512  L1 2.03e-01 1.38e-01 1.38e-01 3.94e-01 | ord  0.56 -0.00 -1.51
rw_t10       hrsp be       C f1     lev 32-64-128-256-512  L1 2.09e-01 1.25e-01 1.76e-01 2.17e+00 | ord  0.75 -0.49 -3.62
rw_t10       hrsp be       C dens   lev 32-64-128-256-512  L1 4.39e-02 1.53e-02 7.81e-03 8.38e-03 | ord  1.52  0.97 -0.10
rw_t10       hrsp be       C velx   lev 32-64-128-256-512  L1 4.01e-02 1.37e-02 7.21e-03 1.16e-02 | ord  1.55  0.92 -0.69
rw_t10       hrsp be       C eint   lev 32-64-128-256-512  L1 4.44e-02 1.55e-02 7.93e-03 8.82e-03 | ord  1.52  0.97 -0.15
rw_t10       hrsp hesdirk2 X e      lev 32-64-128-256-512  L1 1.54e-02 4.10e-03 1.03e-03 1.80e-02 | ord  1.91  1.99 -4.13
rw_t10       hrsp hesdirk2 X f1     lev 32-64-128-256-512  L1 1.68e-01 4.76e-02 4.51e-02 4.14e-01 | ord  1.82  0.08 -3.20
rw_t10       hrsp hesdirk2 X dens   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.78e-04 | ord  2.18  2.13  2.28
rw_t10       hrsp hesdirk2 X velx   lev 32-64-128-256-512  L1 2.74e-02 6.13e-03 1.39e-03 2.16e-03 | ord  2.16  2.14 -0.64
rw_t10       hrsp hesdirk2 X eint   lev 32-64-128-256-512  L1 2.67e-02 5.90e-03 1.35e-03 2.74e-04 | ord  2.18  2.13  2.30
rw_t10       hrsp hesdirk2 T e      lev 1-2-4-8-16  L1 8.48e-04 4.27e-04 2.16e-04 1.50e-04 | ord  0.99  0.98  0.53
rw_t10       hrsp hesdirk2 T f1     lev 1-2-4-8-16  L1 6.42e-04 4.24e-04 2.17e-04 2.48e-04 | ord  0.60  0.96 -0.19
rw_t10       hrsp hesdirk2 T dens   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrsp hesdirk2 T velx   lev 1-2-4-8-16  L1 1.29e-03 6.38e-04 3.17e-04 1.58e-04 | ord  1.02  1.01  1.01
rw_t10       hrsp hesdirk2 T eint   lev 1-2-4-8-16  L1 1.27e-03 6.32e-04 3.16e-04 1.58e-04 | ord  1.01  1.00  1.00
rw_t10       hrsp hesdirk2 C e      lev 32-64-128-256-512  L1 1.43e-02 2.81e-03 1.29e-03 1.58e-02 | ord  2.35  1.12 -3.62
rw_t10       hrsp hesdirk2 C f1     lev 32-64-128-256-512  L1 1.57e-01 2.82e-02 1.96e-02 3.58e-01 | ord  2.48  0.53 -4.19
rw_t10       hrsp hesdirk2 C dens   lev 32-64-128-256-512  L1 2.86e-02 6.61e-03 1.69e-03 4.14e-04 | ord  2.12  1.97  2.03
rw_t10       hrsp hesdirk2 C velx   lev 32-64-128-256-512  L1 2.87e-02 6.43e-03 1.52e-03 1.87e-03 | ord  2.16  2.08 -0.30
rw_t10       hrsp hesdirk2 C eint   lev 32-64-128-256-512  L1 2.87e-02 6.61e-03 1.69e-03 4.09e-04 | ord  2.12  1.97  2.05
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
rw_t1000     hrs  be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrs  be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     hrs  be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     hrs  be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     hrs  be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     hrs  be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     hrs  be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     hrs  be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     hrs  be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     hrs  be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     hrs  be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     hrs  be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     hrs  be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     hrs  be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     hrs  be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     hrs  hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrs  hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     hrs  hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     hrs  hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     hrs  hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     hrs  hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     hrs  hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     hrs  hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     hrs  hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     hrs  hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     hrs  hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     hrs  hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     hrs  hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     hrs  hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     hrs  hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     hrs15 be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrs15 be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     hrs15 be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     hrs15 be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     hrs15 be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     hrs15 be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     hrs15 be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     hrs15 be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     hrs15 be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     hrs15 be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     hrs15 be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     hrs15 be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     hrs15 be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     hrs15 be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     hrs15 be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     hrs15 hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrs15 hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     hrs15 hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     hrs15 hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     hrs15 hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     hrs15 hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     hrs15 hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     hrs15 hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     hrs15 hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     hrs15 hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     hrs15 hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     hrs15 hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     hrs15 hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     hrs15 hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     hrs15 hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     hrs60 be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrs60 be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     hrs60 be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     hrs60 be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     hrs60 be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     hrs60 be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     hrs60 be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     hrs60 be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     hrs60 be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     hrs60 be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     hrs60 be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     hrs60 be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     hrs60 be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     hrs60 be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     hrs60 be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     hrs60 hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrs60 hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     hrs60 hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     hrs60 hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     hrs60 hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     hrs60 hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     hrs60 hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     hrs60 hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     hrs60 hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     hrs60 hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     hrs60 hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     hrs60 hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     hrs60 hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     hrs60 hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     hrs60 hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
rw_t1000     hrsp be       X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrsp be       X f1     lev 32-64-128-256  L1 4.93e-03 1.47e-02 5.91e-04 | ord -1.58  4.64
rw_t1000     hrsp be       X dens   lev 32-64-128-256  L1 8.73e-03 2.04e-03 5.31e-04 | ord  2.10  1.94
rw_t1000     hrsp be       X velx   lev 32-64-128-256  L1 8.40e-03 2.03e-03 5.36e-04 | ord  2.05  1.92
rw_t1000     hrsp be       X eint   lev 32-64-128-256  L1 8.69e-03 2.02e-03 5.30e-04 | ord  2.10  1.93
rw_t1000     hrsp be       T e      lev 2-4-8-16  L1 3.30e-03 1.68e-03 8.50e-04 | ord  0.97  0.98
rw_t1000     hrsp be       T f1     lev 2-4-8-16  L1 2.77e-02 1.69e-03 8.52e-04 | ord  4.03  0.99
rw_t1000     hrsp be       T dens   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.41e-04 | ord  0.98  0.99
rw_t1000     hrsp be       T velx   lev 2-4-8-16  L1 3.31e-03 1.68e-03 8.50e-04 | ord  0.97  0.99
rw_t1000     hrsp be       T eint   lev 2-4-8-16  L1 3.28e-03 1.67e-03 8.43e-04 | ord  0.97  0.99
rw_t1000     hrsp be       C e      lev 32-64-128-256  L1 2.49e-02 1.15e-02 5.58e-03 | ord  1.12  1.04
rw_t1000     hrsp be       C f1     lev 32-64-128-256  L1 2.62e-02 1.31e-02 2.99e-02 | ord  1.00 -1.19
rw_t1000     hrsp be       C dens   lev 32-64-128-256  L1 2.62e-02 1.17e-02 5.62e-03 | ord  1.16  1.06
rw_t1000     hrsp be       C velx   lev 32-64-128-256  L1 2.51e-02 1.17e-02 5.67e-03 | ord  1.10  1.05
rw_t1000     hrsp be       C eint   lev 32-64-128-256  L1 2.56e-02 1.16e-02 5.60e-03 | ord  1.14  1.05
rw_t1000     hrsp hesdirk2 X e      lev 32-64-128-256  L1 8.49e-03 2.03e-03 5.36e-04 | ord  2.06  1.92
rw_t1000     hrsp hesdirk2 X f1     lev 32-64-128-256  L1 5.75e-03 5.87e-03 5.91e-04 | ord -0.03  3.31
rw_t1000     hrsp hesdirk2 X dens   lev 32-64-128-256  L1 8.76e-03 2.04e-03 5.27e-04 | ord  2.11  1.95
rw_t1000     hrsp hesdirk2 X velx   lev 32-64-128-256  L1 8.44e-03 2.04e-03 5.37e-04 | ord  2.05  1.92
rw_t1000     hrsp hesdirk2 X eint   lev 32-64-128-256  L1 8.72e-03 2.01e-03 5.27e-04 | ord  2.12  1.93
rw_t1000     hrsp hesdirk2 T e      lev 2-4-8-16  L1 1.14e-02 5.66e-03 2.81e-03 | ord  1.01  1.01
rw_t1000     hrsp hesdirk2 T f1     lev 2-4-8-16  L1 1.20e-02 5.54e-03 2.79e-03 | ord  1.12  0.99
rw_t1000     hrsp hesdirk2 T dens   lev 2-4-8-16  L1 1.16e-02 5.78e-03 2.88e-03 | ord  1.01  1.00
rw_t1000     hrsp hesdirk2 T velx   lev 2-4-8-16  L1 1.09e-02 5.54e-03 2.78e-03 | ord  0.98  0.99
rw_t1000     hrsp hesdirk2 T eint   lev 2-4-8-16  L1 1.14e-02 5.68e-03 2.83e-03 | ord  1.01  1.00
rw_t1000     hrsp hesdirk2 C e      lev 32-64-128-256  L1 7.99e-02 3.81e-02 1.84e-02 | ord  1.07  1.05
rw_t1000     hrsp hesdirk2 C f1     lev 32-64-128-256  L1 6.32e-02 3.35e-02 1.78e-02 | ord  0.92  0.91
rw_t1000     hrsp hesdirk2 C dens   lev 32-64-128-256  L1 7.95e-02 3.85e-02 1.88e-02 | ord  1.05  1.04
rw_t1000     hrsp hesdirk2 C velx   lev 32-64-128-256  L1 6.43e-02 3.36e-02 1.72e-02 | ord  0.94  0.96
rw_t1000     hrsp hesdirk2 C eint   lev 32-64-128-256  L1 7.92e-02 3.82e-02 1.85e-02 | ord  1.05  1.04
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
hrx0 n=  32  L1 5.936e-03  max|dE/E| 1.536e-02  max(tau<1) 1.536e-02  max(tau>1) 6.389e-03  top 8.197e-04  tau_cell top 4.0e-03 bot 2.8e+01
hrx0 n=  64  L1 1.684e-03  max|dE/E| 1.732e-02  max(tau<1) 1.732e-02  max(tau>1) 5.168e-03  top 2.821e-04  tau_cell top 2.0e-03 bot 1.4e+01
hrx0 n= 128  L1 6.658e-04  max|dE/E| 1.840e-02  max(tau<1) 1.840e-02  max(tau>1) 5.471e-03  top 8.283e-05  tau_cell top 1.0e-03 bot 7.0e+00
hrx0 n= 256  L1 4.304e-04  max|dE/E| 1.893e-02  max(tau<1) 1.893e-02  max(tau>1) 4.958e-03  top 1.520e-05  tau_cell top 5.0e-04 bot 3.5e+00
hrx0 n= 512  L1 3.714e-04  max|dE/E| 1.928e-02  max(tau<1) 1.928e-02  max(tau>1) 4.742e-03  top 3.278e-05  tau_cell top 2.5e-04 bot 1.7e+00
hrx0 Hopf L1 orders 1.82 1.34 0.63 0.21
hrx0 self m1_e 1.77e-03 4.41e-04 1.11e-04 3.91e-05 | ord 2.00 1.99 1.51
hrx0 self m1_f1 1.62e+00 1.02e+01 3.29e+01 9.05e-01 | ord -2.65 -1.69 5.18
hrs n=  32  L1 5.771e-03  max|dE/E| 1.818e-02  max(tau<1) 1.818e-02  max(tau>1) 1.065e-02  top 3.696e-04  tau_cell top 4.0e-03 bot 2.8e+01
hrs n=  64  L1 1.698e-03  max|dE/E| 1.952e-02  max(tau<1) 1.952e-02  max(tau>1) 1.219e-02  top 1.573e-04  tau_cell top 2.0e-03 bot 1.4e+01
hrs n= 128  L1 1.104e-03  max|dE/E| 2.173e-02  max(tau<1) 2.173e-02  max(tau>1) 1.610e-02  top 6.115e-04  tau_cell top 1.0e-03 bot 7.0e+00
hrs n= 256  L1 1.460e-03  max|dE/E| 1.993e-02  max(tau<1) 1.993e-02  max(tau>1) 1.853e-02  top 1.128e-03  tau_cell top 5.0e-04 bot 3.5e+00
hrs n= 512  L1 1.795e-03  max|dE/E| 2.046e-02  max(tau<1) 2.046e-02  max(tau>1) 1.322e-02  top 1.042e-03  tau_cell top 2.5e-04 bot 1.7e+00
hrs Hopf L1 orders 1.76 0.62 -0.40 -0.30
hrs self m1_e 1.88e-03 6.65e-04 4.44e-04 4.97e-04 | ord 1.50 0.58 -0.16
hrs self m1_f1 1.44e+00 1.68e+01 5.47e+01 8.50e-01 | ord -3.54 -1.70 6.01
hrs15 n=  32  L1 5.788e-03  max|dE/E| 3.193e-02  max(tau<1) 3.193e-02  max(tau>1) 1.931e-02  top 6.725e-04  tau_cell top 4.0e-03 bot 2.8e+01
hrs15 n=  64  L1 2.138e-03  max|dE/E| 3.649e-02  max(tau<1) 3.649e-02  max(tau>1) 2.430e-02  top 1.358e-03  tau_cell top 2.0e-03 bot 1.4e+01
hrs15 n= 128  L1 2.352e-03  max|dE/E| 3.327e-02  max(tau<1) 3.327e-02  max(tau>1) 3.247e-02  top 2.409e-03  tau_cell top 1.0e-03 bot 7.0e+00
hrs15 n= 256  L1 3.143e-03  max|dE/E| 2.468e-02  max(tau<1) 2.124e-02  max(tau>1) 2.468e-02  top 2.120e-03  tau_cell top 5.0e-04 bot 3.5e+00
hrs15 n= 512  L1 3.004e-01  max|dE/E| 4.319e+00  max(tau<1) 8.768e-01  max(tau>1) 4.319e+00  top 6.812e-01  tau_cell top 2.5e-04 bot 1.7e+00
hrs15 Hopf L1 orders 1.44 -0.14 -0.42 -6.58
hrs15 self m1_e 2.25e-03 1.13e-03 1.02e-03 2.50e-01 | ord 0.99 0.15 -7.94
hrs15 self m1_f1 1.55e+00 1.79e+01 5.67e+01 2.46e+04 | ord -3.53 -1.66 -8.76
hrs60 n=  32  L1 5.838e-03  max|dE/E| 1.295e-02  max(tau<1) 1.295e-02  max(tau>1) 7.374e-03  top 5.265e-04  tau_cell top 4.0e-03 bot 2.8e+01
hrs60 n=  64  L1 1.609e-03  max|dE/E| 1.319e-02  max(tau<1) 1.319e-02  max(tau>1) 7.749e-03  top 2.454e-05  tau_cell top 2.0e-03 bot 1.4e+01
hrs60 n= 128  L1 7.240e-04  max|dE/E| 1.632e-02  max(tau<1) 1.632e-02  max(tau>1) 9.329e-03  top 4.560e-06  tau_cell top 1.0e-03 bot 7.0e+00
hrs60 n= 256  L1 7.462e-04  max|dE/E| 1.848e-02  max(tau<1) 1.848e-02  max(tau>1) 1.009e-02  top 2.181e-04  tau_cell top 5.0e-04 bot 3.5e+00
hrs60 n= 512  L1 9.325e-04  max|dE/E| 1.979e-02  max(tau<1) 1.979e-02  max(tau>1) 1.165e-02  top 5.354e-04  tau_cell top 2.5e-04 bot 1.7e+00
hrs60 Hopf L1 orders 1.86 1.15 -0.04 -0.32
hrs60 self m1_e 1.79e-03 4.83e-04 2.10e-04 2.05e-04 | ord 1.89 1.20 0.03
hrs60 self m1_f1 1.25e+00 1.65e+01 5.32e+01 7.66e-01 | ord -3.72 -1.69 6.12
hrsp n=  32  L1 6.739e-03  max|dE/E| 2.942e-02  max(tau<1) 2.942e-02  max(tau>1) 6.909e-03  top 1.500e-03  tau_cell top 4.0e-03 bot 2.8e+01
hrsp n=  64  L1 1.983e-03  max|dE/E| 2.377e-02  max(tau<1) 2.377e-02  max(tau>1) 2.184e-03  top 6.249e-04  tau_cell top 2.0e-03 bot 1.4e+01
hrsp n= 128  L1 7.860e-04  max|dE/E| 2.146e-02  max(tau<1) 2.146e-02  max(tau>1) 2.868e-03  top 2.342e-04  tau_cell top 1.0e-03 bot 7.0e+00
hrsp n= 256  L1 4.779e-04  max|dE/E| 2.049e-02  max(tau<1) 2.049e-02  max(tau>1) 3.567e-03  top 6.504e-05  tau_cell top 5.0e-04 bot 3.5e+00
hrsp n= 512  L1 3.892e-04  max|dE/E| 2.006e-02  max(tau<1) 2.006e-02  max(tau>1) 4.051e-03  top 6.607e-06  tau_cell top 2.5e-04 bot 1.7e+00
hrsp Hopf L1 orders 1.77 1.33 0.72 0.30
hrsp self m1_e 2.08e-03 5.78e-04 1.76e-04 6.65e-05 | ord 1.85 1.72 1.40
hrsp self m1_f1 1.72e+00 3.74e+00 1.57e+01 1.55e+00 | ord -1.12 -2.07 3.34
hrp n=  32  L1 6.712e-03  max|dE/E| 2.845e-02  max(tau<1) 2.845e-02  max(tau>1) 6.863e-03  top 1.434e-03  tau_cell top 4.0e-03 bot 2.8e+01
hrp n=  64  L1 1.965e-03  max|dE/E| 2.356e-02  max(tau<1) 2.356e-02  max(tau>1) 2.115e-03  top 6.013e-04  tau_cell top 2.0e-03 bot 1.4e+01
hrp n= 128  L1 7.818e-04  max|dE/E| 2.145e-02  max(tau<1) 2.145e-02  max(tau>1) 3.003e-03  top 2.343e-04  tau_cell top 1.0e-03 bot 7.0e+00
hrp n= 256  L1 4.774e-04  max|dE/E| 2.050e-02  max(tau<1) 2.050e-02  max(tau>1) 3.645e-03  top 6.986e-05  tau_cell top 5.0e-04 bot 3.5e+00
hrp n= 512  L1 3.888e-04  max|dE/E| 2.006e-02  max(tau<1) 2.006e-02  max(tau>1) 4.067e-03  top 4.807e-06  tau_cell top 2.5e-04 bot 1.7e+00
hrp Hopf L1 orders 1.77 1.33 0.71 0.30
hrp self m1_e 2.08e-03 5.67e-04 1.69e-04 6.41e-05 | ord 1.87 1.74 1.40
hrp self m1_f1 1.52e+00 6.04e+00 1.49e+01 1.29e+00 | ord -1.99 -1.30 3.53
```

## RESULTS/arm_minus_cen.txt (raw)

```
pulse_k0.128  hr   hesdirk2 X |hr-cen| 1.43e-05 6.13e-05 2.56e-04 1.05e-03 4.23e-03
pulse_k0.128  hr   hesdirk2 C |hr-cen| 1.05e-04 2.29e-04 4.80e-04 9.82e-04 1.99e-03
pulse_k0.128  hr   hesdirk2 T |hr-cen| 7.83e-03 3.92e-03 1.96e-03 9.82e-04 4.91e-04
pulse_k0.128  hr   be       X |hr-cen| 4.92e-05 2.10e-04 8.77e-04 3.58e-03 1.44e-02
pulse_k0.128  hr   be       C |hr-cen| 3.73e-04 7.89e-04 1.65e-03 3.35e-03 6.76e-03
pulse_k0.128  hr   be       T |hr-cen| 2.67e-02 1.34e-02 6.70e-03 3.35e-03 1.68e-03
pulse_k0.128  hrs  hesdirk2 X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 2.01e-14
pulse_k0.128  hrs  hesdirk2 C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrs  hesdirk2 T |hrs-cen| 5.08e-13 1.24e-14 1.40e-14 0.00e+00 0.00e+00
pulse_k0.128  hrs  be       X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 9.31e-14
pulse_k0.128  hrs  be       C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrs  be       T |hrs-cen| 1.63e-10 1.27e-12 3.21e-14 0.00e+00 0.00e+00
pulse_k0.128  hrs15 hesdirk2 X |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 1.57e-14 2.13e-14
pulse_k0.128  hrs15 hesdirk2 C |hrs15-cen| 3.74e-15 4.67e-15 1.02e-14 1.63e-14 2.48e-14
pulse_k0.128  hrs15 hesdirk2 T |hrs15-cen| 1.23e-10 5.95e-13 1.58e-14 1.63e-14 0.00e+00
pulse_k0.128  hrs15 be       X |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 6.50e-14 8.09e-12
pulse_k0.128  hrs15 be       C |hrs15-cen| 4.84e-15 9.76e-15 2.09e-14 4.38e-14 8.35e-14
pulse_k0.128  hrs15 be       T |hrs15-cen| 4.18e-08 3.25e-10 2.53e-12 4.38e-14 0.00e+00
pulse_k0.128  hrs60 hesdirk2 X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrs60 hesdirk2 C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrs60 hesdirk2 T |hrs60-cen| 9.75e-15 1.33e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrs60 be       X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrs60 be       C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrs60 be       T |hrs60-cen| 6.45e-13 2.36e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrsp hesdirk2 X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 2.17e-14
pulse_k0.128  hrsp hesdirk2 C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrsp hesdirk2 T |hrsp-cen| 7.83e-13 1.38e-14 1.31e-14 0.00e+00 0.00e+00
pulse_k0.128  hrsp be       X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 9.30e-14
pulse_k0.128  hrsp be       C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k0.128  hrsp be       T |hrsp-cen| 1.67e-10 1.30e-12 3.23e-14 0.00e+00 0.00e+00
pulse_k12.8   hr   hesdirk2 X |hr-cen| 3.38e-06 1.21e-05 4.96e-05 2.02e-04 8.20e-04
pulse_k12.8   hr   hesdirk2 C |hr-cen| 2.59e-05 4.55e-05 9.23e-05 1.90e-04 3.87e-04
pulse_k12.8   hr   hesdirk2 T |hr-cen| 1.42e-03 7.37e-04 3.76e-04 1.90e-04 9.55e-05
pulse_k12.8   hr   be       X |hr-cen| 1.12e-05 4.03e-05 1.66e-04 6.78e-04 2.74e-03
pulse_k12.8   hr   be       C |hr-cen| 7.07e-05 1.40e-04 3.02e-04 6.36e-04 1.31e-03
pulse_k12.8   hr   be       T |hr-cen| 4.17e-03 2.34e-03 1.24e-03 6.36e-04 3.23e-04
pulse_k12.8   hrs  hesdirk2 X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 1.21e-14
pulse_k12.8   hrs  hesdirk2 C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrs  hesdirk2 T |hrs-cen| 1.53e-12 1.16e-14 7.51e-15 0.00e+00 0.00e+00
pulse_k12.8   hrs  be       X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 6.63e-14
pulse_k12.8   hrs  be       C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrs  be       T |hrs-cen| 1.76e-11 1.84e-13 1.39e-14 0.00e+00 0.00e+00
pulse_k12.8   hrs15 hesdirk2 X |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 1.06e-14 2.14e-14
pulse_k12.8   hrs15 hesdirk2 C |hrs15-cen| 8.48e-16 2.84e-15 5.46e-15 1.03e-14 1.65e-14
pulse_k12.8   hrs15 hesdirk2 T |hrs15-cen| 3.84e-10 1.99e-12 1.33e-14 1.03e-14 0.00e+00
pulse_k12.8   hrs15 be       X |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 3.08e-14 1.44e-12
pulse_k12.8   hrs15 be       C |hrs15-cen| 2.30e-15 4.96e-15 9.35e-15 1.78e-14 3.50e-14
pulse_k12.8   hrs15 be       T |hrs15-cen| 4.50e-09 4.70e-11 4.52e-13 1.78e-14 0.00e+00
pulse_k12.8   hrs60 hesdirk2 X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrs60 hesdirk2 C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrs60 hesdirk2 T |hrs60-cen| 7.53e-15 7.10e-15 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrs60 be       X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrs60 be       C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrs60 be       T |hrs60-cen| 7.64e-14 2.32e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrsp hesdirk2 X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 1.22e-14
pulse_k12.8   hrsp hesdirk2 C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrsp hesdirk2 T |hrsp-cen| 1.40e-12 1.08e-14 8.17e-15 0.00e+00 0.00e+00
pulse_k12.8   hrsp be       X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 6.59e-14
pulse_k12.8   hrsp be       C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12.8   hrsp be       T |hrsp-cen| 1.84e-11 1.91e-13 1.39e-14 0.00e+00 0.00e+00
pulse_k128    hr   hesdirk2 X |hr-cen| 2.66e-06 5.35e-06 1.07e-05 2.13e-05 4.21e-05
pulse_k128    hr   hesdirk2 C |hr-cen| 1.51e-05 1.78e-05 1.94e-05 2.00e-05 1.97e-05
pulse_k128    hr   hesdirk2 T |hr-cen| 1.27e-04 7.26e-05 3.89e-05 2.00e-05 1.01e-05
pulse_k128    hr   be       X |hr-cen| 7.88e-06 1.58e-05 3.15e-05 6.20e-05 1.19e-04
pulse_k128    hr   be       C |hr-cen| 2.80e-05 4.04e-05 5.14e-05 5.87e-05 6.11e-05
pulse_k128    hr   be       T |hr-cen| 2.28e-04 1.61e-04 1.02e-04 5.87e-05 3.18e-05
pulse_k128    hrs  hesdirk2 X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 3.86e-15 5.09e-15
pulse_k128    hrs  hesdirk2 C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 3.23e-15 5.20e-15
pulse_k128    hrs  hesdirk2 T |hrs-cen| 9.13e-14 1.58e-13 2.61e-14 3.23e-15 0.00e+00
pulse_k128    hrs  be       X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 5.82e-14 1.06e-13
pulse_k128    hrs  be       C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 5.72e-14 1.86e-13
pulse_k128    hrs  be       T |hrs-cen| 9.93e-13 1.15e-13 5.56e-14 5.72e-14 0.00e+00
pulse_k128    hrs15 hesdirk2 X |hrs15-cen| 0.00e+00 0.00e+00 2.02e-15 4.19e-15 5.56e-13
pulse_k128    hrs15 hesdirk2 C |hrs15-cen| 1.14e-16 1.96e-15 2.81e-15 2.80e-15 5.11e-15
pulse_k128    hrs15 hesdirk2 T |hrs15-cen| 2.16e-11 1.33e-12 3.71e-13 2.80e-15 2.68e-15
pulse_k128    hrs15 be       X |hrs15-cen| 0.00e+00 0.00e+00 5.51e-14 5.88e-14 4.18e-12
pulse_k128    hrs15 be       C |hrs15-cen| 3.96e-15 2.87e-14 3.44e-14 1.54e-13 2.45e-13
pulse_k128    hrs15 be       T |hrs15-cen| 2.54e-10 2.43e-11 1.10e-12 1.54e-13 8.92e-14
pulse_k128    hrs60 hesdirk2 X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 5.02e-15
pulse_k128    hrs60 hesdirk2 C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrs60 hesdirk2 T |hrs60-cen| 3.74e-14 3.83e-15 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrs60 be       X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 1.63e-13
pulse_k128    hrs60 be       C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrs60 be       T |hrs60-cen| 5.97e-15 4.27e-14 0.00e+00 0.00e+00 0.00e+00
pulse_k128    hrsp hesdirk2 X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 3.51e-15 5.34e-15
pulse_k128    hrsp hesdirk2 C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 3.01e-15 5.14e-15
pulse_k128    hrsp hesdirk2 T |hrsp-cen| 4.70e-14 1.58e-13 2.87e-14 3.01e-15 0.00e+00
pulse_k128    hrsp be       X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 5.71e-14 9.35e-14
pulse_k128    hrsp be       C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 5.77e-14 1.86e-13
pulse_k128    hrsp be       T |hrsp-cen| 1.74e-14 2.57e-14 5.50e-14 5.77e-14 0.00e+00
pulse_k1280   hr   hesdirk2 X |hr-cen| 4.12e-06 8.16e-06 1.62e-05 3.24e-05 6.49e-05
pulse_k1280   hr   hesdirk2 C |hr-cen| 4.71e-06 9.11e-06 1.74e-05 3.21e-05 5.58e-05
pulse_k1280   hr   hesdirk2 T |hr-cen| 3.72e-05 3.63e-05 3.48e-05 3.21e-05 2.79e-05
pulse_k1280   hr   be       X |hr-cen| 4.58e-06 9.06e-06 1.80e-05 3.60e-05 7.21e-05
pulse_k1280   hr   be       C |hr-cen| 4.74e-06 9.35e-06 1.84e-05 3.59e-05 6.85e-05
pulse_k1280   hr   be       T |hr-cen| 3.74e-05 3.73e-05 3.68e-05 3.59e-05 3.42e-05
pulse_k1280   hrs  hesdirk2 X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 1.53e-14
pulse_k1280   hrs  hesdirk2 C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 5.09e-13
pulse_k1280   hrs  hesdirk2 T |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrs  be       X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 5.90e-14
pulse_k1280   hrs  be       C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 3.75e-15
pulse_k1280   hrs  be       T |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrs15 hesdirk2 X |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 9.02e-13 3.75e-13
pulse_k1280   hrs15 hesdirk2 C |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 5.22e-13 1.34e-12
pulse_k1280   hrs15 hesdirk2 T |hrs15-cen| 1.09e-13 1.90e-13 4.54e-13 5.22e-13 2.59e-13
pulse_k1280   hrs15 be       X |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 2.21e-14 3.20e-13
pulse_k1280   hrs15 be       C |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 7.38e-14 2.19e-13
pulse_k1280   hrs15 be       T |hrs15-cen| 9.06e-14 1.57e-13 6.89e-14 7.38e-14 8.38e-14
pulse_k1280   hrs60 hesdirk2 X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrs60 hesdirk2 C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrs60 hesdirk2 T |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrs60 be       X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrs60 be       C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrs60 be       T |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrsp hesdirk2 X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 1.76e-14
pulse_k1280   hrsp hesdirk2 C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 5.10e-13
pulse_k1280   hrsp hesdirk2 T |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k1280   hrsp be       X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 5.97e-14
pulse_k1280   hrsp be       C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 3.82e-15
pulse_k1280   hrsp be       T |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hr   hesdirk2 X |hr-cen| 4.82e-07 9.52e-07 1.89e-06 3.78e-06 7.56e-06
pulse_k12800  hr   hesdirk2 C |hr-cen| 4.83e-07 9.53e-07 1.89e-06 3.78e-06 7.55e-06
pulse_k12800  hr   hesdirk2 T |hr-cen| 3.79e-06 3.79e-06 3.78e-06 3.78e-06 3.77e-06
pulse_k12800  hr   be       X |hr-cen| 4.82e-07 9.52e-07 1.89e-06 3.78e-06 7.56e-06
pulse_k12800  hr   be       C |hr-cen| 4.79e-07 9.49e-07 1.89e-06 3.78e-06 7.56e-06
pulse_k12800  hr   be       T |hr-cen| 3.76e-06 3.77e-06 3.78e-06 3.78e-06 3.78e-06
pulse_k12800  hrs  hesdirk2 X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs  hesdirk2 C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs  hesdirk2 T |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs  be       X |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs  be       C |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs  be       T |hrs-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs15 hesdirk2 X |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs15 hesdirk2 C |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs15 hesdirk2 T |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs15 be       X |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs15 be       C |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs15 be       T |hrs15-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs60 hesdirk2 X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs60 hesdirk2 C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs60 hesdirk2 T |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs60 be       X |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs60 be       C |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrs60 be       T |hrs60-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrsp hesdirk2 X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrsp hesdirk2 C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrsp hesdirk2 T |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrsp be       X |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrsp be       C |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
pulse_k12800  hrsp be       T |hrsp-cen| 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00
rw_t10        hr   hesdirk2 X |hr-cen| 8.09e-06 1.05e-06 5.90e-05 3.03e-04 8.68e-04
rw_t10        hr   hesdirk2 C |hr-cen| 2.53e-05 1.36e-05 9.20e-05 3.83e-04 8.61e-04
rw_t10        hr   hesdirk2 T |hr-cen| 9.18e-05 7.76e-05 5.76e-05 3.73e-05 2.19e-05
rw_t10        hr   be       X |hr-cen| 8.36e-05 2.74e-04 9.04e-04 2.59e-03 5.17e-03
rw_t10        hr   be       C |hr-cen| 1.51e-03 2.68e-03 4.67e-03 6.04e-03 4.98e-03
rw_t10        hr   be       T |hr-cen| 4.58e-03 2.06e-03 8.48e-04 2.98e-04 7.85e-05
rw_t10        hrs  hesdirk2 X |hrs-cen| 7.56e-07 6.08e-07 4.37e-07 7.82e-07 2.29e-04
rw_t10        hrs  hesdirk2 C |hrs-cen| 2.76e-07 3.77e-08 1.33e-07 8.08e-06 2.08e-04
rw_t10        hrs  hesdirk2 T |hrs-cen| 8.72e-08 3.26e-08 3.92e-07 8.53e-08 4.53e-07
rw_t10        hrs  be       X |hrs-cen| 3.22e-08 2.59e-08 8.25e-07 1.92e-04 6.83e-03
rw_t10        hrs  be       C |hrs-cen| 3.24e-08 1.39e-04 3.25e-06 1.10e-03 6.54e-03
rw_t10        hrs  be       T |hrs-cen| 3.11e-06 1.85e-07 7.60e-07 3.70e-07 1.33e-07
rw_t10        hrs15 hesdirk2 X |hrs15-cen| 7.28e-07 6.20e-07 8.39e-07 1.89e-04 1.61e-03
rw_t10        hrs15 hesdirk2 C |hrs15-cen| 2.77e-07 1.85e-07 2.32e-05 8.36e-04 1.61e-03
rw_t10        hrs15 hesdirk2 T |hrs15-cen| 2.25e-05 4.68e-06 7.24e-07 9.14e-08 1.62e-08
rw_t10        hrs15 be       X |hrs15-cen| 2.36e-08 8.90e-07 2.04e-04 6.94e-03 7.93e-03
rw_t10        hrs15 be       C |hrs15-cen| 5.01e-07 8.72e-05 4.09e-03 1.61e-02 7.66e-03
rw_t10        hrs15 be       T |hrs15-cen| 3.98e-03 1.09e-03 1.77e-04 1.32e-05 2.85e-07
rw_t10        hrs60 hesdirk2 X |hrs60-cen| 7.32e-07 6.65e-07 4.39e-07 2.74e-07 5.78e-06
rw_t10        hrs60 hesdirk2 C |hrs60-cen| 3.22e-07 5.98e-08 9.38e-08 1.62e-07 5.49e-06
rw_t10        hrs60 hesdirk2 T |hrs60-cen| 7.41e-08 2.68e-08 3.73e-07 9.24e-07 1.43e-08
rw_t10        hrs60 be       X |hrs60-cen| 1.81e-07 3.36e-08 3.04e-08 2.94e-07 1.86e-04
rw_t10        hrs60 be       C |hrs60-cen| 1.91e-07 1.40e-04 1.69e-05 3.82e-06 1.70e-04
rw_t10        hrs60 be       T |hrs60-cen| 1.62e-05 4.67e-06 4.57e-08 3.84e-07 1.56e-07
rw_t10        hrsp hesdirk2 X |hrsp-cen| 7.60e-07 6.11e-07 3.96e-07 1.08e-06 2.66e-04
rw_t10        hrsp hesdirk2 C |hrsp-cen| 2.57e-07 5.92e-08 2.65e-07 1.12e-05 2.41e-04
rw_t10        hrsp hesdirk2 T |hrsp-cen| 2.10e-07 3.86e-08 3.32e-07 8.15e-08 4.54e-07
rw_t10        hrsp be       X |hrsp-cen| 3.01e-08 2.54e-08 7.08e-07 1.79e-04 6.69e-03
rw_t10        hrsp be       C |hrsp-cen| 2.64e-08 1.39e-04 2.98e-06 1.07e-03 6.40e-03
rw_t10        hrsp be       T |hrsp-cen| 3.10e-06 4.24e-07 6.74e-07 3.54e-07 8.34e-08
rw_t1000      hr   hesdirk2 X |hr-cen| 8.19e-05 1.64e-04 3.28e-04 6.56e-04
rw_t1000      hr   hesdirk2 C |hr-cen| 8.09e-05 1.65e-04 3.32e-04 6.60e-04
rw_t1000      hr   hesdirk2 T |hr-cen| 3.30e-04 3.28e-04 3.22e-04 3.10e-04
rw_t1000      hr   be       X |hr-cen| 8.33e-05 1.67e-04 3.33e-04 6.66e-04
rw_t1000      hr   be       C |hr-cen| 8.41e-05 1.68e-04 3.34e-04 6.68e-04
rw_t1000      hr   be       T |hr-cen| 3.34e-04 3.33e-04 3.31e-04 3.27e-04
rw_t1000      hrs  hesdirk2 X |hrs-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      hrs  hesdirk2 C |hrs-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      hrs  hesdirk2 T |hrs-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      hrs  be       X |hrs-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      hrs  be       C |hrs-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      hrs  be       T |hrs-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      hrs15 hesdirk2 X |hrs15-cen| 2.24e-07 2.24e-07 2.24e-07 2.23e-07
rw_t1000      hrs15 hesdirk2 C |hrs15-cen| 3.03e-10 1.90e-08 2.34e-09 3.59e-07
rw_t1000      hrs15 hesdirk2 T |hrs15-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      hrs15 be       X |hrs15-cen| 5.35e-07 4.54e-07 1.53e-07 2.08e-08
rw_t1000      hrs15 be       C |hrs15-cen| 1.27e-09 4.97e-10 6.13e-09 1.69e-09
rw_t1000      hrs15 be       T |hrs15-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      hrs60 hesdirk2 X |hrs60-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      hrs60 hesdirk2 C |hrs60-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      hrs60 hesdirk2 T |hrs60-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      hrs60 be       X |hrs60-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      hrs60 be       C |hrs60-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      hrs60 be       T |hrs60-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
rw_t1000      hrsp hesdirk2 X |hrsp-cen| 2.24e-07 2.24e-07 2.24e-07 2.24e-07
rw_t1000      hrsp hesdirk2 C |hrsp-cen| 3.03e-10 1.90e-08 2.34e-09 3.60e-07
rw_t1000      hrsp hesdirk2 T |hrsp-cen| 5.73e-07 2.18e-07 9.88e-08 4.58e-08
rw_t1000      hrsp be       X |hrsp-cen| 5.35e-07 4.54e-07 1.53e-07 2.10e-08
rw_t1000      hrsp be       C |hrsp-cen| 1.27e-09 4.97e-10 6.13e-09 1.75e-09
rw_t1000      hrsp be       T |hrsp-cen| 2.85e-09 2.08e-07 2.07e-07 8.59e-08
```

## RESULTS/gates.txt (raw)

```
g3|/resnick/groups/carnegie_poc/jingze/steep_1010/run/gates/g3_hrsp_c10|/resnick/groups/carnegie_poc/jingze/steep_1010/run/gates/fix_marsh2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm rad_m1/implicit_cfl=10 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak: line 3: 3791947 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1|/resnick/groups/carnegie_poc/jingze/steep_1010/run/gates/g1_hrsp_1e2|/resnick/groups/carnegie_poc/jingze/steep_1010/run/gates/fix_atm2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm rad_m1/implicit_cfl=1e2 time/tlim=2000 output1/dt=2000 output2/dt=2000: line 3: 3791937 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1|/resnick/groups/carnegie_poc/jingze/steep_1010/run/gates/g1_hrsp_1e4|/resnick/groups/carnegie_poc/jingze/steep_1010/run/gates/fix_atm2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm rad_m1/implicit_cfl=1e4 time/tlim=2000 output1/dt=2000 output2/dt=2000: line 3: 3791939 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g3|/resnick/groups/carnegie_poc/jingze/steep_1010/run/gates/g3_hrsp_c1|/resnick/groups/carnegie_poc/jingze/steep_1010/run/gates/fix_marsh2d.athinput|rad_m1/closure=vet_sc rad_m1/vet_tensor=full rad_m1/implicit_flux_faces=all rad_m1/implicit_flux_beam=halfrange rad_m1/implicit_flux=blend rad_m1/implicit_blend=idort_f rad_m1/implicit_blend_tau0=1.0 rad_m1/implicit_blend_alpha=1.0 rad_m1/implicit_blend_flo=0.3 rad_m1/implicit_blend_fhi=0.6 rad_m1/implicit_blend_xthin_mode=steep rad_m1/implicit_blend_xthin=30 rad_m1/implicit_hr_recon=plm rad_m1/implicit_cfl=1 rad_m1/implicit_bc_x1min=marshak rad_m1/implicit_ebath_x1min=1.0e-10 rad_m1/implicit_bc_x1max=marshak: line 3: 3791946 Segmentation fault      (core dumped) OMP_NUM_THREADS=1 nice -n 10 $ATHENA_CPU -i $IN -d $R rad_m1/transport=implicit time/cfl_number=1e6 $K > run.log 2>&1
g1_hrsp_1e4 rc=139  | HOPF: max|dE/E| 1.540e-01 top5 1.540e-01 L1 3.847e-03 (E_top 2.0015 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hrsp_1e2 rc=139  | HOPF: max|dE/E| 1.540e-01 top5 1.540e-01 L1 3.847e-03 (E_top 2.0015 exact 1.7344, tau_top-cell 5.09e-04) 
g3_hrsp_c10 rc=139  | FAIL T6 marshak: L1(E)=1.0000 L1(material)=1.0000 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hrsp_c1 rc=139  | FAIL T6 marshak: L1(E)=1.0000 L1(material)=1.0000 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g1_hrs15_1e4 rc=0 Picard iterations mean=1.741176e+01 | HOPF: max|dE/E| 3.311e-02 top5 7.078e-03 L1 2.335e-03 (E_top 1.7385 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hrs_1e4 rc=0 Picard iterations mean=2.303846e+01 | HOPF: max|dE/E| 2.147e-02 top5 5.162e-03 L1 1.091e-03 (E_top 1.7354 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hrs60_1e4 rc=0 Picard iterations mean=2.201961e+01 | HOPF: max|dE/E| 1.632e-02 top5 4.491e-03 L1 7.240e-04 (E_top 1.7344 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hr_1e4 rc=0 Picard iterations mean=2.150980e+01 | HOPF: max|dE/E| 1.774e-02 top5 5.409e-03 L1 1.946e-03 (E_top 1.7359 exact 1.7344, tau_top-cell 5.09e-04) 
g5_hrsp_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295065e-04 2D=5.208333e-04 ratio=1.016652 (target 1.000000 +- 0.02) 
g5_hrs15_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_hrs15_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.780463e-03 2D=5.208333e-03 ratio=0.533849 (target 1.000000 +- 0.02) 
g5_hrs60_k12.8 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.941056e-02 2D=5.208333e-02 ratio=0.564683 (target 1.000000 +- 0.02) 
g5_hrs_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295065e-04 2D=5.208333e-04 ratio=1.016652 (target 1.000000 +- 0.02) 
g5_hr_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208163e-05 2D=5.208333e-05 ratio=0.999967 (target 1.000000 +- 0.02) 
g5_hrs_k12.8 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.941056e-02 2D=5.208333e-02 ratio=0.564683 (target 1.000000 +- 0.02) 
g5_hrs60_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.780463e-03 2D=5.208333e-03 ratio=0.533849 (target 1.000000 +- 0.02) 
g5_hrs_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.780463e-03 2D=5.208333e-03 ratio=0.533849 (target 1.000000 +- 0.02) 
g5_hrs15_k12.8 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.941056e-02 2D=5.208333e-02 ratio=0.564683 (target 1.000000 +- 0.02) 
g5_hrs_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_hrs60_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_hr_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.756269e-03 2D=5.208333e-03 ratio=0.529204 (target 1.000000 +- 0.02) 
g5_hrsp_k128 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.780463e-03 2D=5.208333e-03 ratio=0.533849 (target 1.000000 +- 0.02) 
g5_hr_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295020e-04 2D=5.208333e-04 ratio=1.016644 (target 1.000000 +- 0.02) 
g5_hrsp_k12800 rc=0 Picard iterations mean=1.023529e+00 | PASS T3 pulse: d(sigma^2)/dt=5.208134e-05 2D=5.208333e-05 ratio=0.999962 (target 1.000000 +- 0.02) 
g5_hrs15_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295065e-04 2D=5.208333e-04 ratio=1.016652 (target 1.000000 +- 0.02) 
g5_hrs60_k1280 rc=0 Picard iterations mean=1.827451e+00 | PASS T3 pulse: d(sigma^2)/dt=5.295065e-04 2D=5.208333e-04 ratio=1.016652 (target 1.000000 +- 0.02) 
g5_hr_k12.8 rc=0 Picard iterations mean=2.000000e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.949960e-02 2D=5.208333e-02 ratio=0.566392 (target 1.000000 +- 0.02) 
g5_hrsp_k12.8 rc=0 Picard iterations mean=2.223529e+00 | FAIL T3 pulse: d(sigma^2)/dt=2.872951e-02 2D=5.208333e-02 ratio=0.551607 (target 1.000000 +- 0.02) 
g1_hr_1e2 rc=0 Picard iterations mean=2.060168e+00 | HOPF: max|dE/E| 1.751e-02 top5 5.118e-03 L1 1.557e-03 (E_top 1.7354 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hrs60_1e2 rc=0 Picard iterations mean=2.062122e+00 | HOPF: max|dE/E| 1.683e-02 top5 4.521e-03 L1 6.613e-04 (E_top 1.7344 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hrs15_1e2 rc=0 Picard iterations mean=2.052745e+00 | HOPF: max|dE/E| 2.109e-02 top5 5.324e-03 L1 1.155e-03 (E_top 1.7357 exact 1.7344, tau_top-cell 5.09e-04) 
g1_hrs_1e2 rc=0 Picard iterations mean=2.048838e+00 | HOPF: max|dE/E| 1.670e-02 top5 4.508e-03 L1 6.670e-04 (E_top 1.7344 exact 1.7344, tau_top-cell 5.09e-04) 
g3_hrs_c10 rc=0 Picard iterations mean=2.526010e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hr_c10 rc=0 Picard iterations mean=2.526977e+00 | PASS T6 marshak: L1(E)=0.0140 L1(material)=0.0140 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hrs15_c10 rc=0 Picard iterations mean=2.526010e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hrs60_c10 rc=0 Picard iterations mean=2.526010e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hr_c1 rc=0 Picard iterations mean=2.006111e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hrs_c1 rc=0 Picard iterations mean=2.006111e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hrs60_c1 rc=0 Picard iterations mean=2.006111e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
g3_hrs15_c1 rc=0 Picard iterations mean=2.006111e+00 | PASS T6 marshak: L1(E)=0.0138 L1(material)=0.0138 (<= 0.02) reference=table /resnick/groups/carnegie_poc/jingze/steep_1010/repo/tests_m1/runs_3a/t6_ref_sn.txt 
```

## RESULTS/beams.txt (raw)

```
shd3b hr rc=0 fatal=0 302.23 s wall | Picard iterations mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=2.062000e+01 max=5.100000e+01 | breakdowns=0.000000e+00 | 
shd3b hrs rc=0 fatal=0 285.84 s wall | Picard iterations mean=2.970000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.700000e+01 max=3.000000e+01 | breakdowns=0.000000e+00 | 
shd3b hrs15 rc=0 fatal=0 285.38 s wall | Picard iterations mean=2.970000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.700000e+01 max=3.000000e+01 | breakdowns=0.000000e+00 | 
shd3b hrs60 rc=0 fatal=0 347.04 s wall | Picard iterations mean=2.970000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.700337e+01 max=3.200000e+01 | breakdowns=0.000000e+00 | 
shd3b hrsp rc=0 fatal=0 611.72 s wall | Picard iterations mean=1.310000e+01 max=2.200000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.478626e+01 max=3.600000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=1.242000e+03
shd3b hrp rc=0 fatal=0 577.55 s wall | Picard iterations mean=1.324000e+01 max=2.200000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.379683e+01 max=3.800000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=1.257000e+03
shd3b cen rc=0 fatal=0 526.86 s wall | Picard iterations mean=3.030000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=7.470627e+01 max=2.000000e+02 | breakdowns=0.000000e+00 | 
cyl hr rc=0 fatal=0 118.11 s wall | Picard iterations mean=5.000000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.678400e+01 max=1.620000e+02 | breakdowns=3.000000e+00 | 
cyl hrs rc=0 fatal=0 105.36 s wall | Picard iterations mean=4.970000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.126157e+01 max=9.600000e+01 | breakdowns=3.000000e+00 | 
cyl hrs15 rc=0 fatal=0 106.60 s wall | Picard iterations mean=4.970000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.125956e+01 max=9.500000e+01 | breakdowns=3.000000e+00 | 
cyl hrs60 rc=0 fatal=0 106.19 s wall | Picard iterations mean=4.970000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.125755e+01 max=9.500000e+01 | breakdowns=3.000000e+00 | 
cyl hrsp rc=0 fatal=0 136.86 s wall | Picard iterations mean=1.300000e+01 max=2.200000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=6.184615e+00 max=1.700000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=1.196000e+03
cyl hrp rc=0 fatal=0 340.63 s wall | Picard iterations mean=3.768000e+01 max=4.800000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=8.670913e+00 max=3.000000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=0.000000e+00
cyl cen rc=0 fatal=0 116.65 s wall | Picard iterations mean=2.180000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=4.687156e+01 max=2.000000e+02 | breakdowns=2.000000e+00 | 
xb20 hr rc=0 fatal=0 311.95 s wall | Picard iterations mean=4.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=3.749750e+01 max=8.900000e+01 | breakdowns=0.000000e+00 | 
xb20 hrs rc=0 fatal=0 269.33 s wall | Picard iterations mean=5.580000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.449821e+01 max=1.800000e+02 | breakdowns=3.000000e+00 | 
xb20 hrs15 rc=0 fatal=0 243.02 s wall | Picard iterations mean=5.560000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.458273e+01 max=1.740000e+02 | breakdowns=3.000000e+00 | 
xb20 hrs60 rc=0 fatal=0 232.70 s wall | Picard iterations mean=5.320000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.484962e+01 max=4.100000e+01 | breakdowns=0.000000e+00 | 
xb20 hrsp rc=1 fatal=1 87.20 s wall |  |  |  | 
xb20 hrp rc=0 fatal=0 394.45 s wall | Picard iterations mean=1.687000e+01 max=2.300000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=9.515116e+00 max=4.700000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=1.654000e+03
xb20 cen rc=0 fatal=0 482.10 s wall | Picard iterations mean=3.720000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=9.485484e+01 max=2.000000e+02 | breakdowns=0.000000e+00 | 
ba0 hr rc=0 fatal=0 360.63 s wall | Picard iterations mean=4.710000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=3.987261e+01 max=1.050000e+02 | breakdowns=0.000000e+00 | 
ba0 hrs rc=0 fatal=0 266.92 s wall | Picard iterations mean=4.720000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.929873e+01 max=5.200000e+01 | breakdowns=0.000000e+00 | 
ba0 hrs15 rc=0 fatal=0 271.85 s wall | Picard iterations mean=4.720000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.921822e+01 max=5.800000e+01 | breakdowns=0.000000e+00 | 
ba0 hrs60 rc=0 fatal=0 253.54 s wall | Picard iterations mean=4.700000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.921702e+01 max=5.300000e+01 | breakdowns=0.000000e+00 | 
ba0 hrsp rc=1 fatal=1 223.84 s wall |  |  |  | 
ba0 hrp rc=0 fatal=0 538.32 s wall | Picard iterations mean=2.039000e+01 max=2.800000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.280726e+01 max=7.800000e+01 | breakdowns=0.000000e+00 | plm positivity fallbacks=2.016000e+03
ba0 cen rc=0 fatal=0 633.23 s wall | Picard iterations mean=4.710000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.057495e+02 max=2.000000e+02 | breakdowns=0.000000e+00 | 
ba20 hr rc=0 fatal=0 409.60 s wall | Picard iterations mean=5.010000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=4.549501e+01 max=2.000000e+02 | breakdowns=0.000000e+00 | 
ba20 hrs rc=0 fatal=0 324.49 s wall | Picard iterations mean=5.830000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=2.127273e+01 max=7.100000e+01 | breakdowns=0.000000e+00 | 
ba20 hrs15 rc=0 fatal=0 329.42 s wall | Picard iterations mean=5.720000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=2.202448e+01 max=2.000000e+02 | breakdowns=1.000000e+00 | 
ba20 hrs60 rc=0 fatal=0 289.94 s wall | Picard iterations mean=5.810000e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=2.127022e+01 max=7.100000e+01 | breakdowns=0.000000e+00 | 
ba20 hrsp rc=1 fatal=1 423.69 s wall |  |  |  | 
ba20 hrp rc=0 fatal=0 1075.08 s wall | Picard iterations mean=3.332000e+01 max=4.600000e+01 NON-CONVERGED=0.000000e+00 | inner iterations mean=1.898559e+01 max=1.290000e+02 | breakdowns=0.000000e+00 | plm positivity fallbacks=3.304000e+03
ba20 cen rc=0 fatal=0 468.13 s wall | Picard iterations mean=4.230000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 | inner iterations mean=7.967612e+01 max=2.000000e+02 | breakdowns=0.000000e+00 | 
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/cyl/hr t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.377  max|E/J-1| 0.571  <|F|/|cH|> 0.673  median dangle 4.6 deg  p90 11.0  <f_M1> 0.502 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.673  median dangle 4.6  face max |F|/cE 0.65  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/cyl/hrs t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.377  max|E/J-1| 0.582  <|F|/|cH|> 0.673  median dangle 4.6 deg  p90 11.0  <f_M1> 0.502 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.673  median dangle 4.6  face max |F|/cE 0.65  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/cyl/hrsp t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.447  max|E/J-1| 0.995  <|F|/|cH|> 0.606  median dangle 4.8 deg  p90 11.9  <f_M1> 0.503 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 0.606  median dangle 4.8  face max |F|/cE 2.11  frac(|F_face| > cE) 0.000
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/cyl/cen t=2.421e-01
  x in 0.3-0.95: <|E/J-1|> 0.414  max|E/J-1| 0.635  <|F|/|cH|> 22.300  median dangle 87.8 deg  p90 143.2  <f_M1> 18.064 <f_ex> 0.462
  cell flux (derived, clipped): <|F|/|cH|> 1.262  median dangle 87.8  face max |F|/cE 130.06  frac(|F_face| > cE) 0.994
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/xb20/hr t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 2.930  sum E/sum J 2.821  |F|/|cH| median 2.508  flux-weighted <|F|>/<|cH|> 2.426
  beam: direction error J-weighted mean 4.8 deg  median 3.8  p90 12.6   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.94
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.1008  max E/Jmax 1.0909
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.811
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.158  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.396
    peaks y: E [0.965 1.035]  J [1.004]   E max/J max 2.047  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.763
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 4.177  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.365
    peaks y: E [0.699 0.879 1.004 1.121 1.301]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.698  width(>50% max) E 0.391 J 0.281
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/xb20/hrs t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 3.019  sum E/sum J 2.898  |F|/|cH| median 2.574  flux-weighted <|F|>/<|cH|> 2.489
  beam: direction error J-weighted mean 4.7 deg  median 3.8  p90 12.7   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.94
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.0805  max E/Jmax 1.1033
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.714
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 3.803  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.331
    peaks y: E [0.965 1.035]  J [1.004]   E max/J max 2.095  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.633
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 3.774  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.285
    peaks y: E [0.879 1.004 1.121]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.778  width(>50% max) E 0.391 J 0.281
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/xb20/hrsp t=0.000e+00  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 0.000  sum E/sum J 0.000  |F|/|cH| median 0.000  flux-weighted <|F|>/<|cH|> 0.000
  beam: direction error J-weighted mean 14.4 deg  median 19.4  p90 26.3   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.00
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.0000  max E/Jmax 0.0000
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.771
  x=0.70 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.891
    peaks y: E []  J [1.004]   E max/J max 0.000  width(>50% max) E 2.000 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.554
  x=0.90 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.773
    peaks y: E []  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 0.000  width(>50% max) E 2.000 J 0.281
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/xb20/cen t=2.421e-01  [xb20, x in 0.45..0.95, beam = J > 0.05 Jmax: 2106 cells]
  beam: E/J median 1.788  sum E/sum J 1.802  |F|/|cH| median 884.853  flux-weighted <|F|>/<|cH|> 1178.540
  beam: direction error J-weighted mean 64.5 deg  median 64.1  p90 160.4   face |F|/cE > 1: frac 1.0000 (beam) 0.9999 (all)  max 3656.19
  dark (J < 1e-3 Jmax, 14108 cells): <E>/Jmax 0.2839  max E/Jmax 0.7855
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.378
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 5.514  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.693
    peaks y: E [0.074 0.965 1.035 1.926]  J [1.004]   E max/J max 1.305  width(>50% max) E 0.188 J 0.109
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.966
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 3.355  E-centroid 1.0000 (exact 1.0000)  E outside the J>5% band / total 0.449
    peaks y: E [1.004]  J [0.887 0.910 0.926 1.004 1.098 1.113]   E max/J max 2.652  width(>50% max) E 0.281 J 0.281
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/ba0/hr t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.073  sum E/sum J 5.599  |F|/|cH| median 7.574  flux-weighted <|F|>/<|cH|> 5.862
  beam: direction error J-weighted mean 5.1 deg  median 1.0  p90 4.9   face |F|/cE > 1: frac 0.0096 (beam) 0.0018 (all)  max 78.88
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.1141  max E/Jmax 6.6123
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.942
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 11.140  E-centroid 0.9806 (exact 0.9732)  E outside the J>5% band / total 0.452
    peaks y: E [0.965]  J [0.980 0.996]   E max/J max 5.894  width(>50% max) E 0.070 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.040
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 12.732  E-centroid 0.9883 (exact 0.9633)  E outside the J>5% band / total 0.516
    peaks y: E [0.957]  J [0.988]   E max/J max 5.715  width(>50% max) E 0.078 J 0.062
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/ba0/hrs t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.066  sum E/sum J 5.597  |F|/|cH| median 7.562  flux-weighted <|F|>/<|cH|> 5.858
  beam: direction error J-weighted mean 5.1 deg  median 1.0  p90 4.9   face |F|/cE > 1: frac 0.0096 (beam) 0.0016 (all)  max 78.88
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.0595  max E/Jmax 10.2512
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.710
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 8.356  E-centroid 0.9605 (exact 0.9732)  E outside the J>5% band / total 0.289
    peaks y: E [0.965]  J [0.980 0.996]   E max/J max 5.892  width(>50% max) E 0.070 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 0.807
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 9.519  E-centroid 0.9355 (exact 0.9633)  E outside the J>5% band / total 0.370
    peaks y: E [0.957]  J [0.988]   E max/J max 5.701  width(>50% max) E 0.086 J 0.062
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/ba0/hrsp t=0.000e+00  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 0.000  sum E/sum J 0.000  |F|/|cH| median 0.000  flux-weighted <|F|>/<|cH|> 0.000
  beam: direction error J-weighted mean 2.8 deg  median 3.6  p90 7.6   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.00
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.0000  max E/Jmax 0.0000
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.922
  x=0.70 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 0.9732)  E outside the J>5% band / total 0.965
    peaks y: E []  J [0.980 0.996]   E max/J max 0.000  width(>50% max) E 2.000 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.876
  x=0.90 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 0.9633)  E outside the J>5% band / total 0.945
    peaks y: E []  J [0.988]   E max/J max 0.000  width(>50% max) E 2.000 J 0.062
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/ba0/cen t=2.421e-01  [ba0, x in 0.45..0.95, beam = J > 0.05 Jmax: 730 cells]
  beam: E/J median 7.967  sum E/sum J 6.898  |F|/|cH| median 2825.598  flux-weighted <|F|>/<|cH|> 3257.458
  beam: direction error J-weighted mean 55.1 deg  median 85.9  p90 169.3   face |F|/cE > 1: frac 0.9521 (beam) 0.9890 (all)  max 9542.13
  dark (J < 1e-3 Jmax, 15608 cells): <E>/Jmax 0.3118  max E/Jmax 5.9405
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.542
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 35.426  E-centroid 0.9595 (exact 0.9732)  E outside the J>5% band / total 0.778
    peaks y: E [0.973]  J [0.980 0.996]   E max/J max 7.487  width(>50% max) E 0.062 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.309
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 23.945  E-centroid 0.9979 (exact 0.9633)  E outside the J>5% band / total 0.662
    peaks y: E [0.965]  J [0.988]   E max/J max 7.892  width(>50% max) E 0.070 J 0.062
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/ba20/hr t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 2.425  sum E/sum J 2.335  |F|/|cH| median 2.145  flux-weighted <|F|>/<|cH|> 1.999
  beam: direction error J-weighted mean 12.2 deg  median 7.0  p90 31.3   face |F|/cE > 1: frac 0.0113 (beam) 0.0015 (all)  max 78.23
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0784  max E/Jmax 3.0565
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.024
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 6.069  E-centroid 0.8618 (exact 0.8959)  E outside the J>5% band / total 0.489
    peaks y: E [0.895]  J [0.918]   E max/J max 2.855  width(>50% max) E 0.078 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.256
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 7.469  E-centroid 0.9191 (exact 0.9590)  E outside the J>5% band / total 0.592
    peaks y: E [0.918]  J [0.980]   E max/J max 2.468  width(>50% max) E 0.156 J 0.070
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/ba20/hrs t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 2.522  sum E/sum J 2.432  |F|/|cH| median 2.203  flux-weighted <|F|>/<|cH|> 2.075
  beam: direction error J-weighted mean 12.5 deg  median 7.3  p90 31.6   face |F|/cE > 1: frac 0.0113 (beam) 0.0012 (all)  max 78.24
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0579  max E/Jmax 2.7333
  x=0.70 shape: L1 |E/sumE - J/sumJ| 0.831
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 5.163  E-centroid 0.8553 (exact 0.8959)  E outside the J>5% band / total 0.378
    peaks y: E [0.895]  J [0.918]   E max/J max 3.034  width(>50% max) E 0.078 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.168
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 6.852  E-centroid 0.8888 (exact 0.9590)  E outside the J>5% band / total 0.515
    peaks y: E [0.910]  J [0.980]   E max/J max 2.866  width(>50% max) E 0.156 J 0.070
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/ba20/hrsp t=0.000e+00  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 0.000  sum E/sum J 0.000  |F|/|cH| median 0.000  flux-weighted <|F|>/<|cH|> 0.000
  beam: direction error J-weighted mean 18.4 deg  median 17.3  p90 22.1   face |F|/cE > 1: frac 0.0000 (beam) 0.0000 (all)  max 0.00
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.0000  max E/Jmax 0.0000
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.914
  x=0.70 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 0.8959)  E outside the J>5% band / total 0.957
    peaks y: E []  J [0.918]   E max/J max 0.000  width(>50% max) E 2.000 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.868
  x=0.90 profile: overlap sum min(E,J)/sum J 0.000  L1 |E-J|/sum J 1.000  E-centroid 1.0000 (exact 0.9590)  E outside the J>5% band / total 0.938
    peaks y: E []  J [0.980]   E max/J max 0.000  width(>50% max) E 2.000 J 0.070
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/ba20/cen t=2.421e-01  [ba20, x in 0.45..0.95, beam = J > 0.05 Jmax: 797 cells]
  beam: E/J median 1.523  sum E/sum J 1.527  |F|/|cH| median 1129.388  flux-weighted <|F|>/<|cH|> 1382.519
  beam: direction error J-weighted mean 95.8 deg  median 74.2  p90 164.3   face |F|/cE > 1: frac 0.9435 (beam) 0.9883 (all)  max 16004.06
  dark (J < 1e-3 Jmax, 15535 cells): <E>/Jmax 0.3362  max E/Jmax 12.1155
  x=0.70 shape: L1 |E/sumE - J/sumJ| 1.875
  x=0.70 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 41.456  E-centroid 1.2620 (exact 0.8959)  E outside the J>5% band / total 0.937
    peaks y: E [0.863 1.504 1.551 1.637]  J [0.918]   E max/J max 6.334  width(>50% max) E 0.203 J 0.047
  x=0.90 shape: L1 |E/sumE - J/sumJ| 1.715
  x=0.90 profile: overlap sum min(E,J)/sum J 1.000  L1 |E-J|/sum J 19.569  E-centroid 0.8809 (exact 0.9590)  E outside the J>5% band / total 0.857
    peaks y: E [0.121 0.309 0.434 0.895]  J [0.980]   E max/J max 2.421  width(>50% max) E 0.266 J 0.070
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/shd3b/hr shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.232  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.593
  x=0.85  E depth 0.000 edge 0.346  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.608
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/shd3b/hrs shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.232  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.593
  x=0.85  E depth 0.000 edge 0.346  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.608
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/shd3b/hrsp shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.000 edge 0.350  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.651
  x=0.85  E depth 0.000 edge 0.503  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.495
/resnick/groups/carnegie_poc/jingze/steep_1010/run/beams/shd3b/cen shd3b.m1_vsc.00001.bin minJfs 2.600e-13
  x=0.65  E depth 0.284 edge 0.017  Jfs depth 0.000 edge 0.205  ex depth 0.000 edge 0.184  Jfs/ex lit 0.626 E/ex lit 0.144
  x=0.85  E depth 0.208 edge 0.011  Jfs depth 0.000 edge 0.328  ex depth 0.000 edge 0.295  Jfs/ex lit 0.528 E/ex lit 0.288
```
