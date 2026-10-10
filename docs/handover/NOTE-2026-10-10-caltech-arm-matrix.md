# NOTE Caltech -> viper: xthinfix-1009 ARM MATRIX (413b38af), 10-10

Answers TASK-2026-10-10-caltech-arm-matrix.md. Run trees stay on Caltech: /resnick/groups/carnegie_poc/jingze/arms_1010
(run/, logs/; REPO = git archive of 413b38af in repo/).

STATUS: first push = binaries + GPU point + early gate finding. CPU battery running (order 10 x 32 cores 4308649-58,
gates 4308659 DONE, beams 4308660, ana 4308661); Debug build + gdb backtrace of a crashing plm gate in 4309595.

EARLY FINDING: every plm arm segfaults in the 2-D gates. G1 (atm2d, cfl 1e2 and 1e4) and G3 (marsh2d, cfl 1 and 10) give
rc=139 for sq, aq and knp (12 of 12 runs); hr and kn03 give rc 0. G5 (pulse2d) runs for all arms. The backtrace follows.

## Binaries

- ATHENA_CPU xthinfix-1009 413b38af, built-in pgens, MPI, Release: `athena_cpu_built_in_pgens_413b38af` md5 `979f0438dab6e2dd29d2dfffc6defdcc`
- ATHENA_GPU 413b38af he_star_m1 H200 nofma (CUDA 12.9, HOPPER90, hpcx, host -ffp-contract=off): `athena_gpu_he_star_m1_413b38af_nofma` md5 `12ce0f6cc19d1f9db258ad74c5c5b5cb`
- ATHENA_BASE rt-integration 5b304cf9 he_star_m1 H200 nofma (reused): `athena_gpu_he_star_m1_5b304cf9_nofma` md5 `a9724eea21f363dc55ef0bc664ab3cf4`

## GPU point (job 4309030, 2 H200; smoke 4308662 sq nlim 3: rc 0, NON-CONVERGED 0, Picard mean 21.3)

Adapted only: #SBATCH header, srun `--mpi=pmix -n 2 -c 8 --cpu-bind=cores gpuwrap.sh`; python with h5py. All 7 arms rc=0 fatal=0.

RESULTS_gpu_arms.txt (raw):
```
base_1    0.463 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
all_1     0.465 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
sq_1      1.901 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.656667e+01', '5.300000e+01', '0.000000e+00')  inner mean 2.947029e+00  plm pos-fallbacks 0.000000e+00  FATAL 0
aq_1      1.825 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.590000e+01', '5.100000e+01', '0.000000e+00')  inner mean 2.946260e+00  plm pos-fallbacks 0.000000e+00  FATAL 0
kn03_1    0.467 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.266667e+00', '8.000000e+00', '0.000000e+00')  inner mean 3.367188e+00  plm pos-fallbacks -  FATAL 0
knp_1     1.752 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.530000e+01', '5.100000e+01', '0.000000e+00')  inner mean 2.902134e+00  plm pos-fallbacks 0.000000e+00  FATAL 0
base_2    0.467 s/cycle (cycles 5-30)  Picard mean/max/NONCONV ('4.333333e+00', '9.000000e+00', '0.000000e+00')  inner mean 3.330769e+00  plm pos-fallbacks -  FATAL 0
-- bin data base_1 vs all_1 (after the parameter header <par_end>), and repeats
base_1 vs all_1: 4 bin files, 0 differ in data
base_1 vs base_2: 4 bin files, 0 differ in data
-- all vs steep / steep+plm by radius (max over angles |a-b|/max|b|), last dump
== all_1 vs sq_1
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 1.3e-06 1.2e-02 6.0e-03 1.3e-04 5.9e-04
agcar3d.hydro_w.00001.bin    velx       9.8e-08 5.0e-08 5.4e-03 3.9e+00 3.3e-02 2.6e-02 1.0e-01
agcar3d.hydro_w.00001.bin    vely       3.8e-06 5.1e-08 2.8e-04 3.6e-01 9.6e-01 1.1e+00 1.0e+00
agcar3d.hydro_w.00001.bin    velz       5.3e-06 2.2e-08 2.8e-04 3.3e-01 9.7e-01 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 4.1e-05 7.6e-02 6.2e-03 8.2e-02 9.5e-02
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 1.7e-04 7.6e-02 9.1e-03 2.7e-02 7.5e-02
agcar3d.m1.00001.bin         m1_f1      8.0e-08 2.0e-08 8.5e-05 7.1e-04 1.9e-03 1.3e-02 5.7e-02
agcar3d.m1.00001.bin         m1_f2      5.4e-06 1.8e-08 2.7e-04 2.1e-01 3.8e-01 1.0e+00 1.0e+00
agcar3d.m1.00001.bin         m1_f3      7.9e-06 2.8e-08 2.8e-04 2.0e-01 3.2e-01 9.5e-01 1.1e+00
== all_1 vs aq_1
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 1.2e-06 1.3e-02 1.5e-03 2.1e-04 1.2e-03
agcar3d.hydro_w.00001.bin    velx       9.8e-08 5.0e-08 5.0e-03 3.1e+00 1.0e-02 3.0e-02 5.1e-02
agcar3d.hydro_w.00001.bin    vely       3.7e-06 3.8e-08 2.8e-04 3.4e-01 1.0e+00 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    velz       3.6e-06 7.0e-09 2.8e-04 3.3e-01 8.7e-01 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 3.8e-05 8.2e-02 1.4e-03 8.1e-02 9.5e-02
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 1.6e-04 8.1e-02 9.4e-03 2.7e-02 7.4e-02
agcar3d.m1.00001.bin         m1_f1      8.0e-08 3.7e-09 8.0e-05 7.1e-04 2.1e-03 1.4e-02 5.6e-02
agcar3d.m1.00001.bin         m1_f2      5.4e-06 1.8e-08 2.7e-04 2.3e-01 4.9e-01 1.0e+00 1.0e+00
agcar3d.m1.00001.bin         m1_f3      5.4e-06 2.8e-08 2.8e-04 2.2e-01 2.7e-01 1.0e+00 1.1e+00
== all_1 vs kn03_1
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 2.1e-07 1.3e-03 7.6e-04 2.0e-07 5.5e-07
agcar3d.hydro_w.00001.bin    velx       9.8e-08 3.5e-08 6.3e-04 9.2e-03 5.9e-03 5.4e-05 3.6e-05
agcar3d.hydro_w.00001.bin    vely       3.0e-06 3.8e-08 1.4e-05 2.2e-01 3.3e-01 5.3e-03 5.9e-03
agcar3d.hydro_w.00001.bin    velz       2.9e-06 4.4e-08 1.3e-05 1.6e-01 3.1e-01 8.5e-03 9.4e-03
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 4.9e-06 1.3e-03 7.6e-04 3.1e-05 2.9e-05
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 1.9e-05 6.0e-04 2.5e-05 1.5e-05 9.7e-06
agcar3d.m1.00001.bin         m1_f1      7.9e-08 2.0e-08 9.2e-06 1.5e-05 1.3e-05 1.0e-05 8.2e-06
agcar3d.m1.00001.bin         m1_f2      4.4e-06 1.8e-08 1.0e-05 2.1e-02 1.2e-03 6.7e-04 2.3e-04
agcar3d.m1.00001.bin         m1_f3      4.5e-06 2.8e-08 1.1e-05 2.2e-02 7.7e-04 4.8e-04 3.2e-04
== all_1 vs knp_1
bands r/R_ph: 0-0.5 0.5-0.8 0.8-0.95 0.95-1.05 1.05-1.3 1.3-2 2-10
agcar3d.hydro_w.00001.bin    dens       0.0e+00 0.0e+00 1.3e-06 1.3e-02 5.1e-03 1.6e-04 1.1e-03
agcar3d.hydro_w.00001.bin    velx       9.8e-08 5.0e-08 5.6e-03 3.1e+00 2.9e-02 3.0e-02 5.5e-02
agcar3d.hydro_w.00001.bin    vely       3.6e-06 5.1e-08 2.8e-04 3.4e-01 9.5e-01 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    velz       3.3e-06 7.0e-09 2.8e-04 3.1e-01 9.8e-01 1.0e+00 1.0e+00
agcar3d.hydro_w.00001.bin    eint       0.0e+00 0.0e+00 4.3e-05 8.2e-02 5.3e-03 8.3e-02 8.8e-02
agcar3d.m1.00001.bin         m1_e       0.0e+00 0.0e+00 1.7e-04 8.1e-02 9.4e-03 2.7e-02 5.8e-02
agcar3d.m1.00001.bin         m1_f1      8.0e-08 2.0e-08 8.8e-05 6.9e-04 1.8e-03 1.2e-02 4.0e-02
agcar3d.m1.00001.bin         m1_f2      4.9e-06 1.4e-08 2.7e-04 2.0e-01 4.0e-01 1.1e+00 1.1e+00
agcar3d.m1.00001.bin         m1_f3      4.9e-06 2.8e-08 2.7e-04 2.0e-01 2.9e-01 1.0e+00 1.1e+00
```
