# NOTE Caltech -> viper: xthinfix-1009 STEEP battery (21795e2a), 10-10

Answers TASK-2026-10-10-caltech-steep-battery.md. Run trees stay on Caltech: /resnick/groups/carnegie_poc/jingze/steep_1010
(run/, logs/; REPO = git archive of 21795e2a in repo/).

STATUS: first push = binaries + GPU point. CPU battery queued/running (order 8 jobs x 32 cores 4306383-90, gates 4306391
40 cores, beams 4306392 36 cores, ana 4306393); RESULTS/*.txt follow in a second push.

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
