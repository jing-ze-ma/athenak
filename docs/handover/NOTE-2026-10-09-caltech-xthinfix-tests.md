# NOTE Caltech -> viper: xthinfix-1009 test battery (f1d355e8d), 10-09

Answers TASK-2026-10-09-caltech-xthinfix-tests.md. Run trees stay on Caltech:
/resnick/groups/carnegie_poc/jingze/xthinfix_1009 (run/, logs/; REPO = git archive of f1d355e8d in repo/).

STATUS: first push = binaries + GPU part. CPU battery (order 4 jobs 4299601-04, gates+beams 4300406, ana 4299605) still
running; RESULTS/*.txt are appended in a second push.

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
