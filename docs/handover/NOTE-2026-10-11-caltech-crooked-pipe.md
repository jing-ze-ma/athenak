# NOTE Caltech: crooked pipe (TASK-2026-10-11-caltech-crooked-pipe), DONE

picked up by Caltech 2026-10-10 21:38 PDT; results 2026-10-11 10:30 PDT. No Raven pipe NOTE existed at posting time.

## Build
- Fresh `git archive 4c033dcf` (kokkos d8e9af03 copied in), CPU, built-in pgens, MPI + OpenMP, Release, gcc 13.2.0 + openmpi 5.0.1, mpicxx.
- Binary /resnick/home/jingze/ATHENAK/builds/bin/athena_cpu_omp_4c033dcf, **md5 f27540864e3689c21c2625be65072746** (Raven's is 1f797264..., different compiler/MPI so different md5).
- pipe.athinput md5 bd95c5a091261844e51e6f311a0efec4.

## Runs
One expansion node (hpc-19-18), job 4336114: 4 runs in parallel, each one process (srun -n1 --exact) x 16 OpenMP threads.
run_pipe.sh used unchanged except that /usr/bin/time (absent on compute nodes) was replaced by date-based timing and the binary launched with srun --mpi=pmix.
fixkeys.py ran (added time2_vet_sc, time2_vstage, implicit_vimp for REC); no FATAL.

| run | rc | FATAL | NON-CONV | Picard mean/max | NOT ADMISSIBLE (hesdirk2 stage -> BE redo) | cycles | wall | s/step |
|---|---|---|---|---|---|---|---|---|
| REC_c220 | 0 | 0 | 0 | 3.01 / 6 | 30 (BE steps 31 of 998) | 1000 | 1270 s | 1.27 |
| RECbe_c220 | 0 | 0 | 0 | 3.02 / 6 | n/a (be) | 1000 | 1411 s | 1.41 |
| REC_c2200 | 0 | 0 | 0 | 3.61 / 7 | 8 (BE steps 9 of 100) | 100 | 384 s | 3.84 |
| REC_c22 | 0 | 0 | 0 | 3.04 / 4 | 29 (BE steps 30 of 998) | 1000 | 561 s | 0.56 |

- time2_vet_sc=predict fallbacks = 0 and vet clips = 0 in all hesdirk2 runs. The NOT ADMISSIBLE events are all "converged=1 min(E,T) ~ 6e-5..1e-4" (early cycles, front entering the pipe).
- Other fallbacks: positivity fallbacks 57 / 40 / 16 / 56 (REC_c220 / be / c2200 / c22); min E from the solve -6.3e-2 / -7.8e-2 / -2.14 / -5.0e-3. bicgstab line_jacobi fallbacks 3-4 per run.
- (sacct shows the job FAILED only because the in-job pipe.py call had python -I, which drops the module's PYTHONPATH; pipe.py was rerun on the login node.)

## pipe.py summary
```
RECbe_c220 t_end=10.00 T_end=0.938 0.855 0.773 0.689 0.555 t(T>=0.5)=0.10 0.10 0.30 0.70 3.90 max|v|=0.0e+00
REC_c22 t_end=1.00 T_end=0.902 0.775 0.654 0.555 0.418 t(T>=0.5)=0.10 0.10 0.40 0.70 >1.0 max|v|=0.0e+00
REC_c220 t_end=10.00 T_end=0.938 0.855 0.773 0.688 0.554 t(T>=0.5)=0.10 0.10 0.30 0.70 4.00 max|v|=0.0e+00
REC_c2200 t_end=10.00 T_end=0.938 0.856 0.774 0.690 0.556 t(T>=0.5)=0.10 0.20 0.60 0.80 3.71 max|v|=0.0e+00
```
Probes: (0.25,0) (2.75,0) (3.5,1.25) (4.25,0) (6.75,0). Dumps every 0.1, so t(T>=0.5) has 0.1 resolution.
**max|vx| = 0 in all runs** (static gas holds).

## dt dependence (max relative T difference over the 5 probes)
| pair | t | max rel diff | at probe | others |
|---|---|---|---|---|
| REC_c2200 vs REC_c220 | 1 | 2.99e-1 | 5 (6.75,0), front arriving | probes 1-4: 0.6-2.0 % |
| REC_c2200 vs REC_c220 | 3 | 7.6e-3 | 5 | |
| REC_c2200 vs REC_c220 | 10 | 4.0e-3 | 5 | |
| REC_c22 vs REC_c220 | 0.5 | 6.4e-2 | 4 (4.25,0) | probe 5 6.1e-2 |
| REC_c22 vs REC_c220 | 1.0 | 1.5e-2 | 5 | |
| RECbe_c220 vs REC_c220 | 1 / 3 / 10 | 1.1e-2 / 2.4e-3 / 1.6e-3 | | |
Note: dt = 1.0026e-2 x 10^k (implicit_cfl from dx = 7/512), so dumps of c220/c2200 sit at t = 1.003, 3.008 while c22 dumps sit at 0.500, 1.000; the comparison takes the nearest dump (offset <= 0.003).
Reading: the late-time state converges in dt (0.4 % at t = 10 for 10x dt); the dt error sits in the front arrival times: c2200 reaches T = 0.5 at probes 2/3/4/5 at 0.2/0.6/0.8/3.71 vs 0.1/0.3/0.7/4.00 for c220; c22 differs from c220 only at probe 3 (0.4 vs 0.3) at this 0.1 dump resolution.

## Files
probes.txt of every run: docs/handover/suite-pipe-1011/caltech-probes/<run>.probes.txt (t, 5 probe T, max|vx|).
Run dir: /resnick/groups/carnegie_poc/jingze/pipe_1011/out. Cost 0 GPU-h, 26 CPU-h.
