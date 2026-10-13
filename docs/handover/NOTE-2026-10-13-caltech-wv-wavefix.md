# NOTE Caltech -> raven (WV): wavefix-1013 IDORT suite — QUEUED, NEW REBUILT at 3972a05d (10-13 13:35 PDT)

## Binaries (plain CPU, same flags as all Caltech CPU builds: gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI on, Release,
PROBLEM built_in_pgens, no Athena_M1_VGD_NMAX)
- OLD 9972a56d: `athena_cpu_omp_9972a56d` md5 cf843e67b88d9d01606d02f4036bf13a
- NEW wavefix-1013 @ **3972a05d** (TASK UPDATE 5c99c030): `athena_cpu_omp_3972a05d` md5 84f8b7e1be477923f15a9813fbf9ede2
- (superseded, not used for results: 41c8b210 build md5 0661dd47cd434a474ae8404f36648b29; the 4 dirs that had run BIT
  with it were deleted and rerun with 3972a05d)

## Run
Per-dir step = the UPDATED run_wavefix.sh step, 4 arms O / a / b (+ implicit_hr_q_face) / z, FULL a and b for their DIFF
dirs (wv_one.sh) with the Caltech fixes (rstonce.py one rst block, both shockref paths,
gtime stand-in, killed runs never cached). Queued behind OD mbband and OT ENTHBE in the same Caltech debug-QoS rounds
(6 nodes per round, shared claim queue). cart/ddf/FIN3_d and cart/ddiff/FIN3_d run in a separate 2-h expansion job
(4456442) at 4 threads. Debug jobs: 4461226 -> 4461227 -> 4461228 -> 4461229 -> 4461230 (more added as needed).
