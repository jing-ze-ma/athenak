# NOTE Caltech -> raven (WV): wavefix-1013 IDORT suite — QUEUED (10-13 12:30 PDT)

## Binaries (plain CPU, same flags as all Caltech CPU builds: gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI on, Release,
PROBLEM built_in_pgens, no Athena_M1_VGD_NMAX)
- OLD 9972a56d: `athena_cpu_omp_9972a56d` md5 cf843e67b88d9d01606d02f4036bf13a
- NEW wavefix-1013 @ 41c8b210: `athena_cpu_omp_41c8b210` md5 0661dd47cd434a474ae8404f36648b29

## Run
Per-dir step = run_wavefix.sh's step (wv_one.sh) with the Caltech fixes (rstonce.py one rst block, both shockref paths,
gtime stand-in, killed runs never cached). Queued behind OD mbband and OT ENTHBE in the same Caltech debug-QoS rounds
(6 nodes per round, shared claim queue). cart/ddf/FIN3_d and cart/ddiff/FIN3_d run in a separate 2-h expansion job
(4456442) at 4 threads. Debug jobs: 4456444 -> 4456533 -> 4456535 (more added as needed).
