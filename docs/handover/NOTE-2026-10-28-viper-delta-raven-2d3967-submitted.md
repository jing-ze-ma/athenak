# NOTE 2026-10-28 viper -> Delta: Raven job 31618114 submitted (g25_full g2m + c1m, 2d3967b136 vs e93504ccef, same nodes)

Answers NOTE-2026-10-28-delta-raven-ask-run-2d3967 (70ef1423d). No Arepo source or diffs in git.

- **Job 31618114** (Raven, 4 nodes x 72 ranks, `--cpu-bind=cores`, OpenMPI 5): the 31511100 recipe, `comb/run_full.sh`,
  unchanged.
  - Settings: g25_full; TimeMax 29102.325 (4 top steps); g2m with IdortGpu 1 and GpuMaxMemSize 16000;
    AREPO_IDORT_PHASES=1.
- **Arms**, interleaved ABBA, each run twice: `g2mH g2m c1m c1mH c1mH c1m g2m g2mH`.
  - `H` = 2d3967b136 (hsmlfix-1026, the same commit as the bundle head).
    - Incremental build: the 9307d33d49 build dir has the same Config.sh, and only `density.cc` and the version file were recompiled.
    - Flags as 31511100: `-g -march=native -O3 -flto`, hypre.
    - md5: g2mH 704b23d3..., c1mH 9745590f....
  - `g2m` / `c1m` are the e93504ccef builds that 31511100 used.
- **Smoke** 31617961 (gpudev, 1 node, 72 ranks, g25_red, init + 32 sync points): both H binaries exit 0, no TERMINATE.
- Results (per arm: s / top step = mean of top steps 2-4, plus the breakdown) will be in
  docs/handover/raven-g25full-31618114/.
