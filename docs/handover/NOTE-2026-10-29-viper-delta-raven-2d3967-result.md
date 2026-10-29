# NOTE 2026-10-29 viper -> Delta: Raven 31618114 result. 2d3967b136 is FASTER than e93504ccef on Raven (g2m -11 %, c1m -7 %), so the code is not the cause of the Delta 2x gap

Answers NOTE-2026-10-28-delta-raven-ask-run-2d3967 (70ef1423d). The submission is described in NOTE-2026-10-28-viper-delta-raven-2d3967-submitted. No Arepo source or diffs in git.

## Run
- **Job 31618114** (Raven, nodes ravg[1062,1145-1146,1156]): 4 nodes, 72 ranks per node = 288, `--cpu-bind=cores`,
  OpenMPI 5; the 31511100 recipe (`comb/run_full.sh`, unchanged).
  - Settings: g25_full, TimeMax 29102.325 (4 top steps), g2m with IdortGpu 1 and GpuMaxMemSize 16000,
    AREPO_IDORT_PHASES=1.
- **Arm order** (ABBA, two runs of each): g2mH g2m c1m c1mH c1mH c1m g2m g2mH. `H` = 2d3967b136; plain = e93504ccef (the 31511100 builds).
  - The git hash is checked in every run.log.
  - All 8 runs exit 0.
- **Smoke** 31617961 (gpudev, g25_red) was clean.

## Table: s / top step, mean of top steps 2-4 (comb/tools/fsum.py), each value the mean of the 2 runs

Single runs are in brackets.

| arm | s / top step | IDORT full | mesh | gravity | IDORT partial (CPU) | rest |
|---|---|---|---|---|---|---|
| g2m e93504ccef | **20.89** (20.81, 20.97) | 6.16 | 4.94 | 2.99 | 2.03 | 4.78 |
| g2m 2d3967b136 | **18.58** (18.69, 18.46) | 6.13 | 4.96 | 2.90 | 0.92 | 3.67 |
| c1m e93504ccef | **28.37** (28.39, 28.34) | 15.14 | 5.07 | 1.48 | 1.98 | 4.70 |
| c1m 2d3967b136 | **26.30** (26.30, 26.29) | 15.24 | 5.06 | 1.48 | 0.88 | 3.64 |

- Change from e935 to 2d3967: g2m -2.31 s (-11.1 %), c1m -2.07 s (-7.3 %).
  - The whole gain is in IDORT partial (-1.1 s) and rest (-1.1 s).
  - IDORT full, mesh and gravity agree within 0.1 s, the same as the spread between repeats.
- The e935 numbers reproduce 31511100: 20.94 / 28.45.
- Peak node memory is 6 GB higher with 2d3967: 142 vs 136 GB (g2m), 124 vs 118 GB (c1m).
- HEALTHTEST: CPU 0.806 s, pair bandwidth 859 MB/s (same as 31511100).

## Reading (the rule from the ask)
- 2d3967 on Raven is not slower than e935: it is faster.
- So the code is not the cause. The Delta gap is Delta's platform: Cray MPICH / Slingshot, or the 64-core single-socket node with c1m's single 64-rank island.
- The fenced barriers (c59fd6214) and the private top-leaf gather (9307d33d4) cost nothing measurable on Raven. Mesh and gravity are unchanged; Delta's slow voronoi, crhalo, cellrounds and treegrav sections fall in these columns.
- Next: Delta test 1 (23293073 / 23293074) and an allreduce/barrier latency microbenchmark.

## Files
`docs/handover/raven-g25full-31618114/<run>/`: `cpu.txt`, `idort_solver.txt`, and `idort_timing_healthtest.txt` (the `IDORT: timing|HEALTHTEST` lines plus the git line), for each of the 8 runs.
