# NOTE viper -> caltech (USER): FINAL reply to TASK 8513c5b56 (g25_full on Raven): 4-node results only; the 8-node job did not run

Follows the interim NOTE 01508663c. No Arepo source in git.

## What ran
- **8-node job 31511099** (576 ranks, g2m/c0m/c1m): **did not run**. sacct shows it CANCELLED before it started (submitted 10-25 18:23, cancelled 10-26 16:45). The reason is not in viper's notes. So there are **no 8-node / 576-rank numbers.**
- **Two 4-node runs** (288 ranks, 72/node, 4 A100/node, `srun --cpu-bind=cores`, OMP 1, OpenMPI 5): job script `/raven/ptmp/jinma/arepo2_gpu/comb/run_full.sh`, logs `comb/logs/full_<job>.out`, runs `comb/runs/full_<arm>_4n_*_<job>`.
- s / top step: top step = 7275.58 s sim = 32 steps of 227.36 s, with 1 full IDORT solve. Mean over top steps 2-4 (min-max), from `comb/tools/fsum.py`.

| job | code (git) | arm | s / top step | IDORT full | mesh | gravity | IDORT partial (CPU) | rest | imbalance max/mean | peak node mem GB |
|---|---|---|---|---|---|---|---|---|---|---|
| 31511100 (10-26) | e93504ccef | **g2m GPU recipe** (IdortGpu 1) | **20.94** (19.31-22.45) | 6.15 | 5.12 | 2.91 | 2.01 | 4.74 | 9.3 | 135 |
| 31511100 | e93504ccef | c1m CPU node-local | 28.45 (26.59-29.77) | 15.23 | 5.08 | 1.47 | 1.98 | 4.69 | 6.1 | 117 |
| 31546612 (10-28) | 3405e3ab39 | g2 GPU recipe | 27.44 (25.41-30.16) | 9.09 | 6.98 | 3.42 | 2.60 | 5.35 | 9.3 | 148 |
| 31546612 | 3405e3ab39 | c0 CPU classic | 34.55 (32.75-37.91) | 16.58 | 7.00 | 2.98 | 2.84 | 5.14 | 6.3 | 128 |
| 31546612 | 3405e3ab39 | c1 CPU node-local | 32.82 (31.39-34.83) | 16.48 | 6.87 | 1.56 | 2.75 | 5.17 | 6.1 | 123 |

- g2 and g2m have the same Config.sh options; they differ only in the code commit. 31546612 is the older code (3405e3ab39, the comb branch with the g25 direct-gravity fix). 31511100 is the newer e93504ccef.
- The newer code is faster in every section, mesh included (5.1 vs 7.0 s). Run-to-run and node variation is not separated from the code change: different days, different nodes, one run each.
- Peak GPU memory is 117 GB per node (of 160) in both GPU arms.
- **Against Orion** (35.6 s on 560 cores): the 4-node GPU recipe is 1.7x faster (20.94 s). The 4-node CPU arms are 1.0-1.25x faster on 288 ranks.

## Delta, for completeness (see the Delta notes; same g25_full)
- Delta 23248921 (4 A100 nodes, 256 ranks, hsmlfix 2d3967b136): g2m 55.8 / c1m 62.9 / c0 67.2 s per *Delta* top step. But its g2m ran IDORT on the CPU (IdortGpu unset), and its top step held 0.8 full solves.
- Delta's CPU IDORT is 2.8x Raven's and MPI-bound.
- A like-for-like rerun (IdortGpu 1, Raven's TimeMax) is queued on Delta: 23269011/12.

## Not done (from the TASK)
- 8 nodes / 576 ranks (the cancelled job).
- The extra 4-node scaling arm exists only as the 288-rank runs above.
- Whether to submit an 8-node job now is the USER's call; nothing is queued.
