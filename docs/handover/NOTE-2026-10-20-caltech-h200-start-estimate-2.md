# NOTE 2026-10-20 caltech -> raven: fresh H200 start estimates (reply to a44616c3) -> Caltech is NOT competitive

`sbatch --test-only` 10-19 21:43 PDT (10-20 06:43 CEST), partition gpu, --exclude=hpc-sm-01-09, 120 GB host mem/GPU. 12 h and 24 h give identical estimates.
- 8 H200 (2 nodes x 4): 10-23 22:14 PDT  (~96 h from now)
- 6 H200 (2 nodes x 3): 10-23 14:07 PDT  (~88 h)
- 4 H200 (1 node x 4):  10-22 19:12 PDT  (~69 h)
None is earlier than ~10 h; Raven (~15 h) / Viper are faster. Reason: partition crowded (857 pending) and our fair-share factor is low (~0.22).
Our 2-H200 acc3d check 5087543 is itself estimated ~15 h out. No build TASK needed for Caltech.
