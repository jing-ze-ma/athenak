# NOTE 2026-10-19 caltech -> raven: H200 determinism check submitted (reply to TASK 29dba835)

- Build a2713be0 H200 nofma md5 e31ea39ea5e05148e2dbea0b95328b27.
- 5001051: production deck + production layout (2 nodes x 4 H200, 8 ranks, 3 x 591x64x64 per GPU): C 200 cyc, C2 rerun, R from C cycle-100 rst. Slurm estimate start ~10-20 19:26, i.e. after link 1 (est 10-19 16:47).
- 5001065: reduced smoke layout (cs N=48, 24 x 591x24x24, 1 node x 4 H200 = 6 blocks/GPU), same C/C2/R, should run earlier.
- Production links 4993305-07 unchanged (READY=1). On any DIFF Caltech holds them (READY=0) and reports.
