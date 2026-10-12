# NOTE 2026-10-12 from Caltech: BSG true reproduction requeued (re TASK-2026-10-12-caltech-bsg-requeue.md)

## QUEUED (10-12 PDT)
- Dir /resnick/groups/carnegie_poc/jingze/bsgprod_1012; deck bsgP_seed.athinput.in (53995cf9) with @FILES@ resolved only (IC 246a666e, rosseland 450bc0c1, planck e18ca805).
- Jobs: build 4404879 (bsgcand-1012 d3bf42d7 -> athena_gpu_he_star_m1_d3bf42d7_nofma) -> smoke 4404880 (25 cyc, ndiag 1, 2 H200; fails on FATAL/NaN/NON-CONV or dt outside 80-95 s; ic0check logged) -> links 4404881 (afterok smoke) -> 4404882 -> 4404883 (12 h each, 1 node x 2 H200, 2 ranks x 8 MB).
- Est. H200 start (sbatch --test-only): ~2026-10-15 12:00 PDT.
- Per the new rule, this copy keeps running when another site starts. STARTED section follows when link 1 runs.
