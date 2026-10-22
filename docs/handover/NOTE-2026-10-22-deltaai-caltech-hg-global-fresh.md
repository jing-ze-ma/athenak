# NOTE 2026-10-22 deltaai -> caltech (USER): HG global FRESH (TASK c2c8bf6e) -- NOT RUN on DeltaAI (USER decision)

USER (10-22 ~10:40 CDT): do not run. No job submitted, 0 SU. Caltech's HG copy is the only one.

Estimates given to the USER (sbatch --test-only, 10-22 10:34 CDT, balance bivj-dtai-gh 833 SU):
- ghx4 3 nodes x 4 GH200, 48 h (12 SU/h, 576 SU/link): est start 10-27 16:44 CDT (~5 d) -- loses to Caltech tonight.
- ghx4 2 nodes, 48 h (8 SU/h): est start 10-27 16:30 CDT.
- ghx4-interactive 2 nodes, 2 h (x2 charge, 16 SU/link): est start ~35 min. Layout would be 8 GH200 x 3 MeshBlocks
  (704x64x64), ~68 GB/GPU expected (scaled from the hg145 wedge 21.8 GB/block-GPU and RSG E memcut 52.9 GB at 3 x 548).
Prepared (kept, unused): deck with @FILES@ substituted at /work/hdd/bivj/jma20/hgglob_1022 (input files md5 f9678c2a,
e2868b75, 23ee543f verified); incremental build of hesnf_gpu @ ff1f52bc started on the login node (no SU).
