# NOTE 2026-10-21 caltech -> raven (+ all sites): H200 run-to-run reproducibility fix = UCX_RNDV_SCHEME=put_ppln (with GDR off)

Screens on RSG E (deck E). C vs C2 = two identical runs from t=0; "IDENT" = every hst row bitwise.

1 node x 4 H200 (deck A N=32, 200 cycles, 9a56b4a1):
- default, GDR off, CUDA_LAUNCH_BLOCKING=1, pml ob1/smcuda, ^gdr_copy: DIFF from cycle ~57 (so not a missing sync in our code).
- UCX_RNDV_THRESH=inf: IDENT but ~4x slower. UCX_CUDA_IPC_CACHE=n: >4x slower. UCX_TLS=^xpmem,cuda_ipc: crashes.
- UCX_RNDV_SCHEME=put_zcopy: IDENT, +5 % (rd3_g 5227350; restart triplet 5228230: C2==C, R2==R, R!=C from cycle 101).

2 nodes x 4 H200 (prod deck + layout, 63a7d28b, GDR off, 100-120 cycles):
- GDR off alone (production until now): C2 != C from cycle 57.
- put_zcopy: one run MPI_ERR_INTERN after 37 s, the other OK but 2.7x slower -> unusable across nodes.
- get_ppln: DIFF from cycle 57.
- **put_ppln: IDENT over 101 cycles at full speed (2.98e7 vs 3.03e7 zone-cycles/s for get_ppln). Job 5229321, dir rsg_E_1020/ucx2n/pp.**
- am, RNDV_THRESH=inf: still running.

Applied: UCX_ENV="UCX_IB_GPU_DIRECT_RDMA=no UCX_RNDV_SCHEME=put_ppln" in rsg_E_1020/run63.cfg (RSG E, from link 2 on;
link 1 5227533 runs without it) and hgglob_1021/run.cfg (HG global link 1 5224982, released now).
Remaining: the real restart bug (R != C from the first cycle after a restart, with a deterministic env) is being worked on (branch rsg-rst-fix-1021).
Suggest testing UCX_RNDV_SCHEME=put_ppln (or the site equivalent) on Raven / DeltaAI / Viper multi-node runs too.
