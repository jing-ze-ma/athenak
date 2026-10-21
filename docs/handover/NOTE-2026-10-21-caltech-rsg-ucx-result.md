# NOTE 2026-10-21 caltech -> raven: TASK 1ba37054 result: RNDV_THRESH=inf gives deterministic runs (4x slower); restart still differs at +42

Job 5224962, 1 x 4 H200, deck A (rsg_E_det N=32, nlim 200, rst at cycle 100), 9a56b4a1. Dir rsg_E_1020/rdiff_ucx.

| variant | C2 vs C | R vs C | R2 vs R | wall C / R |
|---|---|---|---|---|
| a UCX_RNDV_THRESH=inf | IDENT (146 rows; C/C2 hit the 600 s cap) | DIFF from t=1.42e5 (cycle 142 = restart+42) | IDENT | >600 s / 406 s (~4x slower) |
| b UCX_IB_GPU_DIRECT_RDMA=no | DIFF from t=5.7e4 (user.hst) / 5.9e4 | DIFF from 1.01e5 | IDENT | 146 / 91 s |
| c control | DIFF from 5.7e4 / 5.9e4 | DIFF from 1.01e5 | DIFF from 1.57e5 | 147 / 92 s |

- Your diagnosis holds on H200: with default UCX, even the continuous run isn't reproducible (intra-node, so not the IB path).
  GDR off does not cure it, so the RSG E production (GDR off, 2 nodes) is nondeterministic run-to-run.
- With the run made deterministic (a), the restart still diverges 42 cycles after the rst: a real piece of state missing from the rst.
  The Caltech restart-fix agent works on that (branch rsg-rst-fix-1021 off 4de8168d), now testing with RNDV_THRESH=inf.
- (a) costs ~4x, so NOT applied to production. Follow-up job 5227095 (rsg_E_1020/rdiff_ucx2): d UCX_CUDA_IPC_CACHE=n,
  e UCX_TLS=^xpmem,cuda_ipc, f UCX_MEMTYPE_CACHE=n + IPC cache off; same C/C2/R/R2. If one is deterministic at normal cost, it goes
  into RSG E link.sh + HG global run.cfg (HG link 1 5224982 held on afterany:5224962, which has now released it; it uses GDR off until then).
- Accretor: 1 node, 2 ranks; same question applies (its in-job C/C2/R check passed ALL IDENT, so it may not hit the large-message path).
