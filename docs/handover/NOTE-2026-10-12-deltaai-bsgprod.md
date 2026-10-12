# NOTE DeltaAI: REAL BSG true reproduction, final method, chained 1 h links (TASK-2026-10-12-deltaai-bsgprod) -- STARTED (job 3367085)

**Picked up** 11:16 CDT (TASK 40558c5d, user GO). bsgcand-1012 **d3bf42d7** built: `athena_hes_gpu_d3bf42d7dc27`,
md5 63528bf7f170366e4fd5cbc05f4b74b6 (deck keys implicit_thin_js_face, vet_gd_rotate_slope, vet_gd_halo_ipc, vet_gd_twin_lowmem present).
- Files: the 10-09 BSG files (md5 ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805). Deck bsgP_seed.athinput.in
  (md5 e193df54, after the 53995cf9 fix) with @FILES@ resolved only.
- Layout 1 node x 4 GH200, 4 ranks x 4 MeshBlocks, plain srun, MPICH_GPU_SUPPORT_ENABLED=1.
- Link 1 = job **3367085** (ghx4-interactive, 1 h): 25-cycle smoke (nlim 25, ndiag 1, no rst; rc / FATAL / NaN / NON-CONVERGED -> STOP),
  then FRESH production with `-t` = wall left - 5 min. ic0check run on the smoke dumps from the login node while production runs.
- Later links: newest rst + rst_info.py outputN/last_time keys; STOP on rc != 0 / FATAL / NaN / DIVERGED / NON-CONVERGED > 0; DONE at t >= 4.96e6.
  A login-node chain submits each link when the previous one ends (1 job per user on ghx4-interactive); it stops if a bsgprod
  CANCEL / WITHDRAWN NOTE appears. Cap 12 links.
- Dir /work/nvme/bivj/jma20/bsgprod_1012 (run/, smoke/, logs/link.<job>.out). STARTED line follows as soon as link 1 runs.

**STARTED**: link 1 job **3367085** RUNNING on gh016 since 11:23:07 CDT (smoke 25 cycles, then fresh production in the same job). Smoke numbers and ic0check follow.
