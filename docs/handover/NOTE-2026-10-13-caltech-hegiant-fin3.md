# NOTE Caltech -> Raven coordinator: He giant FIN3 production (TASK-2026-10-13-caltech-hegiant-fin3)

## SMOKE PASS + CHAIN QUEUED (10-13 11:45 PDT)
- Binary `athena_gpu_he_star_m1_dad00933_nofma` (dad00933acc5..., incremental from the 9972a56d tree), **md5 ff7b7f11170340bee7d53dacec851ae9**
  (same binary as the Caltech E1 copy).
- Files: hegiant_1009/files (ic_giant_own ddc72692, rosseland c1ef063c, planck b699c60d). Deck hegiant_fin3.athinput.in (md5 ba10350a) -> @FILES@ sed only
  -> `hegiant_1013/hegiant_fin3.athinput` (read-only) **md5 02f8fd7ac8fc063f67dfa35b9552767b**. vet_gd_twin_lowmem NOT set (fits).
- Layout: 1 node x 2 H200, 2 ranks x 2 MeshBlocks (897x64x64), CUDA_VISIBLE_DEVICES=$SLURM_LOCALID wrapper, --cpu-bind=cores.
- **Smoke job 4452306** (hpc-sm-01-05), keys as TASK sect. 4:
  - A fresh 24 cycles: rc 0, FATAL 0, NaN 0, DIVERGED 0, NON-CONVERGED 0; dt 20.0 s on every cycle (25 lines);
    Picard mean 5.33 max 10 (88 solves); "vet_rotate_avg = 4" line present; he_ic_equil = false so no equil line;
    he_ic_balance max |rho/rho_col-1| 5.1e-3, face residual 9.5e-15.
  - **s/cycle (cycles 3-24) = 2.74** (elapsed 63.5 -> 121.2 s). Init ~55 s. GPU peak **85.7 GB** per H200.
  - B (-r A/rst/hegiant.00001.rst, cycle 12 -> 24): rc 0, clean, Picard 5.00/6.
  - **Restart BITWISE**: final hydro_w / m1 / m1_vet dumps SAME (data after the header); user.hst 12/12 and hydro.hst 12/12 B rows verbatim in A.
- **Production chain** run dir `hegiant_1013/fin3` (fresh t = 0, tlim 2.592e6): links **4454171 -> 4454172 -> 4454173 -> 4454174 -> 4454175 -> 4454176**
  (12 h, afterany). Link script = BSG link_bsg.sh pattern (+ DIVERGED stop), prune keep 3 (prune_rst.sh), DONE at rst t >= tlim. Start estimate: .
- **Link count**: at 2.74 s/cycle and dt 20 s, 1.3e5 cycles is ~99 h = **~9 links**, not 5-6 (2 H200 x 2 blocks is about 2x slower than 4 A100 x 1 block).
  6 are queued; more will be added after link 1 measures the real rate (user decides).
