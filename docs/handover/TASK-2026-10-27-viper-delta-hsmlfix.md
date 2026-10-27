# TASK viper -> Delta (USER): test the Hsml race fix hsmlfix-1026 on g25_full at your failing layout, one short job

**Wait for the bundle.** The user copies it into your Arepo bundle directory (next to the d9c0152bce one):
- `arepo_hsmlfix_1026.bundle`, md5 **f547eaaf285bd307c72a32c5483801e2**, 6652 bytes.
- Head: branch `hsmlfix-1026` = **2d3967b136**. Prerequisite: **d9c0152bce** (you already built it for job 23219344).
- Run `git bundle verify` first. If the bundle is not there, or verify fails, stop and report. Do not run anything else.

**What the fix is (3 commits on d9c0152bce):**
- **9307d33d49, the suspected fix:** a race on a shared buffer in the top-leaf exchange (`ngb_exchange_topleafdata`).
  - Every rank of an island did `MPI_Allgatherv` straight into one node-shared (`mymallocshared`) buffer.
  - Now each rank gathers into private memory, and the island master alone publishes it (fenced), via the new
    `myMPI_Allgatherv_to_shared`.
  - Sites: ngbtree.cc:846, forcetree.cc:2257, sidm.
  - This bug class matches the DeltaAI island-sort hang (fixed there by `MPICH_SHARED_MEM_COLL_OPT=0`), and your V4 tests it independently.
- **83c1a8573c + 2d3967b136:** the env-gated init diagnostic `AREPO_HSML_DIAG`, off by default.
  - `=1` prints every cell whose converged Hsml is > 0.2 box: initial guess, NumNgb, Left/Right.
  - `=2` also prints every boundary cell.
  - Lines start with `HSMLDIAG`.
- **Raven gates on the fix, all PASS:**
  - g25_full CPU init at 256/192/128 ranks: the ngb iterations, VCELL build checksums and Sync-Points are identical to d9c0. Hsml is 0.062 box, 17 iterations.
  - IDORT I4/A2/B3/I1/A4/A7, CPU vs GPU: solver logs ALL SAME.
  - GPU smoke: PASS.
  - Raven never reproduced the abort, so only your machine can show whether the fix cures it.

**Run (one job, the job 23219344 layout, about 20 min):**
1. Build **g2m** and **c1m** from 2d3967b136, using exactly the recipe and modules of job 23219344 (gcc 13 + CUDA 12.9, sm_80).
2. Use 4 gpuA100x4 nodes, 256 ranks, 64 per node, 16 per island. Do **not** set `MPICH_SHARED_MEM_COLL_OPT`.
   Set `AREPO_HSML_DIAG=2`.
3. g2m first, then c1m. Each stops at the first Sync-Point (PASS) or TERMINATE (FAIL).

**Reply NOTE:**
- PASS/FAIL per arm.
- The number of `HSMLDIAG` cells with Hsml > 0.2 box, and the converged Hsml of the 14 IDs from d6aa24a9.
- The density-loop iteration count.
- Your V1-V4 results, if they have run (V4 = d9c0 with `MPICH_SHARED_MEM_COLL_OPT=0`).
- SU used.

Then stop. Leave out the VCELL seed-cap workaround. No Arepo source or diffs in git.
