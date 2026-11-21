# NOTE 2026-11-21go orion -> caltech (coordinator): PRE-AUTHORISED follow-ups to NOTE gn (H200CHK) -- do them without waiting for us. Test nodes only.

CLOSE confirms that the **full production layout** (256^2 lateral, 128^2 blocks, 4 ranks) **streams** on Orion CPU: pass 0 is 9.63e-2, with no 4-pass convergence. So the trap comes from either the **H200 binary** or the **Caltech multi-GPU environment**. Depending on the tiny 1-GPU smoke of NOTE gn:

## (A) If the tiny 1-GPU smoke TRAPS (pass 0 ~1e-6, converges in ~4 passes) -> suspect the build
Rebuild **e688cb7bc** in a **fresh, separate build dir** (not B_v203's, and not copied from any slot) with the Raven integ3x_cu recipe, using Hopper instead of Ampere:
```
cmake -G Ninja -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=<kokkos>/bin/nvcc_wrapper \
  -D CMAKE_CXX_FLAGS=-ffp-contract=off -D Athena_ENABLE_MPI=ON \
  -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_HOPPER90=On -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF \
  -D Athena_M1_VGD_NMAX=512 -D PROBLEM=he_star_m1
```
- Use **no Athena_EXTRA_DEFS, no ATHENA_CUDA_BODY_GLOBAL, and no NCCL**. Use the Kokkos the source tree pins (Raven used Kokkos 4.6.00).
- At runtime, set `Kokkos_CUDA_MALLOC_ASYNC=OFF KOKKOS_CUDA_MALLOC_ASYNC=OFF OMP_NUM_THREADS=1`.
- **Note: Raven builds with `-ffp-contract=off`.** Please tell us whether B_v203 was built with it.
- Name the result **B_v203r**, record its md5 and provenance, and rerun the tiny 1-GPU smoke with it. Report the two clean_init lines.

## (B) If the tiny 1-GPU smoke STREAMS (pass 0 ~9.5e-2, not converged after 50) -> suspect the multi-GPU environment
With the same tiny deck on **4 GPUs / 4 ranks** (CLI `mesh/nx2=32 mesh/nx3=32 meshblock/nx2=16 meshblock/nx3=16 time/nlim=2`, i.e. 2x2 lateral blocks of 16^2):
1. B_v203, as normal.
2. B_v203 + `rad_m1/vet_gd_halo_ipc=false`. If that key is not in the deck, copy the deck with it added under <rad_m1>.
3. Like 2, plus UCX GPU-direct off: `UCX_IB_GPU_DIRECT_RDMA=0` and `UCX_TLS=^cuda_ipc`, or the closest equivalent in your MPI environment.
Report both clean_init lines for each of the three runs, plus the env (UCX_* / OMPI_* variables) used for each.

## Optional, if easy: invalid-closure fraction
`docs/handover/caltech_runs_1121/close/inval.py <m1_vet bin> [stride] [path to vis/python]` prints the fraction of cells per shell whose m1_vet off-diagonals d12, d13 and d23 are exactly 0, i.e. the isotropic path. Please run it on one m1_vet bin of a full-res H200 BSG smoke or production run, plus any m1_vet bin from the tiny runs.

Commit everything to `caltech_results/H200CHK/`.

-- orion (coordinator)
