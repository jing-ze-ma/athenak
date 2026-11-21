# NOTE 2026-11-21gn orion -> caltech (coordinator, URGENT and decisive): H200-only closure check + B_v203 build provenance -- AHEAD OF EVERYTHING

CLOSE finding: the f = 1/3 trapped-radiation atmosphere appears **only on the Caltech H200 runs**. The same decks at rev e688cb7bc on Orion CPU and on Raven A100 stream correctly, and CPU and A100 agree bitwise. All BSG / T19 Caltech arms ran on B_v203, so this check decides how to read every one of them. On **ONE free test node, ~1 min, before any other arm**:

## 1. Tiny BSG smoke on 1 GPU / 1 rank
- Binary **B_v203** = athena_v203_e688cb7bc (md5 7df2b648).
- Deck **bsg_wngl_smoke.athinput** (md5 84dbac02).
- CLI overrides (all keys exist in the deck): `mesh/nx2=32 mesh/nx3=32 meshblock/nx2=32 meshblock/nx3=32 time/nlim=2`.
- **NR=1** (one rank, one GPU). If a CLI key is rejected, copy the deck with those values instead.
- Report the log lines **`clean_init steady pass 0 ...`** and **`clean_init steady: ...`** verbatim, plus the whole log head up to cycle 1.
- Expected if healthy (Orion CPU / Raven A100):
  - **pass 0 = 9.534e-02**
  - **NOT converged after 50 passes**
  - **E -> 4.1604261309e+45**
  - **F/c -> 1.3591e+42**
- **Also reported, as a reference only:** in the Caltech smokes so far, clean_init converged in 3 passes.

## 2. B_v203 build provenance on H200
Please post:
- the .prov file, or the git rev + source digest of the source tree that built athena_v203_e688cb7bc, and `git status` / `git diff --stat` of that tree;
- **CMakeCache.txt** entries: Kokkos arch (e.g. Kokkos_ARCH_HOPPER90), Kokkos_ENABLE_CUDA, CMAKE_BUILD_TYPE, CMAKE_CXX_FLAGS / CMAKE_CUDA_FLAGS, **ATHENA_CUDA_BODY_GLOBAL**, **VGD_CS_LB**, **Athena_EXTRA_DEFS** and any other -D defines, plus `-arch` / sm_90;
- compiler + CUDA versions (`nvcc --version`, host compiler, the module list used);
- the build dir path, and whether that build dir was ever copied from another slot.

Commit the results to `caltech_results/H200CHK/` (log head + cmake grep + versions), one commit.

-- orion (coordinator)
