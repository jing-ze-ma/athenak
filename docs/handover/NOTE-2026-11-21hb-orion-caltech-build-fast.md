# NOTE 2026-11-21hb orion -> caltech: BUILD FAST — incremental builds from the B_v203n tree (USER: top priority now)

The fastest path is incremental make in the already-built CUDA 12.9 tree, not fresh configures:
- B_icbaln: 4a7f72334 differs from e688cb7bc in ONE file (src/pgen/he_star_m1.cpp). Copy the B_v203n binary aside, then in
  the B_v203n worktree: `git checkout 4a7f72334 && make -C build_n -j32` (same CUDA 12.9 module env) -> only he_star_m1.cpp
  recompiles + relink, ~1-3 min. Copy out as athena_icbaln_4a7f72334; then `git checkout e688cb7bc` back if needed.
- B_v204n: 7cbc13437 differs in 14 src files; build it the same way in a SECOND copy of the tree (cp -a the B_v203n worktree
  + build_n dir to a new dir, fix nothing else, `git checkout 7cbc13437 && make -j32`) on another node — only those TUs rebuild.
  If cmake complains about the copied build dir, re-run cmake in place with the cached flags (CMakeCache.txt is copied).
- Enable ccache with CCACHE_DIR=/resnick/groups/carnegie_poc/jingze/ccache_cu129 for these (never the 12.2 cache).
- Smoke-gate each binary (tiny BSG 1 GPU, pass 1 = 4.020e-06), then start the R_* arms of NOTE ha with B_icbaln immediately.
Report build seconds + md5.

-- Orion coordinator
