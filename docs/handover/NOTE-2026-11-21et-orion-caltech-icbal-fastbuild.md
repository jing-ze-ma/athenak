# NOTE 2026-11-21et orion -> caltech (USER, URGENT): build B_icbal INCREMENTALLY, then start MB_bsg62m + F_* at once

The USER wants the icbal build fast. Between e688cb7bc (B_v203) and 4a7f72334 (icbal) **only src/pgen/he_star_m1.cpp differs** (git diff --stat: 1 file). Please rebuild incrementally instead of from scratch:
1. If a full icbal build is already more than about half done, let it finish. Otherwise stop it.
2. Incremental build. Use the source tree whose build dir produced athena_v203_e688cb7bc, so the CMake cache points at the right sources.
   - Keep the existing B_v203 binary where it is; do not overwrite athena_v203_e688cb7bc.
   - `git checkout 4a7f72334` in that same source tree.
   - `make -j<all node cores>` in its existing build dir: CUDA arch sm_90 only, PROBLEM=he_star_m1 only (as configured for B_v203).
   - Expect exactly **1 object (he_star_m1.cpp) recompiled + relink**.
   - Copy the binary to bin/athena_icbal_4a7f72334.
   - If the B_v203 build dir is no longer usable, a copied build dir is OK **only** if you re-run cmake against the 4a7f72334 source path. A copied CMake cache pointing at another tree once produced a stale binary for us, which is why we ask for these checks.
3. **Provenance checks** (all required):
   - .prov / version string shows rev 4a7f72334.
   - The build log shows he_star_m1.cpp compiled; 0 TUs compiled = stale.
   - The md5 differs from athena_v203_e688cb7bc.
   - The run log of the first arm shows the new keys being read (he_ic_balance_aref / the discrete-force re-march message).
4. Start **MB_bsg62m + the 5 F_* smokes (NOTE eq) immediately** after the build, ahead of everything (NOTE es).
Please post the **build start and finish times** (and the number of objects compiled) in your next line.

-- orion (for the USER)
