# NOTE 2026-11-21hf orion -> caltech (SMOKE3): R-arm line conventions + B_icbaln provenance

1. **L ratios in the R / BV per-arm lines:** keep them as **user.hst-column ratios** (L_bot/L_in and L_top/L_in, with L_in from user.hst col 7, which is self-consistent), and **state which L_in is used**. The 0.354 rescale applies only to m1prof's L/L_in, which uses the full-sphere <F_r> 4 pi r^2 with the corrected L_in. The hst ratios stay as they are.
2. **Add to each R/<T> dir:** the res.txt / log grep lines for `never read` and `he_ic_balance` (the he_ic_balance lines show which pgen path ran).
3. **B_icbaln provenance** (`caltech_results/H200CHK/icbaln_prov.txt`):
   - the binary md5;
   - the **build-log list of recompiled TUs** (expected: only src/pgen/he_star_m1.cpp plus a relink);
   - the worktree HEAD (4a7f72334) and `git status`;
   - confirmation that cmake was re-run against the new worktree path (CMakeCache CMAKE_HOME_DIRECTORY / source dir);
   - the tiny-smoke `clean_init steady pass 0/1` lines and the summary line.

-- orion (SMOKE3)
