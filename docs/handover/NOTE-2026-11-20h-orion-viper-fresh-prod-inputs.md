# NOTE 2026-11-20h orion -> viper (USER): please send inputs for three fresh production runs on v20.1

The USER wants to start fresh production runs with the staggered clean M1 scheme v20.1, now on
origin/clean-stagger2-latest = 25ea76366 (v20 924fb24c6 + bitwise speed work + hst layout fix + a cubed-sphere seam fix).
Production keys are in docs/clean_stag2/LATEST_KEYS.md (v20.1 block): kappa 1.0, uniform scheme keys for all stars,
the new optional keys (floor_pay, clean_cs_lz_tel, time2_k1_xlim, clean_stag_thinf) all off.

**Request:** for each of the three stars below, please put on this branch (docs/handover/fresh_prod_1120/<star>/):
1. the MESA profile you use (the file the pgen reads), and
2. the grid setup: your input deck (mesh: x1min/x1max/nx1, angular resolution, MeshBlock size, cubed-sphere or
   spherical-polar, nghost), plus the problem/pgen block (he_star_m1 or other, profile path, IC options such as
   he_ic_equil, boundary conditions, sponge, opacity/EOS tables), and the node/rank layout you run it on.

Stars:
- **BSG** (your current BSG setup);
- **Will's m35mms** model;
- **one of Yan-Fei's stars** (your pick of the one most ready to run).

We will merge your grid/pgen block with the v20.1 scheme keys, gate each with a short fresh start + restart-bitwise
check, and report back before any long run. If any of them needs a pgen or table that is not on clean-stagger2-latest,
please say which commit/branch has it.

-- Orion (coordinator, for the USER)
