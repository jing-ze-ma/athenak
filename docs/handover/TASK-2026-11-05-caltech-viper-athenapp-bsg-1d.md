# TASK 2026-11-05 Caltech -> viper (USER): did we run Linhao's Athena++ BSG (1-D / column) before?

The user remembers that we ran **Athena++ itself** (Linhao's bsg_star.cpp, implicit DO radiation) in 1-D or as a thin column
before, to understand the BSG behaviour. Caltech finds no record (no handover note, memory entry or commit mentions an Athena++
run; only AthenaK column checks: viper's 10-06 4x4 column in MA2026_REPRO_IC.md and the bsg_m1_column step-3 tests).

Please check viper/Raven (and any DeltaAI/orion notes) and reply with a NOTE on this branch:
1. Was Athena++ (stock or Linhao's fork) ever built and run for the BSG? Where (paths, dates), which setup (athinput,
   1-D/column/3-D, angles, grid, IC), and what it showed (near-base Mdot, L_in/L_bot, R_ph(t), shrinking).
2. If yes: push the athinput(s), the build recipe, and the key results/plots under docs/handover/athenapp-bsg/ (paths if large).
3. If no: say so. Caltech will then build stock Athena++ (PrincetonUniversity/athena @053bc738, nr_radiation implicit) with
   bsg_star.cpp (minus the non-public rad_source_term output) and run 1-D columns of his exact setup (r^+2 atmosphere and
   invsq variants) next to matching AthenaK 1-D arms.

Context: our Linhao-exact 3-D runs show the paper's early near-base flow (3.4/yr at 0.5 d) but then only 0.3-0.6/yr sloshing
(paper: "often exceeding 4/yr" for the first ~10 d) and no shrinking (r_photo,3D 51-60 Rsun vs the paper's ~43). The 10-14
bsgli run (his r^+2 gravity + old M1 face) had a sustained 1.7-2.6/yr at the base. Caltech is running AthenaK 1-D arms now.
