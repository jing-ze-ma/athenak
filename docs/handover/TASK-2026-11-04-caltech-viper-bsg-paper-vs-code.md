# TASK 2026-11-04 Caltech -> viper (USER): Ma+2026 BSG — paper vs Linhao's code

The user asks: "what is the difference between Linhao's code and the description in his paper?"
Caltech cannot find the material to answer this properly. Please locate it and push it here
(branch bsg-files-1009, under `docs/handover/bsg-paper-code-1104/`), or reply with the paths if it is too large.

## Please send
1. **TRUEREPRO.md** (written for bsg-truerepro-1008, ~10-08): the side-by-side of the Ma, Bildsten & Jiang 2026 paper
   setup vs our reproduction, incl. the "DIFFERENT-can-match" list. Not on any branch Caltech sees
   (truerepro2-1009 / bsg-files-1009 / bsglinhao-1014 / fin4-linhao-bc have only decks and notes).
2. **Linhao Ma's Athena++ source** `bsg_star.cpp` (+ its athinput and any other problem files / opacity readers),
   i.e. what fin4-linhao-bc (eb1e032c7) and bsglinhao-1014 (d32aebaff) were ported from.
3. **The Ma+2026 paper PDF** (and arXiv id), if you have it.
4. Any existing paper-vs-code comparison (e.g. from the bsglinhao-1014 / fin4-linhao-bc work, 10-14).

## If you have time: the comparison itself
A table, paper text (section/equation) vs code (file:line), for: IC profile and its extension above r_max
(code: isothermal gas-only with g = g(rmax)(r/rmax)^2), inner boundary (paper text vs bd_data ghosts +
ghost-bath radiation), outer boundary, perturbation/seed (code: rho *= 1 + 0.05 (U - 0.5) where rho > 1e-11),
EOS (code P = rho T), gravity (point mass?), opacity tables (absorption/scattering split, padding/clamping), RT
(implicit DO; angles), grid/SMR, floors, run length / time units. Mark each: same / different / paper silent.

## What Caltech knows now (for reference)
- Caltech BSG runs: bsg_1009 truerepro2 (zero-contact = our reading of the paper), bsgR/RC 10-11, bsgprod_1012
  bsgP_seed_fin3 (39 d, 1e-6 T seed), bsgli_1014 (his IC + he_bc_linhao, 18.6 d).
- he_bc_linhao / he_ic_linhao are NOT in prod-1102 (only on fin4-linhao-bc / bsglinhao-1014).
- New plan (USER 11-04): BSG wedge = Ma+2026 settings except OUR inner wall, WB, and Y_lm seed (being prepared; not launched).

Reply with a NOTE on this branch.
