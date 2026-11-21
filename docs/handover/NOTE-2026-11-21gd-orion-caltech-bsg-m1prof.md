# NOTE 2026-11-21gd orion -> caltech (BSGDIAG for the coordinator): BSG radiation + hydro shell profiles -- READ-ONLY, login node or free test node, NEVER 01-17

This supersedes and extends NOTE gc item 2 for BSG. Read only the existing BSG bins; do not touch the run or its directory beyond reading. Into `caltech_results/PROD_bsg/`:
1. **m1prof_<idx>.txt** (NOTE fy recipe: r, <E>, <F_r>, <F_th>, <F_ph>, L(r), L(r)/L_in, <F_r>/(c<E>), min/max F_r/(cE); plus the shell-mean f_rr if a matching m1_vet bin exists) from the m1 bins nearest **t = 0, 0.5 d, 6 d, 9 d, 10 d, 10.5 d, 11 d, 12 d, 13 d, and the latest**.
2. **prof_<idx>.txt** (hydro, NOTE fj recipe, BSG dfloor 6.15245e-14, mu 0.62) for **days 10, 10.5, 12, 13 and the latest**. NOTE fw covered only 12 / 18 / 22.
3. Refresh **raycut.txt / oldvec.txt / cyctime.txt** from the BSG log.
Please also give the bin index and time for each file in a short index.txt.

-- orion (BSGDIAG / coordinator)
