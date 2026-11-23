# NOTE 2026-11-23e orion -> Caltech: USER - build + smoke now, but WAIT for our GO NOTE before the production launch

USER (11-23 ~02:25 PST): the three fresh runs (bsg, m35mms, t9l62; original ICs; cfl 0.3) should use the FINAL keys of the FINAL
interim scheme. We are closing the last items now (one IDORT warning row and two uniformity arms), so:
1. Build caltech-ia9-77551e823 (CUDA 12.9 + g++13) and run your smoke as in 23c/23d. That build may be superseded by one more commit.
2. Do NOT launch production yet. Within the next few hours we post a GO NOTE with: the final branch/commit, and the exact <rad_m1>
   block to paste into all three decks (keys in the deck, not on the command line, per 23d).

Our own fresh smokes of the same three decks on IA9 (viper, 55 cycles each) passed: 0 FATAL / 0 CLEAN WARNING / 0 fallbacks /
0 dt_retry, dt equal to your production's (bsg 165.5, m35mms 12.0, t9l62 96.9), L_top/L_in close to your PROD1121b hst
(bsg 1.057 vs 1.052, m35mms 1.002 vs 1.000, t9l62 1.007 vs 0.994). Note: a fresh run writes an rst at t=0 even with
<output rst>/last_time=0 (bsg 22 GB, m35mms 29 GB, t9l62 88 GB) - plan disk accordingly.
Details: Orion /orion/ptmp/jinma/clean_1115/stag2/calsmoke (decks with keys in <rad_m1> in stage/*.athinput).

-- orion (coordinator)
