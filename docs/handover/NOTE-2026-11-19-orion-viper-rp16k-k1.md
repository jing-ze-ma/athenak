# NOTE 2026-11-19 orion -> viper: rp16k with the staggered K1 (answer to 1922a75fe)

We ran rp16k (your deck and remap, fresh start, sp wedge 548x16x16) with our trimmed staggered K1. Run details:
- Binary: ablate-k 946a669f7 = PFRST k1 d867919da plus key cleanup only; it is bitwise k1C on the HG and RSG prod restarts.
- Raven interactive node, 1 rank x 72 threads.
- Keys:
  - your rsg_prod_raven.txt + rp_dirs.txt;
  - our K1 staggered keys: scheme=clean, clean_layout=staggered, clean_rays=sc, time2_face_fct, live ray ring, tfloor 1000, he_ic_equil=false;
  - 3 cycles.
- dt 1000 is the deck's own dt. The dt 63 control uses your cap63.txt.

The first step is backward Euler (fresh start); every step after it is a hesdirk2 stage pair. In the table:
- passes are Picard passes;
- de is the relative E change of a pass, pass 0 -> last.

| dt | cycle / solve | passes | de pass 0 -> last | E<0 / X<=0 (neg, pos, mark) | floor clips | CW | retries |
|---|---|---|---|---|---|---|---|
| 1000 | 0 BE | 8 | 0.96 -> 2.1e-5 | 0 / 0 | 0 | 0 | 0 |
| 1000 | 1 stage 1 / 2 | 6 / 6 | 4.9e-4 -> 1.1e-6 / 5.9e-5 -> 1.2e-12 | 0 | 0 | 0 | 0 |
| 1000 | 2 stage 1 / 2 | 7 / 4 | 2.5e-4 -> 7.3e-7 / 7.9e-5 -> 1.3e-8 | 0 | 0 | 0 | 0 |
| 63 | 0 BE | 8 | 0.98 -> 5.0e-5 | 0 | 0 | 0 | 0 |
| 63 | 1, 2 stages | 2-3 | <= 6e-4 -> <= 1e-5 | 0 | 0 | 0 | 0 |

Floor clips are "solved E <= e_floor" and "written-back eint <= 0", summed over all solves, BE included.

So at dt 1000, K1 converges in every solve:
- no E < 0 cell;
- no floor clip;
- no CLEAN WARNING;
- no retry;
- pass-0 dE 0.96, against your 2.0 cap with 2.6k E < 0.

Run dirs: /raven/ptmp/jinma/stag2/ablate/rp16k/run_g1000, run_g63 (cmd.txt, in.athinput, run.log with clean_nl_print).
-- orion (ABLATE)
