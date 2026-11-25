# NOTE 2026-11-24zzzzc3 caltech -> orion: production c4 (t9l62) on caltech-c4-a000d803c (your 24zzzzc3; USER go relayed)

- Binary athena_c4_a000d803c md5 3ae130e97b13e5bdfa4a7f9a14370b9f (checked by the launcher before srun; contains the ray-hole cure); deck candidate_1124 IA123_t9l62 md5 f2e3b57d239f67e20b88ab89987e3323, copied unchanged (re-read from the branch: same md5). **Deck diffs: none.** No command-line keys. thinf = true, as in the deck. Deck outputs, rst cadence and final rst as in the deck.
- Node 01-15, holder 7333036 (ends 11-26 23:27:36; END 23:07; 53 h margin at start, the launcher refuses a fresh start below 34 h), 4 H200, 4 ranks, NEW dir prod_c4_1124/t9l62. The node was clean (c4_m35mms and its analysis step had ended). bsg, m35mms, M35ZAMS on 95813cf85 untouched.
- Start line (chunks.txt):
```
2026-11-24 17:52:04 chunk 0 STAR=t9l62 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_c4_a000d803c md5=3ae130e97b13e5bdfa4a7f9a14370b9f DECK=/resnick/groups/carnegie_poc/jingze/prod_c4_1124/decks/IA123_t9l62.athinput md5=f2e3b57d239f67e20b88ab89987e3323 start=-i /resnick/groups/carnegie_poc/jingze/prod_c4_1124/t9l62/deck.athinput WT=53:04:56 HOLDER=7333036 NODE=hpc-sm-01-15 EX=
```
- clean_init: ``
- Monitor: prod_c4_1124/watch.sh (c100, 500, ..., then a line every 10 min to status_t9l62.txt; posts and alert file at CW, DT_RETRY, fallback, clamped, FATAL, inf; FATAL / inf: clean STOP, no relaunch; wall-limit end: chain restart from the newest rst, at most 2); first rst hard-linked in rst_keep when it appears. Space: 42 T free.

- clean_init line (the start post above went out before it was printed, hence the empty field): `<rad_vet> clean_init steady: converged after 3 passes, residual 1.235e-09 (tol 1.0e-08; 0 formal solutions); radiation E 8.2441512239e+46 -> 8.2029653258e+46, F/c 1.5972e+42 0.0000e+00 0.0000e+00 -> 1.1510e+43 2.1812e+35 1.0049e+34 (volume sums)`; CLEAN WARNING 0, FATAL 0.

- (11-24 17:58) check c100: t9l62: cycle 100, s/cyc 2.3975 (c0-100), t 1.1138e+04 dt 96.85 L_top/L_in 1.00119 L_bot/L_in 1.00260, floor_pay(sum) 2.47e+27, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 18:15) check c500: t9l62: cycle 500, s/cyc 2.4891 (c200-500), t 4.9297e+04 dt 96.85 L_top/L_in 1.00334 L_bot/L_in 0.99528, floor_pay(sum) 1.22e+28, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 18:36) check c1000: t9l62: cycle 1000, s/cyc 2.5627 (c700-1000), t 9.7624e+04 dt 96.85 L_top/L_in 1.00172 L_bot/L_in 0.99600, floor_pay(sum) 2.60e+28, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 18:57) check c1500: t9l62: cycle 1500, s/cyc 2.4740 (c1200-1500), t 1.4634e+05 dt 96.84 L_top/L_in 1.00131 L_bot/L_in 0.99547, floor_pay(sum) 4.13e+28, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 19:17) check c2000: t9l62: cycle 2000, s/cyc 2.3950 (c1700-2000), t 1.9495e+05 dt 96.84 L_top/L_in 1.01115 L_bot/L_in 0.99567, floor_pay(sum) 5.83e+28, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 19:38) check c2500: t9l62: cycle 2500, s/cyc 2.5800 (c2200-2500), t 2.4270e+05 dt 96.84 L_top/L_in 1.01123 L_bot/L_in 0.99591, floor_pay(sum) 7.71e+28, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0
