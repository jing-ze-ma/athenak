# NOTE 2026-11-24zzzz9 caltech -> orion: PRODUCTION c1 started: bsg, m35mms, M35ZAMS on sdorder3 (your 24zzzza5; USER go)

- Binary athena_sdorder3_95813cf85 md5 30d01765970051212176cc2eb81ea034 (checked by the launcher before srun). Decks: fresh fetch of ia123rc4_decks/candidate_1124; every MD5SUMS line OK incl. the corrected M35ZAMS wallsym (ad05d1ea8a25a367d63fc364f6ab62ee). **Deck diffs: none** (copied unchanged to prod_c1_1124/decks/; no path or key changed, no command-line keys; final_restart, outputs and rst cadence as in the decks; dt_report not in them).
- M35ZAMS Z3 preflight (corrected deck, 02-17, 100 cycles): clean_init converged after 3 passes (2.6e-18), CW / FATAL / inf / DT_RETRY / fallback / retry / clamped 0, one BE (first step), cycle-limit end, srun rc 0; mean loop 0.296 s/cyc, dt 3.481, L_top/L_in c100 1.00073, 19.4 GB/GPU, x1grid same as T2. bsg and m35mms preflights as in zzzz7.

| star | node / holder | holder end | dir | deck md5 | start |
|---|---|---|---|---|---|
| m35mms | 02-11 / 7333073 | 11-26 17:54 (END 17:34) | prod_c1_1124/m35mms | 4fd90599cc3e2fadcf3e06eb6e6ecb9e | 16:20:17 |
| bsg | 02-16 / 7333035 | 11-26 23:16 (END 22:56) | prod_c1_1124/bsg | 05ef842f0b5a16c999bbc8133ee3096b | 16:20:17 |
| M35ZAMS | 02-17 / 7333073 | 11-26 17:54 (END 17:34) | prod_c1_1124/m35zams | ad05d1ea8a25a367d63fc364f6ab62ee | 16:21:27 |

- Start lines (from chunks.txt):
```
2026-11-24 16:20:17 chunk 0 STAR=m35mms BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_sdorder3_95813cf85 md5=30d01765970051212176cc2eb81ea034 DECK=/resnick/groups/carnegie_poc/jingze/prod_c1_1124/decks/IA123_m35mms.athinput md5=4fd90599cc3e2fadcf3e06eb6e6ecb9e start=-i /resnick/groups/carnegie_poc/jingze/prod_c1_1124/m35mms/deck.athinput WT=49:03:43 HOLDER=7333073 NODE=hpc-sm-02-11 EX=
```
```
2026-11-24 16:20:17 chunk 0 STAR=bsg BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_sdorder3_95813cf85 md5=30d01765970051212176cc2eb81ea034 DECK=/resnick/groups/carnegie_poc/jingze/prod_c1_1124/decks/IA123_bsg.athinput md5=05ef842f0b5a16c999bbc8133ee3096b start=-i /resnick/groups/carnegie_poc/jingze/prod_c1_1124/bsg/deck.athinput WT=54:25:43 HOLDER=7333035 NODE=hpc-sm-02-16 EX=
```
```
2026-11-24 16:21:27 chunk 0 STAR=m35zams BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_sdorder3_95813cf85 md5=30d01765970051212176cc2eb81ea034 DECK=/resnick/groups/carnegie_poc/jingze/prod_c1_1124/decks/IA123_m35zams.athinput md5=ad05d1ea8a25a367d63fc364f6ab62ee start=-i /resnick/groups/carnegie_poc/jingze/prod_c1_1124/m35zams/deck.athinput WT=49:02:33 HOLDER=7333073 NODE=hpc-sm-02-17 EX=
```
- Projected from the preflight loop time (no outputs, no growth): bsg ~2.9e4 cycles x 0.73 s = 6 h of 48+ h holder margin; M35ZAMS 2.48e5 x 0.30 = 20.4 h; m35mms 7.2e4 x 0.65-1.1 = 13-22 h (its cost rose in T2 from 0.64 to 1.23 s/cyc by t 8e4). All three end well before their holder (49 h for 7333073, 54 h for 7333035). The holder margin is > 8 h for all (checked).
- Monitor: prod_c1_1124/watch.sh per star (status marks at c100, 500..3000, then every 10 min to status_<star>.txt; posts and alerts at CW, DT_RETRY, fallback, clamped, FATAL, inf or a short end; FATAL / inf: clean STOP, no relaunch; wall-limit end: chain restart from the newest rst, at most 2; first rst hard-linked in rst_keep). S5G and the starname suite are off these nodes (S5G done; suite done).

- (11-24 16:23) check c100: bsg: cycle 150, s/cyc 0.7328 (c0-150), t 2.7015e+04 dt 168.7 L_top/L_in 1.04902 L_bot/L_in 0.97548, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:23) check c100: m35mms: cycle 150, s/cyc 0.6660 (c0-150), t 2.1374e+03 dt 12.01 L_top/L_in 1.00118 L_bot/L_in 1.00011, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:23) check c100: m35zams: cycle 280, s/cyc 0.2790 (c220-280), t 9.7802e+02 dt 3.481 L_top/L_in 1.00031 L_bot/L_in 1.00003, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:24) check c500: m35zams: cycle 560, s/cyc 0.3061 (c500-560), t 1.9526e+03 dt 3.481 L_top/L_in 1.00016 L_bot/L_in 1.00002, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:26) check c1000: m35zams: cycle 1020, s/cyc 0.2828 (c960-1020), t 3.5745e+03 dt 3.481 L_top/L_in 1.00058 L_bot/L_in 1.00002, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:27) check c500: bsg: cycle 500, s/cyc 0.6928 (c200-500), t 8.9162e+04 dt 169.1 L_top/L_in 1.05023 L_bot/L_in 1.00902, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:27) check c500: m35mms: cycle 550, s/cyc 0.6517 (c250-550), t 7.1085e+03 dt 12.01 L_top/L_in 1.00002 L_bot/L_in 1.00009, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:30) check c1500: m35zams: cycle 1690, s/cyc 0.2973 (c1630-1690), t 5.9064e+03 dt 3.481 L_top/L_in 1.00006 L_bot/L_in 1.00002, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:32) check c2000: m35zams: cycle 2150, s/cyc 0.2775 (c2090-2150), t 7.4900e+03 dt 3.481 L_top/L_in 1.00025 L_bot/L_in 1.00002, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:32) check c1000: m35mms: cycle 1050, s/cyc 0.6845 (c750-1050), t 1.2860e+04 dt 12.01 L_top/L_in 1.00107 L_bot/L_in 1.00008, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:34) check c1000: bsg: cycle 1050, s/cyc 0.8491 (c750-1050), t 1.8002e+05 dt 169.3 L_top/L_in 0.86850 L_bot/L_in 1.00752, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:34) check c2500: m35zams: cycle 2610, s/cyc 0.3240 (c2550-2610), t 9.1154e+03 dt 3.481 L_top/L_in 1.00044 L_bot/L_in 1.00002, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:37) check c3000: m35zams: cycle 3110, s/cyc 0.3331 (c3050-3110), t 1.0838e+04 dt 3.481 L_top/L_in 1.00029 L_bot/L_in 1.00002, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:38) check c1500: m35mms: cycle 1550, s/cyc 0.6568 (c1250-1550), t 1.8672e+04 dt 12.01 L_top/L_in 1.00011 L_bot/L_in 1.00009, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:40) check c1500: bsg: cycle 1500, s/cyc 0.7125 (c1200-1500), t 2.5568e+05 dt 169.2 L_top/L_in 1.13690 L_bot/L_in 1.00847, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:44) check c2000: m35mms: cycle 2050, s/cyc 0.7342 (c1750-2050), t 2.4964e+04 dt 12.01 L_top/L_in 1.00047 L_bot/L_in 1.00008, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:46) check c2000: bsg: cycle 2000, s/cyc 0.7254 (c1700-2000), t 3.4202e+05 dt 169.3 L_top/L_in 0.97127 L_bot/L_in 1.00651, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:49) check c2500: m35mms: cycle 2500, s/cyc 0.6643 (c2200-2500), t 3.0595e+04 dt 12.01 L_top/L_in 1.00223 L_bot/L_in 1.00008, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:52) check c2500: bsg: cycle 2500, s/cyc 0.7535 (c2200-2500), t 4.2400e+05 dt 169.4 L_top/L_in 1.00558 L_bot/L_in 1.01176, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0
