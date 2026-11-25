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

- (11-24 16:56) check c3000: m35mms: cycle 3050, s/cyc 0.6991 (c2750-3050), t 3.6923e+04 dt 12.01 L_top/L_in 1.01155 L_bot/L_in 1.00008, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 16:59) check c3000: bsg: cycle 3050, s/cyc 0.8296 (c2750-3050), t 5.1854e+05 dt 169.4 L_top/L_in 0.97028 L_bot/L_in 1.00739, floor_pay(sum) 0.00e+00, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

## m35mms c1 (sdorder3, candidate keys, mg 4) vs the T2 run (mg 2) at the same model time, read-only (11-24 17:17)
c1: t 5.71791e+04, 4822 hst rows; T2 run to t 7.97864e+04
| model time [s] | c1 cycles | dt min/mean c1 | dt min/mean T2 | s/cyc c1 | s/cyc T2 | its t2s2/t2s3 (solves) c1 | T2 | L_top/L_in c1 / T2 | L_bot/L_in c1 / T2 |
|---|---|---|---|---|---|---|---|---|---|
| 0-6000 | 0-500 | 12.01/12.01 | 12.01/12.01 | 0.655 | 0.638 | 8.2 (2.00) / 9.8 (2.15) | 9.2 (2.00) / 10.7 (2.17) | 1.0001 / 1.0001 | 1.00010 / 1.00010 |
| 6000-12000 | 500-1000 | 12.01/12.01 | 12.01/12.01 | 0.661 | 0.641 | 7.9 (2.01) / 10.6 (2.28) | 8.7 (2.00) / 11.9 (2.26) | 1.0010 / 1.0010 | 1.00008 / 1.00008 |
| 12000-18000 | 1000-1500 | 12.01/12.01 | 12.01/12.01 | 0.672 | 0.627 | 8.2 (2.00) / 10.8 (2.39) | 8.7 (2.00) / 9.7 (2.10) | 1.0003 / 1.0003 | 1.00009 / 1.00009 |
| 18000-24000 | 1500-2000 | 12.01/12.01 | 12.01/12.01 | 0.710 | 0.694 | 8.1 (2.00) / 15.5 (2.75) | 9.0 (2.00) / 9.7 (2.09) | 1.0003 / 1.0003 | 1.00008 / 1.00008 |
| 24000-30000 | 2000-2500 | 12.01/12.01 | 12.01/12.01 | 0.682 | 0.627 | 8.2 (2.00) / 11.8 (2.53) | 10.0 (2.00) / 9.1 (2.04) | 1.0019 / 1.0019 | 1.00008 / 1.00008 |
| 30000-36000 | 2500-3000 | 12.01/12.01 | 12.01/12.01 | 0.714 | 0.658 | 8.8 (2.00) / 14.7 (2.86) | 10.9 (2.00) / 11.9 (2.40) | 1.0091 / 1.0091 | 1.00008 / 1.00008 |
| 36000-42000 | 3000-3500 | 12.01/12.01 | 12.01/12.01 | 0.728 | 0.716 | 10.8 (2.00) / 16.0 (2.78) | 12.1 (2.00) / 19.2 (2.81) | 1.0214 / 1.0213 | 1.00008 / 1.00008 |
| 42000-48000 | 3500-4000 | 12.01/12.01 | 12.01/12.01 | 0.710 | 0.792 | 11.5 (2.00) / 12.2 (2.56) | 13.1 (2.00) / 17.0 (2.86) | 1.0099 / 1.0099 | 1.00008 / 1.00008 |
| 48000-54000 | 4000-4500 | 10.69/11.93 | 10.71/11.93 | 0.774 | 0.747 | 13.6 (2.30) / 13.1 (3.20) | 15.3 (2.30) / 14.6 (3.15) | 1.0471 / 1.0471 | 1.00007 / 1.00007 |
| 54000-57179 | 4500-4800 | 9.10/9.92 | 9.11/9.92 | 0.809 | 0.798 | 18.2 (3.46) / 9.5 (2.51) | 22.8 (3.53) / 11.4 (2.42) | 1.0891 / 1.0890 | 1.00002 / 1.00002 |
c1: first dt < 11.9 at t 5.3398e+04; T2: 5.3410e+04
- **The dt drop is the flow's, not the candidate's: it starts at the same model time** (first dt < 11.9 at t 5.340e4 in c1, 5.341e4 in T2) and follows the same values (4.8-5.4e4: 10.69 / 10.71 min; 5.4-5.6e4: 9.10 / 9.11). L_top/L_in agree to 1e-4 at every bin end.
- s/cyc is the same within a few % at equal model time (0.77 now); Krylov its t2s2 17 vs 22 (c1 has mg 4, T2 mg 2), passes 3.07 vs 3.31.
- **The limiting cell is not known here:** dt_report is off in both runs at this time (no DTREP lines; c1's deck has none, T2 dropped it after its first rst), so the cause of the limit is not measured and I do not infer it. A short side arm with hydro/dt_report=10 from a c1 rst would name it, only if you want it; nothing is run on the production.
- T2 after this point, for the projection: dt 8.6-12.0 with mean ~10.6, s/cyc 1.0-1.25 (t 6e4-8e4), see zzzz4.

- (11-24 22:23) **bsg finished** (t 4.9248000000000000e+06 of 4.9248e6): bsg: cycle 29650, s/cyc 0.7227 (c29350-29650), t 4.9248e+06 dt 163.6 L_top/L_in 1.00774 L_bot/L_in 1.00597, floor_pay(sum) 5.94e+17, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-25 12:15) **m35mms finished** (t 8.6400000000000000e+05 of 8.64e5): m35mms: cycle 73250, s/cyc 0.9436 (c72950-73250), t 8.6400e+05 dt 5.353 L_top/L_in 1.02089 L_bot/L_in 1.00070, floor_pay(sum) 1.54e+25, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-25 12:27) **m35zams finished** (t 8.6400000000000000e+05 of 8.64e5): m35zams: cycle 248240, s/cyc 0.5097 (c248180-248240), t 8.6400e+05 dt 0.4957 L_top/L_in 0.99996 L_bot/L_in 1.00001, floor_pay(sum) 1.43e+21, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0.000000e+00
