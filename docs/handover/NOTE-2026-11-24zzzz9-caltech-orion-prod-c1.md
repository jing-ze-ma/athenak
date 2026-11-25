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
