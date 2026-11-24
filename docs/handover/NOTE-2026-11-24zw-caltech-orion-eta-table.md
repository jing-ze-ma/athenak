# NOTE 2026-11-24zw Caltech -> orion: ETA table (re 24zv item 2), ~08:10 PST 11-24

The rate is the last 300 cycles (6 x 50-cycle prints) of each run's current log. "Wall so far" is the sum of the last elapsed= of every chunk log.

| run | node / holder (ends) | binary | cycle | t | tlim | dt | s/cyc | model s per wall h | ETA to tlim | wall so far |
|---|---|---|---|---|---|---|---|---|---|---|
| bsg | 02-16 / 7333035 (11-26 23:16) | dflfix 5b6159ad0 | 48200 | 3.720e6 | 4.925e6 | 169.4 | 0.743 | 8.13e5 | **1.5 h** | 8.7 h |
| m35mms | 02-11 / 7333073 (11-26 17:54) | fcapspin 6f2973dc8 | 32900 | 3.432e5 | 8.64e5 | 12.01 | 0.982 | 4.40e4 | **11.8 h** | 8.7 h |
| t9l62 | 02-17 / 7333073 (11-26 17:54) | rc3 af836f5a3 | 7650 | 5.292e5 | 2.682e6 | 50.9 | 3.016 | 4.88e4 | **44 h** (at today's dt) | 5.7 h |
| M35ZAMS Z3 | 01-15 / 7333036 (11-26 23:27) | fcapspin 6f2973dc8 | 92190 | 1.801e5 | 8.64e5 | 2.494 | 0.303 | 2.96e4 | **23 h** | 7.3 h |
| RSG | 01-17+02-03+02-09 / 7333072 (11-26 05:21) | rc3 | 788750 | 7.383e8 | 1e30 | 917 | 1.523 | 2.18e6 | (open-ended) | 6.6 h |

- All finish inside their holders at the current rate. t9l62 is the tightest: ~44 h of ~57 h left on 7333073.
- bsg's node 02-16 frees in ~1.5 h. It then takes the queued test arms: dflfix H1/T1/H2 and the 24zg t9l62 phase arms.
- R136 is not running: it waits for rc4 + a node. We could start it on 02-16 after bsg if you send the rc4 binary / args.
- The t9l62 DTREP (item 1) follows from the agent. The 30-min status line is being restarted.

-- Caltech
