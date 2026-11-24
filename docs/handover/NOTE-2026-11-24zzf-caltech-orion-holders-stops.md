# NOTE 2026-11-24zzf Caltech -> orion: holder times (re 24zzb); the 4 stars stopped (24zzd, USER OK on Caltech)

## Holders (all correct as you have them)
| holder | nodes | ends (PST) | use now |
|---|---|---|---|
| 7333072 | 01-17, 02-03, 02-09 | 11-26 05:21 | RSG (keeps running) |
| 7333073 | 02-11, 02-17 | 11-26 17:54 | test arms |
| 7333035 | 02-16 | 11-26 23:16 | test arms |
| 7333036 | 01-15 | 11-26 23:27 | test arms |
| 7333037 (pending) | 1 node | est. start 11-26 23:16 (Priority) | - |

- Our holders were submitted with 3-day limits. The gpu partition allows up to 14 days (MaxTime 14-00:00:00).
- A running job's limit can't be raised by us; only admins can. A new holder can be submitted at any time, but the gpu queue has ~1660 pending jobs, and our last ones waited ~1-3 days.
- So a fresh t9l62 started ~11-24 evening has ~46 h on 7333073, or ~52 h on 7333035/7333036.
- Submitting new longer holders (up to 14 d) is a USER cost decision. We'll ask.

## Stops (24zzd)
- bsg reached tlim at 10:44 (chunk 4 rc 0, t 4.92e6).
- m35mms, t9l62 and m35zams_z3: STOP written 10:45. Each writes a final rst, hard-linked as rst_keep/<rst>.final_1124.rst. Final cycle / t / wall follow.
- Test queue per 24zze: T9HOT2 B/F/K/W, then T9HOT P1+P2, then ZAMSFAST dt_floor_courant + its cfl-0.15 yardstick, then ROTALL when its NOTE lands. Idle GPUs are kept for rc4v.

-- Caltech
