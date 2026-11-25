# NOTE 2026-11-26k caltech -> orion: M35C coarse m35mms (your 26k): launch and results

## Before the run
- Deck: production deck (4fd90599) with exactly this diff (deck 10590265):
```
25c25
< nx2                      = 256                    # [PAPER] S22: 256 theta cells
---
> nx2                      = 128                    # M35C (NOTE 26k): 128 theta cells (production 256)
30c30
< nx3                      = 256                    # [PAPER] S22: 256 phi cells
---
> nx3                      = 128                    # M35C (NOTE 26k): 128 phi cells (production 256)
45c45
< layout_tile = 128   # IA12 LAYOUT: = this deck's lateral MeshBlock (128 x 128); must divide meshblock nx2, nx3
---
> layout_tile = 64   # M35C (NOTE 26k): = this deck's lateral MeshBlock (64 x 64); production 128
48,49c48,49
< nx2                      = 128                    # 4 blocks = 1 per GPU on a 4-GPU node (Raven A100); viper 2 GPUs: 2 blocks/GPU
< nx3                      = 128                    # 
---
> nx2                      = 64                     # M35C (NOTE 26k): 2 x 2 blocks of 336 x 64 x 64, 1 per GPU (production 128)
> nx3                      = 64                     # M35C (NOTE 26k) (production 128)
56c56
< tlim = 8.64e5   # PROD 11-21 (orion ig/ih), replaced: 8.64e5
---
> tlim = 3.5e5   # M35C (NOTE 26k) (production 8.64e5)
337c337
< dt                       = 1800.0                 # 30 min (USER 11-23, Caltech 23x)
---
> dt                       = 7200.0                 # M35C (NOTE 26k): 2 h (production 1800)
341c341
< dt                       = 1800.0                 # 30 min (USER 11-23, Caltech 23x)
---
> dt                       = 1.0e30                 # M35C (NOTE 26k): off (production 1800.0)
347c347
< dt                       = 8.64e4   # 1 d (USER 11-24: rst cadence; was 21600.0)
---
> dt                       = 1.0e30   # M35C (NOTE 26k): only the final rst (production 8.64e4)
351c351
< dt = 43200   # FRESHGO (USER 11-20): every 4th output2 output (output2 dt 10800.0)
---
> dt                       = 1.0e30                 # M35C (NOTE 26k): off (production 43200)
355c355
< dt = 43200   # FRESHGO (USER 11-20): every 4th output2 output (output2 dt 10800.0)
---
> dt                       = 1.0e30                 # M35C (NOTE 26k): off (production 43200)
```
- Binary c4speed d979d5268 (3627c9a8), as the current production chunk; production's first 2.0e5 s (16340 cycles) were run by sdorder3 95813cf85, so the 256 run up to 2e5 s is another binary than M35C.
- Seed: he_seed_band = n2 is not a random field per cell: it is a sum of real Y_lm, l 16..64 (capped at 64; pi r / H_p 185), all m, fixed rng 1234, amplitude 1e-3, in the radial band r 8.98e11-9.35e11 (29 cells), the same radial band and l window on both grids (256 and 128 cells per wedge: 64 cells per l = 16 wave at 128, and l = 64 needs 22 cells per wave at 128 against 45 at 256). So the same smooth field sampled on coarser cells, not a new random realisation. Will be checked in the start-up prints (log line he_star_vet: seed = ...).
- mg levels stay 4: the code halves the tile 64 -> 32 -> 16 -> 8 -> 4 (x2, x3 only, tile's cells >= 4 and even); no reduction.
- Same: radial grid 336, box, rays GL 7 x 18 (vet_gd), opacity tables, wb_x1 arrays (radial), IC, dt_retry etc. Different by construction: lateral cell 2 x larger, blocks 336 x 64 x 64, layout_tile 64. Silent points checked: he_seed band, mg levels, the outputs (hydro_w 7200 s; vet, vet_vet, vet_face only at t = 0 and the end; rst only the final one; hst every cycle). dt to be reported.
- Size estimate: hydro_w 336 x 128 x 128 x 5 x 4 B = 110 MB each, 50 dumps: 5.5 GB; vet / face at 0 and the end about 1 GB; final rst about 7.5 GB (29 GB / 4 plus ghosts); total about 14-15 GB (38 TB free on the group volume).

## Running
- M35C started 10:47 PST on hpc-sm-02-03 (step 7333072.820, verified in squeue), 4 ranks, fresh directory /resnick/groups/carnegie_poc/jingze/m35c_1126/M35C; start-up: he_star_vet seed = real Y_lm l 16..64 (3969 terms, all m), amplitude 1e-3, rng 1234, r 8.9756e11-9.34835e11 (29 cells): identical to production's print, same warnings (only output1/dt not read, in production too); no FATAL, no CLEAN WARNING.
- Correction to my seed line above: in cells per wave l = 16 is 144 / 288 cells (128 / 256 grid) and l = 64 is 36 / 72 cells; the numbers I wrote there were wrong.
- After 200 cycles: s/cycle 0.200 (production 256: 0.64 at the start of its fresh run; 1.03 now), dt 12.1269 s (production 12.0076 at the start), estimate 28860 cycles x 0.20 = 5800 s = 1.6 h wall for 3.5e5 s; guard every 5 min (4 h, FATAL).
- Silent points: keys that differ from the defaults as production (no difference in that list); 4 blocks, layout_tile 64, mg levels 4; the dt differs by +1.0 % from production's start (lateral cells twice as wide), otherwise the same limiter.

## End (12:52:11)
- rc=0, Terminating on time limit at 3.5e5 s, 28980 cycles; wall 7475 s (2 h 5 min) in the main loop; CLEAN WARNING 0, FATAL 0, inf 0, DT_RETRY 0, stage fallbacks 0, backward-Euler steps 1 (the last, cut step).
- dt 12.126 s (12.1257 to 12.127; production 256 12.0066-12.0076); s/cycle 0.20 at the start, 0.26 mean over the run; dumps: 50 hydro_w (every 7200 s from 0 to 3.456e5, and the end), vet at 0 and the end; rst at t = 0 and the end only; 21 GB.
- Dumps: /resnick/groups/carnegie_poc/jingze/m35c_1126/M35C/bin/m35mms.hydro_w.*.bin, hst /resnick/groups/carnegie_poc/jingze/m35c_1126/M35C/m35mms.user.hst; grid file m35mms.x1grid.txt. The fixed readings A / B / C (band |k|, axis share, growth rate against the 256 run at 0.72, 0.74, 0.76 R_ph) go to the agent with the analysis scripts (a0dc2503); they follow in its NOTE.

- Correction: dt is 12.126 s in most cycles but dips: 1468 of 28980 rows have dt below 12.0 s (the two last, cut rows excluded, see the numbers below), so the range I gave before is wrong.
- Minimum 9.9619 s at t = 1.73000e+05 s; rows below 11.5 / 11 / 10 / 9 s: 949 / 601 / 20 / 0; DT_RETRY 0 (so the dt limiter, not a retry).
