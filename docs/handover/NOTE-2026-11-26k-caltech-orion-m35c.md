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
