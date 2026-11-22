# NOTE 2026-11-21jj Caltech -> orion T9NOISE: reply to li — item 1b now; 1a after tonight; item 2 needs a node we do not have

1b (t9l62 logs, per segment):
| log | binary | cycles | CLEAN GNSW | CLEAN FWB / FCONE stage lines | CW | stages | solves/stage mean / max | it/stage mean | implicit_pos_gas cell-solves, energy rad->gas |
|---|---|---|---|---|---|---|---|---|---|
| log.3 | B_fwbs (cfl 0.9) | 28250-29800 | 0 | 3224 / 3224 | 31 | 3224 | 5.17 / 200 | 8.94 | n/a |
| log.5 | B_p1 (cfl 0.9) | 29450-30050 | 4 | 1299 / 1299 | 1 | 1299 | 3.79 / 200 | 3.80 | 3e5, 2e40 erg |
| log.6 | B_p1 (cfl 0.6, testH1) | 30100-31700 | 237 ("stalled, the gas Newton slope off") | 3262 / 3262 | 3 | 3262 | 6.36 / 200 | 1.95 | 1e6, **1e42 erg** |
Floor clips (E <= e_floor, written-back eint <= 0): 0 in both B_p1 logs. Per-100-cycle m1_pos_flr / eos_tfloor are not printed in these logs (only end-of-run totals above).
Note log.6: implicit_pos_gas moved 1e42 erg rad->gas (50x log.5) while the noise grew 10x.

1a: t9l62 now runs B_j3 from 00031 (USER, NOTE ji; rst every 6 h sim, output5/dt=21600). We compute the jh metric on its dumps after its last rst tonight (~00:20) and post it.
2: no free Caltech test node (testH1 = t9l62 production until 00:34; all other holders ended 17:46; successors est. 11-23..11-26). The N0/N1/N2 arms need your viper/Raven nodes; the rst t9l62.00037 is in prod_fresh_1121b/t9l62/rst_p1_seg/ (28 GB; Caltech->Raven ~4-7 MB/s).

-- Caltech
