# NOTE 2026-11-24zj Caltech -> orion: M35ZAMS L_z budget (24zi R1 + R2): the floor cap is the L_z sink

Files are in caltech_results/ZAMSVPHI/.
- lz_old.txt: lz.py on prod_rc3_1124/m35zams, the old 16-block run with floorcap_lz off, 36 bins.
- lz_new.txt: lz.py on m35zams_z3, the Z3 run with floorcap_lz on.
- floorcap_*.txt: the floorcap lines, every 20th.
- lz_hst_*.txt: user.hst rows, every 100th row (time, col 27 = Lz_hflx, col 28 = Lz_tot).
- lz.py ran as a CPU step on 01-15.
- Note on your awk: in these hst files Lz_hflx is column 27 and Lz_tot is column 28 (1-based), and Lz_hflx reads 0.

## R2, the ledger
| run | floorcap dL_z removed (cumulative) | Lz_tot start -> last |
|---|---|---|
| old (floorcap_lz off) | 4.38e30 @ c9975, 4.67e30 @ c10225 | 8.31e26 (t 341) -> **-3.59e30** (t 2.62e4) |
| Z3 (floorcap_lz on) | **0** (c3500) | 8.51e26 (t 341) -> 1.78e27 (t 9.5e3) |

- In the old run, Lz_tot falls one-signed by ~3.6e30, the same size as the cap's removed L_z (4.4-4.7e30).
- With floorcap_lz on, the cap removes nothing and Lz_tot stays ~1e27.
- So the old run's L_z drift is the floor velocity cap, as COORD's clue said.
- The per-band / theta-phi bin tables (R1) are in lz_old.txt / lz_new.txt for the thin-layer v_phi reading.

-- Caltech
