# NOTE 2026-11-20s orion -> caltech (USER): v20j GPU gate, 4 arms, top priority

New build **B_v20j** = caltech-runs/v20j, commit 6b53883159f696bca903674103fc8dd217608907. It is v20g3 plus the hst budget columns (lzbook) and the H200 speed
commits (h200-p4 / h200 from batch 4, plus two more). Same cmake line, without the MB2 define.
Arms: `docs/handover/caltech_runs_1121/arms_v20j.txt` (V20J_rsg rt NH=2, V20JP_rsg poison, V20JR_rsg repro, V20J_hg rt NH=2).
Default io (romio321), right after the V20G3 arms.

Please add these xchecks to the result lines:
- V20JP_rsg and V20JR_rsg hst == V20J_rsg/a, bitwise.
- **user hst columns 1-26 and the hydro hst == V20G_rsg/a (RSG) and V20G_hg_ompio/a (HG)**, bitwise. v20j appends 4 user
  columns at the end (M_hflx E_hflx Lz_hflx Lz_tot), so compare only the first 26.
- s/cycle next to V20G's (2.86 RSG / 4.20 HG). The speed commits should cut host syncs.

-- orion (for the USER)
