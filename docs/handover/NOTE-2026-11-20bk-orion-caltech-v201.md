# NOTE 2026-11-20bk orion -> caltech (USER): v20.1 history-fix gate, 4 arms (new build v201), HIGH priority

New build **B_v201** = caltech-runs/v201, commit a3b72fdcc2793ed62c793ccca68dacea4b62d4c6 (v20n plus one commit: exact history sums; the state is unchanged). Plain cmake line.
Arms are in `arms_v201.txt`:
- **LI4b_rsg / LI8b_rsg**, exactly as LI4/LI8 (double bins, bins and final rst kept). Run `licmp.py` on LI4b vs LI8b, and also on LI4b vs LI4 (v20n),
  into `caltech_results/LI_cmp_v201.txt`. Expected: LI4b vs LI8b hst identical; bins identical to LI4.
- **V201_rsg / V201_hg**: rt NH=2 (romio). Xcheck vs V20N_*/a: give the max relative hst difference. Only round-off in sum columns is expected.

-- orion (for the USER)
