# NOTE 2026-11-20cx orion -> caltech (USER): v20.2 GPU gate (6 short arms)

New build **B_v202** = caltech-runs/v202, commit f4cefb18115b438349f293faffdce5a4a5965308. It is v20.1 plus the theta-wall fix and the rhse vcs_max MPI_MAX fix; both should be bitwise for production decks.
The arms are in `docs/handover/caltech_runs_1121/arms_v202.txt` (4 cycles each, the V201D recipe and restarts):
- V202_rsg (rt), V202P_rsg (poison), V202_hg (rt): hst should equal V201D_*/a bitwise.
- V202F_rsg (rt), V202FP_rsg (poison), V202F_hg (rt): the same, plus hydro/floor_pay=rad_column.
**Priority:** after the BEAMFIX arms and G_on/G_off (NOTE cv), ahead of LF/LX and everything else. They are short, so take any free single node.
Please report, as for V201D:
- the rt verdict;
- poison vs plain;
- hst vs V201D_*/a (the V202_* arms);
- CW / FATAL;
- s/cycle.

-- orion (for the USER)
