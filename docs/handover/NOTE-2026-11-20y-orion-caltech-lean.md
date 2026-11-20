# NOTE 2026-11-20y orion -> caltech (USER): LEAN, 3 arms, LOWEST priority (fill-in)

New build **B_lean** = caltech-runs/lean, commit 425f2548fa3d6880f523c6c6b596bbb4ff2ca877. It is v20g with the cell-centred clean layout removed; the staggered path
should give bitwise the same results. Arms: `docs/handover/caltech_runs_1121/arms_lean.txt`:
- LN_rsg: rt NH=2.
- LNP_rsg: poison.
- LN_hg: rt NH=2 with **ompio**, because this base has no 1 GiB rst fix.
Xchecks to add to the result lines:
- LN_rsg/a, LNP_rsg hst == V20G_rsg/a.
- LN_hg/a hst == V20G_hg_ompio/a.
- rt within each arm.
If a run stops with "### FATAL ... belonged to the removed cell layout", send the line (it names a key) and skip that arm.
Run these only when nothing else is waiting; they are lowest priority.

-- orion (for the USER)
