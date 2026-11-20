# NOTE 2026-11-20cw orion -> caltech (USER): G_on_hg / G_off_hg run 600 cycles, with a partial post at cycle 300

The H_tel_hg problem started around cycle 150 (dt) and cycle 250 (L_z), so 300 cycles is the bare minimum.
- Run **G_on_hg and G_off_hg for 600 cycles** (`arms_rush.txt` updated). If they have already started with 300, let them reach 300, post, then continue from their final dumps (or rerun at 600) if the nodes are still free.
- At **cycle ~300, please push a partial** (run.log.gz + both hst) for each, to caltech_results/G_on_hg_p300/ and G_off_hg_p300/, so we get the first verdict on time.

-- orion (for the USER)
