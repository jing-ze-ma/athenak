# NOTE 2026-11-20m orion -> caltech (USER): CONSERV energy-remainder ablation, 5 short RSG arms (6 cycles each)

All five arms use your existing **B_tie-prof** build (2ae22f2c4) with the same setup as TQS_rsg_tie, plus one key each.
Arms: `docs/handover/caltech_runs_1121/arms_conserv_ab.txt`. Priority within the group: AB1, AB3, AB2, then AB4 / AB5.
Place the group after the v20g2 arms and the CONV RSG/HG trio, ahead of batch 2/3.

AB3 uses `rad_m1/dbg_hydro_skip_t2=true`. If the command line refuses that key (unknown parameter), say so in its result line and
skip it. We will send another way to set it.

Send back as before: run.log.gz plus the user/hydro hst files.

-- orion (for the USER)
