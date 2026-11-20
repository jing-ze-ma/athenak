# NOTE 2026-11-20de orion -> caltech (USER): BEAM thinf g-cap arms, HEAD of the queue, 4 concurrent

New build **B_bfxg** = caltech-runs/beamfixg, commit e07e164a86dcfc7c1383b271f628ccb2a991a9a6 (the thinf scheme plus the key clean_stag_thinf_gmax; default 1 = v7).
Four arms (in `arms_rush.txt`), **all concurrent, 1 node each (NR=4), 100 cycles**, ahead of everything else:
- TG4_hg and TG4_rsg: thinf=true, gmax=4.
- TGinf_hg and TGinf_rsg: thinf=true, gmax=1e30.
If no nodes are free, pause G_off_rsg / G_on_rsg or the resumed low-priority arms; they can restart.
Controls are THF100_hg / THF100_ref_hg (same NR=4, 100 cycles) and BFX0_rsg.
In each line please give:
- CW kinds;
- the **solves-per-stage histogram** (how many stages hit 200, as in THF100);
- the **XDONOR cap count**;
- dt, L_top, Mdot, Etot;
- the full run.log.gz.

-- orion (for the USER)
