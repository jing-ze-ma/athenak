# NOTE 2026-11-20dr orion -> caltech (USER): TGF_hg / TGF_rsg at the head; drop TGinf_* and TG4_rsg if not started

- New build **B_bfxf** = caltech-runs/beamfixf, commit 6ad6817aa5a124e1c18b95307f615ca5e05ad175.
- **TGF_hg and TGF_rsg** (in `arms_rush.txt`): thinf=true, thinf_gmax=-1 (the face-reconstruction cap), clean_maxit=800. NR=4, 1 node each, 100 cycles, no xcorner. They run concurrently at the **head of the queue**, with the same report items as TG4_hg.
- **Drop TGinf_hg, TGinf_rsg and TG4_rsg if they have not started:** g > 1 is rejected as non-realisable. If they are running, let them finish.
- XCP* and the v20.3 gate (NOTE dp) stay as queued.

-- orion (for the USER)
