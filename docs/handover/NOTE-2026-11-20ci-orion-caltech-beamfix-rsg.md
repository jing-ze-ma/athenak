# NOTE 2026-11-20ci orion -> caltech (USER): BEAMFIX: drop BFX7_hg, add RSG arms of forms 12/11/5

- **Drop BFX7_hg**: form 7 died on RSG (BFX7_rsg FATAL at cycle 50). If it is already running, let it finish.
- The HG order is BFX12_hg, BFX11_hg, BFX5_hg.
- **New arms:** BFX12_rsg, BFX11_rsg and BFX5_rsg (B_bfx12, NR=4, 300 cycles). They go to the head of the single-node pool now, in parallel with the HG arms. They are appended to `arms_beamfix.txt`.
- Please post BFX4_rsg's result line when it is done.

-- orion (for the USER)
