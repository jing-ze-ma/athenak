# NOTE 2026-11-20cg orion -> caltech (USER): BEAMFIX forms 12/11/5 (HG), new build B_bfx12

New build **B_bfx12** = caltech-runs/beamfix12, commit 8ba22894106c88e372bb7e6a28c0a7390b1f73d5. It is B_bfx7 plus forms 11, 12 and 5; form 0 is still v20.1 bitwise.
Three HG arms are appended to `arms_beamfix.txt`: BFX12_hg, BFX11_hg and BFX5_hg (NR=8, 300 cycles, the same layout as BFX4_hg).
**BEAMFIX order:** BFX7_hg, BFX12_hg, BFX11_hg, BFX5_hg, then BFX7_rsg, with the BFX0 controls as already scheduled. BFX5 is a reference that is expected to fail like form 4.
If an arm goes FATAL early, just post its line; no rerun is needed. RSG arms for the surviving forms will follow later.

-- orion (for the USER)
