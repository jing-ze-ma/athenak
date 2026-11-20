# NOTE 2026-11-20cl orion -> caltech (USER): BEAMFIX diagnostic forms 16/15 (HG); drop BFX11_hg if not started

New build **B_bfx16** = caltech-runs/beamfix16, commit 4ed45827f4198af91d0d1b25932e5e8fc073d690. It is B_bfx13 plus forms 15 and 16.
- **New HG arms:** BFX16_hg, then BFX15_hg (NR=8, 300 cycles; 60 cycles is enough if slots are tight). A FATAL by step ~12 is the signal.
- **HG order:** G_on_hg keeps its place (USER), then BFX16_hg, BFX15_hg, BFX5_hg, then the rest as in NOTE ck.
- **Drop BFX11_hg if it has not started.** BFX13_hg runs if it is already going. The RSG arms of forms 11/13/5 stay.

-- orion (for the USER)
