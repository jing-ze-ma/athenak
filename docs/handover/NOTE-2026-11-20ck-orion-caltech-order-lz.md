# NOTE 2026-11-20ck orion -> caltech (USER): new order, the L_z fix moves up

The USER says the L_z fix matters as much as the beam forms. The order from now on:
1. BEAMFIX forms 13 / 12 / 11 / 5, on both stars (as queued).
2. **LZTEL3 G_on_hg (1000 cycles, NR=8) and G_rsg (1000 cycles, NR=4) on B_lztel3f.** Run them **in parallel with the beam forms on any free node**, ahead of everything else. G_rsg may already be running; keep it.
3. LF_* / LX_* long single-key runs.
4. D_off_hg.
5. LEAN2 rows (LN2_*, V201Dk_*).
6. The L2_* seeds.
No preemption of running arms. Keep the 2-node reservation from 11-21 13:00 PST.

-- orion (for the USER)
