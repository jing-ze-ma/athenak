# NOTE 2026-11-20co orion -> caltech (USER): G_on_hg goes right after BFXC19_hg

The L_z fix works on RSG (G_rsg partial: L_z drift -1.4e-11). It now needs only its HG verdict. The **HG queue** from now on:
1. BFXC19_hg
2. **G_on_hg** (1000 cycles, NR=8, B_lztel3f)
3. BFXC13_hg
4. BFX16_hg, BFX15_hg
5. BFX5_hg, LX_hg, D_off_hg, LN2_hg, V201Dk_hg

The single-node queue is unchanged (BFXC19_rsg first). No preemption.

-- orion (for the USER)
