# NOTE 2026-11-20cr orion -> caltech (USER): G_off_hg, the speed partner of G_on_hg (plain B_lztel3)

This firms up NOTE cq. G_on_hg runs on the plain **B_lztel3** (caltech-runs/lztel3 74c8a8bb1), with no ENVX. Add its key-off partner:
- **G_off_hg**: B_lztel3, HG, NR=8, 1000 cycles, key off (appended to `arms_lztel3.txt`). Its hst should equal L1_hg_off rows 0-999 bitwise. Please state whether it does.
- Run the pair **at the same time on the two 2-node slots, or back to back on the same nodes**, whichever starts sooner. Same node type and 2-node layout.
- **HG queue:** BFXC19_hg, G_on_hg, G_off_hg, BFXC20_hg, BFX16_hg, BFX15_hg, then the rest.
- Lines: s/cycle, mean dt, CW/FATAL, plus the bitwise verdict for G_off.

-- orion (for the USER)
