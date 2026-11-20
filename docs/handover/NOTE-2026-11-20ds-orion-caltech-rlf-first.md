# NOTE 2026-11-20ds orion -> caltech (USER): R_LF on the FIRST free node (critical path)

R_LF is on the critical path: RA_on / RA_off need its 592690 dump, about 2.5 h after R_LF starts. The single-node order from now on:
1. **R_LF**, the v20.3 gate (V203L_rsg, V203L_hg, V203LP_rsg, V203_rsg, V203_hg) and TGF_hg / TGF_rsg, together on the first free nodes. R_LF takes the very first.
2. The 14 CONVERR arms (file order).
3. R_off, then LFa / L0a / LFb / L0b.
4. The EW pairs.
RA_on / RA_off start as soon as R_LF's 592690 dump exists. No preemption.

-- orion (for the USER)
