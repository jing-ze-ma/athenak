# NOTE 2026-11-20cq orion -> caltech (USER): G_on_hg on the plain fix binary (if not started)

If G_on_hg has **not started**, run it on a new build, **B_lztel3** = caltech-runs/lztel3, commit 74c8a8bb1cf0a3dfc237db1cdc82f9a8ff5ef40c. It has the same tree as the fix and no LZ3 prints. Use no ENVX (no LZ3, no ATHENA_DT_REPORT) and the H_tel_hg / L1_hg_off layout (NR=8 on 2 nodes). This gives a clean speed number along with the gate. The arm line is updated in `arms_lztel3.txt`.
**If the rebuild would delay G_on_hg, keep the study binary B_lztel3f as posted;** the gate still holds. If it has already started, let it run.
Queue position is unchanged: right after BFXC19_hg.

-- orion (for the USER)
