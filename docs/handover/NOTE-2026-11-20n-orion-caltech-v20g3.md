# NOTE 2026-11-20n orion -> caltech (USER): v20g3 replaces v20g2 (restart READ fix); 2 arms, top priority

New build **B_v20g3** = caltech-runs/v20g3, commit 35a5b37e683d963c051eca2582c111a3d25106db. It is v20g2 plus one commit. v20g2's 1 GiB pieces can make the HG
restart read-ahead a single slab, and then a buffer that was never resized FATALs at the restart READ (leg b of the rt arm).
v20g3 sizes that buffer to one slab per MeshBlock from the start.
- Arms: `docs/handover/caltech_runs_1121/arms_v20g3.txt`, V20G3_hg and V20G3_rsg, rt NH=2, **default io (romio321)**. V20G3_rsg
  hst must equal V20G_rsg; please add that xcheck.
- If V20G2_* has not started yet, drop it and run V20G3_* instead. If V20G2_* already ran, send its result anyway.

-- orion (for the USER)
