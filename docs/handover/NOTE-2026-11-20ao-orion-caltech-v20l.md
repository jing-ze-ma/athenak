# NOTE 2026-11-20ao orion -> caltech (USER): v20l gate, 2 rt arms + speed pair, HIGH priority

New build **B_v20l** = caltech-runs/v20l, commit e067945c2cc1188ecab1e3283102ec66a1a97504. It is v20k (the lzcs-tel build, key default off) plus one CMake commit that
turns the H200 P5 device-memory bodies off by default. **Use the plain cmake line**: don't add any ATHENA_CUDA_BODY_GLOBAL or
EXTRA_DEFS. Arms: `docs/handover/caltech_runs_1121/arms_v20l.txt`:
- V20L_rsg and V20L_hg: rt NH=2, romio. Xcheck: hst == V20J_rsg/a and V20J_hg/a (all 30 user columns plus hydro).
- SPL_g_1 / SPL_l_1: together on 2 same-type nodes (not 02-16); then SPL_g_2 / SPL_l_2 with the nodes swapped. Report s/cycle and
  node; the hst must equal SPG_g_1 on cols 1-26.
Priority: right after T2K_rsg / T3K_hg, ahead of the remaining batch-4 / LEAN arms.

-- orion (for the USER)
