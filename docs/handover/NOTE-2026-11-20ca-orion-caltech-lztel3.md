# NOTE 2026-11-20ca orion -> caltech (USER): LZTEL3 cs_lz_tel study, 3 HG arms (+1 RSG if a slot is free), LOW priority

Two new builds:
- **B_lztel3s** = caltech-runs/lztel3-study, commit 3b22ddd2fa776df474269221726c3de6e19723dd (v201e plus env-gated LZ3 prints; no state change).
- **B_lztel3f** = caltech-runs/lztel3-fixstudy, commit 90c4715106f10abf520a5b2a6b767d88e00bd930 (candidate cs_lz_tel fix plus the prints).

Arms are in `docs/handover/caltech_runs_1121/arms_lztel3.txt`: D_on_hg 500, D_off_hg 300, G_on_hg 1000 (same layout as H_tel_hg, NR=8 on 2 nodes), and G_rsg 1000 (NR=4) only if a slot is free.
**Priority 5:** run them after BEAMFIX, the attribution set, LEAN2 and the L2 seeds. Do not preempt anything.
Results go to caltech_results/<TAG>/ as usual. Please include the full run.log.gz, because we grep the `^LZ3` and dt-report lines, plus the hst files.
Provenance: D_on_hg hst should equal H_tel_hg rows 0-499 bitwise, and D_off_hg should equal L1_hg_off rows 0-299.

-- orion (for the USER)
