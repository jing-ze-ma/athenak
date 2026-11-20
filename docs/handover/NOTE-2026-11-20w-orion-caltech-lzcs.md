# NOTE 2026-11-20w orion -> caltech (USER): LZCS, 2 short arms (new build)

New build **B_lzcs** = caltech-runs/lzcs, commit 8ea56557cc313f184244d9d3742203ed578e3548 (v20g + lzbook + 2 output-only study commits).
Arms: `docs/handover/caltech_runs_1121/arms_lzcs.txt`. L1_rsg and L2_hg run 6 cycles each, env LZCS=1, with the same extra keys
as CS_*_ctl. They need about 1 GB/rank more GPU memory than CS_*_ctl (HG was at 122 GB/GPU, so it should still fit).
Priority: after AB7/AB8, before batch 2/3. Send back run.log.gz plus the hst files.

-- orion (for the USER)
