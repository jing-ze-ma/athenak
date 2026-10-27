# TASK viper -> deltaai (USER): after probe 3466043, test the Hsml race fix hsmlfix-1026, ONLY if the probe reproduces the bad Hsml

**Bundle:** the user copies it to `/work/nvme/bivj/jma20/arepo2_gpu_1025/bundle/`:
- `arepo_hsmlfix_1026.bundle`, md5 **f547eaaf285bd307c72a32c5483801e2**.
- Head: `hsmlfix-1026` = **2d3967b136**. Prerequisite: **d9c0152bce**, which you have.
- Run `git bundle verify` first. If the bundle is missing or verify fails, stop and report.
- The bundle has 3 commits:
  - **9307d33d49**, the suspected fix: the top-leaf exchange no longer does `MPI_Allgatherv` straight into the node-shared buffer. Each rank gathers privately, then the island master publishes (fenced).
  - The env-gated diagnostic `AREPO_HSML_DIAG`. `=1` is the same print as your probe binary; `=2` also prints every boundary cell.
- Raven gates all PASS: g25_full init bitwise at 256/192/128 ranks; IDORT CPU vs GPU ALL SAME.

**Keep probe 3466043 as it is** (d9c0 + print, md5 31b1e582). Send its NOTE first.

**Then decide on what it shows:**
- **The probe shows cells with Hsml > 0.2 box, or a VCELL TERMINATE:**
  1. Build g2m from 2d3967b136 with the same recipe as the probe binary.
  2. Run one job in your interactive slot with the probe's layout: 2 GH200 nodes; 576 ranks first, then 256 if time allows.
  3. Set `AREPO_HSML_DIAG=2`. Each run stops at the first Sync-Point (PASS) or TERMINATE (FAIL).
     - Do **not** set `MPICH_SHARED_MEM_COLL_OPT`. The point is to see whether the code fix alone cures it.
  4. Keep it to 20 min or less (about 3-6 SU). Exclude gh082 and gh011.
- **The probe is clean (no such cells, both runs PASS):** do not run the fix. Report only.
  The bug is then not on DeltaAI, and Delta answers alone.

**Reply NOTE:**
- PASS/FAIL per run.
- `HSMLDIAG` counts of cells with Hsml > 0.2 box, and their IDs/Hsml.
- The density-loop iteration count.
- SU used.

Then stop. Short-test rules apply: 1 job, 2 nodes or fewer, 1 h or less. job3 stays cancelled; it is resubmitted only on viper go + USER go.
No Arepo source in git.
