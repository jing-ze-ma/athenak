# TASK 2026-11-17 viper -> Orion (USER): stag2 production-state tests on HG and RSG, in the idle Raven pHG262 holder (24 A100, now until ~05:00 PST)

USER: use the idle Raven allocation for staggered-scheme (stag2) tests on the HG and RSG production states.

## The allocation
- **Raven job 32307709 (pHG262), 6 nodes x 4 A100**: ravg1031, 1134, 1135, 1206, 1207, 1211. Running now and idle.
- **Hard stop: every step must END by 05:00 PST 11-17 (14:00 CET).** RSG production moves into this job when its current link in pHG
  32319677 ends (~05:34 PST).
  - Our watcher /raven/ptmp/jinma/prod_1107/jobs/rsg_in_phg_1116.sh waits until the job has no step other than batch/extern, then probes
    and links RSG.
  - Steps still running then delay RSG production. Never leave anything running past 05:00 PST.
- **Use overlap steps only:** `srun --jobid=32307709 --overlap ...`, step names `orion_*`, A100 recipe as in NOTE 8077e78a6
  (gcc/13 cuda/12.6 openmpi_gpu/5.0, -ffp-contract=off so --fmad=false).
- **Do NOT touch:**
  - /raven/ptmp/jinma/hg262_1113/STOP. It holds the job; never remove it and never write FIXED there.
  - the RSG claim files or the watcher;
  - pHG 32319677 (RSG production);
  - GPU 0 of pool 32367525 (W189 twin).

## States (read only; hard-link or copy, never write into production dirs)
- **RSG E (cubed sphere, 24 ranks x 4 per node in production):**
  - newest rst /raven/ptmp/jinma/raven_gap_1105/rsgE/rst/rsg.01592.rst (production continues writing there, so link a stable one, not the
    newest being written);
  - production keys /raven/ptmp/jinma/raven_gap_1105/rsgE/run.cfg; binary prod-1109 864d120d.
- **HG (old He giant, cubed sphere, viper production, 12 MI300A):**
  - newest rst is on viper /viper/ptmp2/jinma/prod_1107/hg/rst/hegiant.00171.rst (~72 GB);
  - Raven compute nodes can't see /viper, so copy with `dd oflag=direct` from the viper login and check md5;
  - an older Raven HG dir exists, /raven/ptmp/jinma/raven_gap_1105/hg; check whether its rst is recent enough;
  - keys /viper/ptmp2/jinma/prod_1107/hg/run.cfg.

## What to test (your design; suggestions)
1. stag2 (`clean_layout = staggered`, your recommended set) restarted from each production state, short runs (tens of cycles):
   - stability (0 CLEAN WARNING / FATAL);
   - s/cycle and peak VRAM vs production on the same GPUs and decomposition;
   - the first-cycle transient vs production.
2. The same with the cell layout (clean-cc 5785ee51d + the restart fix 5eb7adf3f, see NOTE 86a83eced) for a three-way comparison, if time
   allows.
3. Report as NOTE-2026-11-17-orion-viper-stag2-hg-rsg.md: per state, scheme vs production s/cycle, VRAM, warnings, and any blow-up.

Note: the clean restart fix 5eb7adf3f (stale vgd_cs below vlat_icut) matters for any clean run restarted from a vet_gd production state.
Take it if stag2 reads vgd_cs in the FV path.
