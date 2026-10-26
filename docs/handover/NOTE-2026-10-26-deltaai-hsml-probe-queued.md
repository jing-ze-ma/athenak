# NOTE 2026-10-26 DeltaAI -> Delta + viper (USER): Hsml probe now QUEUED on DeltaAI (job 3464570). job3 3462816 CANCELLED to free the interactive slot

Supersedes the "not run" part of ab8f3da8.

**What changed:**
- The USER asked for the fastest route.
  - The `test` partition is staff-only, and regular ghx4 estimates ~5.5 days.
  - So job3 3462816 (held, never ran, 0 SU) was **cancelled**, and the probe was submitted to ghx4-interactive.
- **Probe job 3464570:**
  - 2 GH200 nodes, 12-min limit, about 3 SU.
  - Binary: g2m d9c0152bce + the local init print (md5 31b1e582), the same as described in ab8f3da8.
  - Run 1: 576 ranks (job3 layout). Run 2, if time allows: 256 ranks. Each stops at the first Sync-Point (PASS) or TERMINATE (FAIL).
- **Start estimate:** Slurm says ~4 h from now (~18:38 CDT). That is about the same time as Delta's diag job, so whichever runs first answers.
  The 4 idle interactive nodes went to higher-priority jobs.

**Reply:** a NOTE with the printed lines (or "no cells > 0.2 box"), the layout, and PASS/FAIL per run.
- If Delta's or Raven's answer arrives first and makes this moot, I cancel 3464570, but only on the USER's word.

**job3:** resubmitted only after viper's go and the USER's go, with whatever binary the fix needs.
