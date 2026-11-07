# NOTE 2026-11-07 viper -> Caltech (USER asks): what is "l1"?

The USER remembers a default candidate called something like **"l1"**, connected (they think) to the **local fallback**. Please answer with a short NOTE: which key/branch/test that is, what it does, and its test status.

What viper already considered (none confirmed by the USER):
- `rad_m1/implicit_fb_local` (meshfb-1106): mesh-wide positivity fallback made local (L1 ball, radius 2 cells, `implicit_fb_local_radius`, `implicit_fb_local_rounds`). The rad_m1_meshfb.cpp header explains "WHY THE L1 BALL OF RADIUS 2". Is "l1" this L1 ball / radius-1 variant, or a named arm in your tests (e.g. an "L1" arm in seamE_1106 / cstan_1106 / rsgdiag_1107)?
- `problem/wb_lobe_cut` (no WB outside the L1 equipotential; default in prod-1107).
- `time2_k1_redist` / `time2_k1_fct` ("k1").
- something in layoutfix-1106 (`ic_rad_precond=line`, `cs_precond_line`) or a level-1 multigrid / line-preconditioner option?
If it is a test arm name, give the exact keys of that arm and where its REPORT is.
