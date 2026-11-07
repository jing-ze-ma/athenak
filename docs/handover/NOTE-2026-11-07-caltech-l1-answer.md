# NOTE 2026-11-07 Caltech -> viper: what "l1" is

"l1" is a **test arm name** in seamE_1106 (HG real state). It is not a separate key. Its local part is `implicit_fb_local`, which is the "local fallback" connection.

- **r1** (binary sd10) = `rad_m1/cs_cross_mono=true`, `rad_m1/implicit_opac_newton=false`, `mesh/cs_seam_rho_guard=4`, `rad_m1/implicit_predictor_tcap=1.1`, `rad_m1/time2_k1_redist=0.1`.
- **l1** (binary sd12) = r1 + `rad_m1/implicit_fb_local=true`, `rad_m1/implicit_fb_local_rounds=4`. The radius is the default, 2 cells (the L1 ball in rad_m1_meshfb.cpp; why the arm is called "l1" was not recorded).
- Arm l1 ran from HG rst 139 (meshfb_1106/rst/hegiant.00139.rst) to tlim 50284800 on 6 H200, with cfl 0.3. The full command line is in `$J/seamE_1106/runs/l1/arm.info`; the base keys are the HG prod set listed there.
- Later, **r7** (sd19) is the same key set on a newer build (donor-order fix). It agrees with l1 to 0.05 %. So "r7 keys" = "l1 keys" = the r7 default candidate in the handover.

**Status (AUDIT `$J/seamE_1106/AUDIT.md`, REPORT `$J/seamE_1106/REPORT.md`)**
- l1 vs control:
  - NOT POSITIVE 12 (to t0+130000 s) vs 321 (to t0+88000 s);
  - energy residual/L <= 7e-6 vs up to 1.3e-3;
  - vmax 135 vs 543 km/s;
  - vlat fallback per cycle 0.002 vs 1.5.
- l1 was the key-off reference for the cs full-tensor HG arm (cstan_1106/REPORT.md).
- fb_local adds nothing to L_top above noise. The +1 % L_top comes mostly from cs_cross_mono, and redist adds +0.3 %.
- Still open before default:
  - layout excess of guard 4 / Newton off (opac lag; passes=3 cuts it ~9x);
  - fb_local layout and restart gate (no layout case existed);
  - recommendation: keep time2_k1_redist opt-in;
  - on RSG, fb_local raised the floor count; not explained (rsgdiag_1107).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
