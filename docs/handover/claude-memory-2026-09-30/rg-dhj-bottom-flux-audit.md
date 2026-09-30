---
name: rg-dhj-bottom-flux-audit
description: 09-15 audit of red giant + dhj for the B-star bug classes (bottom flux on one twin path; RKL1 halo window) - none live; dhj has a LATENT double count of the internal flux if the tau blend is off (ck_int_at_cut default true + rad_flux_inner)
metadata:
  type: project
---

Mode 3 / rt_bottom_flux / RTCol3 exist only on rg-box-pcr (bench/wt_m3acc), not on implicit-transverse-raddiff or polar-average-perf, so bug (A) cannot be live in rg/dhj. Inner wall flux rad_flux_inner is added on exactly ONE path (conduction.cpp:1075-1079, explicit face is); the implicit x1 tridiagonal (conduction.cpp:1317) and the RKL1 stencil (conduction_transverse.cpp:109) exclude the physical faces by design - consistent. rg forces rt_int_at_cut=false and fatals if the blend is off (red_giant.cpp:981-987, 1018).
LATENT dhj: deep_hot_jupiter_rt.cpp:442 ck_int_at_cut defaults TRUE and :1580 sets rad_flux_inner = sigma T_int^4; with rad_tau_mode=false the internal flux would enter twice (Iint at is AND the wall). Fix = red giant's fatal at :1577 (or force int_at_cut=false when rad_flux_inner != 0). Dead legacy solvers with hard-coded Tint 100/500 at :2711, :2936 (never called) - delete.
Bug (B) window: rg/dhj refuse rad_implicit_ang on curvilinear meshes (conduction.cpp:227), so a8c69bc9 is irrelevant to them.
Mode-3 port to rg: shared two_stream_rt.hpp inherits the 1fdce2b0 fix; BUT the mode-3 gate is grey_on && rt_split (covers rg, NOT dhj's rt_ck path) and rt_bot_flux is wired only from box_convection.cpp (rg needs its own rad_flux_inner -> rt_bot_flux transfer, which zeroes rad_flux_inner). polar-average-perf is ahead on RT physics (rt_ali_diag, rad_x1_uform/kiter) vs implicit-transverse-raddiff -> merge order first. See [[he-mode3-flux-collapse]], [[port-rt-fixes-to-rg-dhj]].
