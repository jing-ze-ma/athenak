---
name: dhj-kink-horizontal-rad-candidate
description: Kink-cure candidate (user asked to note 09-24): horizontal radiative conduction already exists (<mhd>/rad_angular, off in prod4/hyd4) but acts only in the thick blend layers, while the kinks sit at p < 1e-7 bar
metadata:
  type: project
---

For the night-side kink decision (open since the 09-24 handover; options: resolution test (recommended), HLLE,
contact smoothing), the user asked to note horizontal radiative transport as a candidate.

Facts checked 09-24: the code already has HORIZONTAL radiative conduction: <conduction>/isotropic_conduction =
radiative with <mhd|hydro>/rad_angular (default true; implicit transverse options rad_implicit_ang / rad_sts_all /
rad_ang_solver in src/diffusion/conduction.cpp). prod4 and hyd4 set rad_angular = false (change 3 vs prod3; its dt term
was a concern). It acts only where the tau blend gives conduction weight (rad_tau_lo 30 / rad_tau_hi 300, thick deep
layers); the ck two-stream has no sideways transport at all.

**Why:** candidate cure for the kinks; but the kinks sit at p < 1e-7 bar (optically thin), where radiation cools
locally to space rather than diffusing sideways, so thick-layer horizontal diffusion is unlikely to reach them.

**How to apply:** if tested, a cheap A/B: rad_angular = true (implicit transverse) vs false from a production restart,
count kinks by pressure band; expect an effect only if kinks originate deeper. Related: [[next-prod-ck-c2]].

**RESOLVED 09-24 (real prod4 profile, /viper/ptmp2/jinma/hdiff_0924):** NOT a candidate. Where horizontal diffusion beats advection (p < 0.01-0.03 bar) cells are horizontally thin (diffusion invalid); where thick, advection (1e3-1e4 s) and vertical cooling are faster; blend weight is exactly 0 above column tau 30 (~0.3 bar), and at 1 bar t_diff,h = 500 rotations (128 grid) / 8.5 (1024 grid). At p<1e-7 bar rad_angular adds exactly zero flux. Kinks must be traced to the two-stream or the top boundary.
