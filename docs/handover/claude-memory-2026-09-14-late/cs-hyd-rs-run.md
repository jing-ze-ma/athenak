---
name: cs-hyd-rs-run
description: cs_hyd_rs (bench/cs_hyd_rs/{hllc,lhllc,ausmpup}), the PURE-HYDRO twin of cs_mhd_prod2 run three times with different Riemann solvers to see whether the jet is the same (launched 09-07 08:25, 3 h leads before the maintenance, chained 24 h jobs after); f6f0d6f4 GPU binary
metadata:
  type: project
---

**What:** bench/cs_hyd_rs/base.athinput = cs_mhd_prod2's input with `<mhd>` -> `<hydro>`
(MHD-only params dropped), bbot = 0, output hydro_w; everything else (cs 128x32x32,
stretch, WB polytropic + rot_potential, radiative diffusion tau blend + ck table, floors,
general EOS, rst every 2 rot) identical. Arms: hllc; lhllc; ausmpup with ausm_mcut_p = 1
and ausm_wall_hllc = true. CPU smoke (16x16 panels, 20 cycles) passed for all three.
Jobs: hllc 11529253 (lead 3 h) -> 11529254 -> 11529255; lhllc 11529256 -> 57 -> 58;
ausmpup 11529259 -> 60 -> 61; ppmx (hllc + reconstruct = ppmx, nghost 3; horizontal order 4 on cs, radial stays PLM on the stretch) 11529330 -> 31 -> 32, launched 08:40. All three solvers go through the same curvilinear face
machinery (PrimFace rotations, gnomonic flux correction, WB face states, sp x3 shift)
which sits outside the solver branch in hydro_fluxes.cpp; no grid guard on rsolver.

**Compare (dhjcs.py pipeline, hydro_w):** zonal-mean jet u(lat, p) at matched rot,
KE3 history, deep v_r floor and its vertical odd-even content (the decoupling
signature), mass loss. Expect: lhllc/ausmpup less dissipative -> possibly stronger jet
and more deep residual; whether the JET differs is the question.
