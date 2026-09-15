---
name: lhlld-solver
description: 2026-09-13 LHLLD ALIGNED with Minoshima & Miyoshi 2021 (branch lhlld, wt_lhlld, HEAD ac578d99, NOT pushed): momentum weight ONE (base star flux already = F(U*); energy keeps the Batten weight S_M S_a/(S_a-S_M)), chi = min(1, max(c_u)/max(c_f)) with c_u the fast-speed formula with |u| in place of c_s -> MAGNETIC FLOOR; the Mach-1e-6 PPM Alfven CHECKERBOARD IS GONE, lhlld tracks hlld to 0.1 % in the whole lwave suite (skips removed). Bitwise = hlld at Mach 5.7. Literature tests (tst/lhlld_bench/): CPAW 2nd order both; Orszag-Tang null pass; Gresho beta 100 gain only 1.03-1.15 (vortex not an MHD equilibrium, field winds up) -> BALSARA VORTEX pgen DONE 09-14 (balsara_vortex, inputs/mhd/balsara_vortex.athinput, test_mhd_balsara_vortex_cpu.py): at M_Alf 10 lhlld error is MACH-INDEPENDENT (0.010-0.012 vs hlld 0.02->0.13, Ek retained 0.87 vs 0.30, hlld needs 1.2x (M 0.1) to 3x (M 0.01) resolution) = Leidi+22 reproduced; advantage vanishes at beta_K>=1 (magnetic floor, by design). OPEN: BLOW-UPS (no NaN, energy x100-1000) at vortex Mach 1e-3 with beta_K (mag/kin) 1e2 for BOTH solvers and beta_K 1 for lhlld; onset 0.2-0.4 crossings; also rk3/ppmx at M 1e-2 beta_K 1 grows -> an instability of low-Mach strong-field flow, undiagnosed, matters for an MHD FeCZ box near equipartition; KH growth inconclusive (rk3+wenoz reconstruction-dominated; redo plm 16/32^2). theta shock factor not implementable (needs transverse stencils). No GPU/MPI run.
metadata:
  type: project
---
Tests: tst/test_suite/mhd/test_mhd_lhlld_lowmach_cpu.py (beta 100, 64^2, gates hlld>=0.20, lhlld>=0.30, gain>=1.15).
chi = min(1, max|u_n|/max c_f) uses the interface-normal velocity, matching lhllc_hyd.hpp; the paper may use |u|.
See [[fecz-box-projects]].
