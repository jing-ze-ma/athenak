---
name: cs-stretched-source-term-bug
description: SECOND stretched-grid cs bug, FIXED 13a97399 -- the angular-momentum curvature source used the index spacing as dr and an unstretched r, off by 0.57x-2.15x across the production grid; plus the full sp-vs-cs source audit (what else differs, what was checked clean) and the cs_test stretch fix
metadata:
  type: project
---

**Found 2026-09-06 by auditing for the [[cs-stretched-resistive-rcm-bug]] pattern** (a radial
position or spacing taken from CellCenterX/LeftEdgeX/mb_size.dx1, i.e. INDEX space, on a
stretched grid). `SrcTermsGnomonicEquiangleImpl` converted the radial-face flux of each
angular momentum with dr/(2r), dr = mb_size.dx1, r = CellCenterX -- neither stretched. On
the production dhj grid (nx1 128, c1..c4, width ratio 4.07) the true/used factor runs
0.57 to 2.15, outside 0.8-1.25 in 94 of 128 cells. sp takes (r_r-r_l)/(r_r+r_l) from the
stretched faces; cs now does the same (13a97399). Uniform grid: identical in every digit.
Affects EVERY stretched cs dhj run, hydro and MHD (cs_prod_hyd_rot included) -- so the
5-25 % cs-vs-sp hydro disagreement in [[cs-dhj-production-retry]] was measured with this
bug in.

**cs_test was itself unstretched** (all 35 radial positions): every stretched cs_test
number before 13a97399 compared the code against an IC/BC/reference on a DIFFERENT grid
(the 4a16f07b message, the 2.8e-2 curl residual, the 6e-5 "ideal mass drift"). Now
ApplyRStretch (moved to coordinates/grid_stretch.hpp) is applied everywhere. After both
fixes on the 8-cell stretched grid: iprob=9 L1(v) 3.4e-3 -> 1.1e-4, L1(B) 3.8e-4 -> 4.4e-5;
iprob=11 curl OPERATOR residual 2.8e-2 -> 2.5e-4 (= uniform).

## The audit (2026-09-06): what else was checked
CLEAN: sp dual mesh (x1v_/xx1f_ arrays); Coordinates::dx1 (stretched); resistive dt
(pcoord dx); mhd_ct / hydro+mhd_update / newdt .dx1 uses (Cartesian branch only);
current_density.hpp size.dx1 (Cartesian branch); BuildWBGeometry (areas/volumes, angular
only from CellCenterX); mhd_corner_e + derived_variables CellCenterX (GR/Cartesian paths);
deep_hot_jupiter_rt.cpp: all 13 radial positions ApplyRStretch'd, RT uses pcoord dx1,
gravity/rotation via areas/volumes; stellar_tide on cs is REFUSED at startup (not silent).
DESIGN DIFFERENCES sp vs cs (not bugs): cs evaluates g(r), omega^2 r at the stretch-mapped
midpoint, sp at the volume centroid -- 7e-5 in g on the production grid; sp src3 (cot
term) is in flux form, cs in state form; sp subtracts the WB static background pressure
pwb in the geometric source when `wellbalance_static` is on, cs NEVER does (the cs WB
path passes false by design) -- LATENT: `wellbalance_static=true` on cs would be
inconsistent; production has it off on both grids.
See [[cs-stretched-resistive-rcm-bug]], [[radial-grid-stretch]], [[validate-the-instrument]].

## MEASURED IMPACT on the hydro production run (2026-09-06, cs_prod_hyd_rot rerun on 979edada)
Against the archived aa6ddb9b run and the sp hydro control sp_dhj_hyd at rot 5-80: mass and
tot-E now match sp to 4 digits at every rotation (3.454/3.454e26 at rot 5, 6.355/6.355e38)
where the old cs run was 0.1 % low in both from rot 5 on (3.451e26, 6.349e38) -- the
budget offset was the source-term bug. KE1 (radial) is still 20-30 % below sp (6.6e32 vs
8.4e32 at rot 40, old 7.1e32) and horizontal KE within 5-20 % either way -- unchanged, i.e.
turbulent scatter plus the 32x32-per-panel vs 64x128 resolution difference, not a budget defect.
