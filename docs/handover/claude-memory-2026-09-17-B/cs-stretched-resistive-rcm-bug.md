---
name: cs-stretched-resistive-rcm-bug
description: cs RESISTIVITY on a STRETCHED radial grid was anti-diffusive (FIXED c5c85e3b) -- the dual-mesh centre-to-centre length mixed a stretched r_c with an UNSTRETCHED r_cm, changing sign near the top; exactly the production dhj configuration that builds the deep toroidal sheet
metadata:
  type: project
---

**Found 2026-09-06.** cs_test iprob=11 (constant eta 0.5, 8 stretched radial cells,
`mesh/use_grid_stretch_r_poly=true`) blew up by t = 0.011 of 0.2 on BOTH the pre- and
post-4a16f07b binaries: mass 29.3 -> 21.8, tot-E 43 -> 1.6e5, then dt = 0 grinding.
Discriminators (all on the stretched grid): eta = 0 CLEAN; eta 0.05 blows up later
(t ~ 0.1) and STILL blows up at cfl 0.1 -> not a timestep problem; uniform grid clean at
eta 0.5. So: a consistency defect in the resistive operator, present only with the stretch.

**Cause** (`src/coordinates/coordinates.cpp`, the dual mesh consumed by
`AddEMFGnomonicResist` via `dxface.x1f` and `areaedge.x2e/x3e`): `r_c` was put through
StretchR/StretchRPoly, `r_cm = CellCenterX(i-1)` was NOT (since 4bfacdd8, the original cs
resistivity). On the production stretch (max/min width 2.45 over 8 cells) r_c - r_cm was
wrong everywhere and NEGATIVE near the top -> J with the wrong sign -> anti-diffusion.
**Fix c5c85e3b:** r_cm gets the same stretch. Stretched iprob=11 runs to tlim, Ohmic
ratio 0.978; uniform grid bitwise identical.

**Why it matters:** the production dhj cs MHD runs are stretched + eos resistivity, and the
from-scratch ablation ([[cs-deep-toroidal-sheet]]) found the eos resistivity BUILDS the
deep toroidal sheet while sp with the same eta has none. This is the cs-specific defect in
the resistive path that [[cs-deep-toroidal-sheet]] predicted. Every cs MHD production
result with resistivity before c5c85e3b is suspect (cs_prod_mhd_rot, cs_both_chain, the
ablation arms). NOT YET SHOWN: that the fixed binary does not build the sheet -- needs a
HIP rebuild at HEAD and a from-scratch cs_prod_mhd_rot rerun.

**Residuals still open on the stretched test:** the curl OPERATOR residual is 2.8e-2 on 8
stretched cells vs 2.5e-4 uniform (the centre-to-centre Stokes loop is off-centred on a
non-uniform grid); and the ideal (eta 0) stretched run drifts 6e-5 in mass where the
uniform one is exact (a stretched radial-BC/pgen effect, not resistive).
See [[cubed-sphere-resistivity]], [[resistivity-audit]], [[validate-the-instrument]]
(the resistive test was gated on the UNIFORM grid only).
