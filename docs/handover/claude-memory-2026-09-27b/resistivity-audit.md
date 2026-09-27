---
name: resistivity-audit
description: 2026-09-04 audit of src/diffusion for bugs that could jeopardise sp or cs resistive runs -- ONE live gap found (polar EMF average excludes the resistive EMF, fixed behind <mesh>/use_polar_average_eresist); eta_b ghosts, resistive dt, and the pole-face curl are all clean
metadata:
  type: project
---

Asked: "is there any other bug in resistivity that might jeopardise sp or cs runs?" Checked,
each with the file:line, so the next reader need not redo it.

**CLEAN**
* `eta_b` (per-cell EOS/Perna resistivity) is filled over the FULL ghost-inclusive range:
  `SetResistivity(..., 0, n1m1, 0, n2m1, 0, n3m1)` from both `mhd_tasks.cpp:659` and the RKG
  path `resistivity_rkg_tasks.cpp:172`. The face EMF reads `eta_b` at i-1/j-1/k-1, so this
  matters, and it is right. On cs the seam ghosts' eta comes from the halo-filled w0, hence
  consistent with the ghost state.
* Resistive timestep (`Resistivity::NewTimeStep`, ~line 960-1010): per-cell `eta_b` with the
  PHYSICAL `pcoord->dx1/2/3` on both curvilinear grids. The index-spacing variant
  `NewTimeStepConstantGridResist` (mb_size.dx*, which are angles on [-1,1] on cs) exists but
  its only call, line 950, is commented out -- DEAD.
* Resistivity builds NO cell-centred field of its own; Poynting flux and Ohmic heating
  consume `bc` = bcc0 from RaiseVelMHD, so the stretched-grid interpolation (ff99c522)
  flows through.
* The curl at the POLE FACE: theta starts at EXACTLY 0 on the production mesh and the
  geometry forces `sinl = 0` there, so `area_.x2f` at the pole face is 0 -- but
  `CurrentDensity` does NOT use it. It uses `areaedge.x2e` and `dxface.x3f`
  (`current_density.hpp:40-60`), both built from the CELL-CENTRED theta
  (`coordinates.cpp` ~1466-1471: `... * fabs(sin(x2v_(m,j))) * ...`), so at the pole face
  they carry sin(dtheta/2): finite, and the same factor in numerator and denominator.
  The 3D kernel evaluates j = js..je+1 including that face, and the value is finite and
  of ordinary size. NOT a 0/0 and NOT a floored-huge eta*J. The polar B blow-up is not
  sourced by the curl geometry.
* The RKG STS sub-loop applies the same `ApplyPhysicalBCs -> BFieldBCs`, so the polar B
  average (`use_polar_average_b`) would cover it; production runs explicit anyway
  (`use_rkg_sts = false`).

**THE ONE LIVE GAP (fixed, flagged)**: `MHD::EField` runs `CornerE` (which ends with the
polar average of the ideal E_r) BEFORE `AddResistiveEMFs` into the same array, so the
resistive E_r at the axis was never averaged. `<mesh>/use_polar_average_eresist` re-applies
it after the resistive add. Default OFF; see [[sp-polar-field-blowup]].

**Traps met while gating** (see [[validate-the-instrument]]): the polar blast with constant
eta collapses dt to 4e-9 at cycle 1 (dx^2/eta with dx -> 0 at the pole) and grinds to NaN,
so any gate on it compares NaN with NaN; a from-scratch dhj gate cannot show the polar
flags acting because the seed field is axisymmetric at the axis; and restarting
sp_dhj_ctl's rot-10 file needs BOTH a fresh `-t` (it carries the old wall limit and exits
after one cycle "on wall clock limit") AND `nlim` above its true cycle 630647.
