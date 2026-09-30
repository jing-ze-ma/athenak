---
name: sp-pole-edge-area-zero
description: sp POLE BUG, FIXED 2026-09-06 (coordinates.cpp) -- areaedge.x1e at the polar r-edge was |cos x2v(js) - cos x2v(js-1)| = EXACTLY 0 (mirror centroids), so the resistive J_r at the pole divided by zero and any resistive run with a tangential field at the pole went NaN in cycle 1 on a UNIFORM theta grid; production survived only because f_stretch_theta breaks the mirror symmetry
metadata:
  type: project
---

**The bug (Coordinates::CoordSphericalPolar, kernel "spcoord_add").** The dual face around
an r-edge spans theta from the centroid of cell j-1 to that of cell j:
`areaedge.x1e = r^2 |cos x2v(j) - cos x2v(j-1)| dphi`. At the pole the ghost cell's
centroid is the exact mirror of the active one (x2v(js-1) = -x2v(js), e.g. +-0.130815 at
dtheta = pi/16), cos is even, so the area is IDENTICALLY zero. CurrentDensity's sp branch
then computes J_r = circulation/area1 at j = js (and je+1) -> inf/NaN -> the whole north
polar row NaN after ONE cycle. Reproduced with sp_test iprob=3 (rigid rotation + uniform
axial field, J = 0 exactly) plus `ohmic_resistivity = constant, eta_ohm_const = 1e-12`.
Fix: at a polar edge the dual face is two cap slices, area = r^2 (2 - |cos a| - |cos b|) dphi.

**Why production never saw it.** The dhj sp inputs use `use_grid_stretch_theta = true,
f_stretch_theta = 3`: the stretched ghost centroid is no longer the mirror image, the area
is tiny but nonzero, J_r at the pole is huge but finite. On any UNIFORM-theta sp grid
(every test input) a resistive run with B_theta or B_phi at the pole dies immediately --
which is also why no sp resistive test ever existed. Cf. [[resistivity-3d-curl-missing-terms]]
(found in the same session, same test), [[resistivity-audit]] (said the pole-face curl was
clean: it checked the edge LENGTHS, not this AREA), [[sp-polar-field-blowup]].
