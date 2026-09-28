---
name: rt-rad-force-transverse-cs-bug
description: rt_rad_force transverse term (Prad grad w) divided by the COORDINATE spacing size.dx2/dx3, not the arc length -> 2e11x too large on the cubed sphere, r-times too large on spherical polar; fixed 4bcdc855 on he4-presn-global (DX2/DX3 helpers); box bitwise. dhj/rg sp runs with rt_rad_force carried it
metadata:
  type: project
---

Found 09-17 in the 1-D He4 gate (tests_3d/README.md): transverse KE amplified 1e10 per cycle on a
spherically symmetric state; driver = rt_rad_force. In src/utils/two_stream_rt.hpp the
`Prad grad w` term used `(rhoN(j+1)-rhoN(j-1))/(2 size.dx2)` where on the cubed sphere x2 is the
gnomonic angle in [-1,1] (dx2 = 0.5 at nx2 = 4) instead of the arc length r dxi ~ 1e11 cm. The
radial component used X1V differences (a length) and was correct, so the plane-parallel box never
saw it. Fixed 4bcdc855: DX2/DX3 = pcoord arc lengths on curvilinear meshes, size.dx on Cartesian
(box bitwise). u0(IM2/IM3) are covariant on cs, so no metric inverse is needed.

**Why it matters beyond the He star:** spherical polar has dx2 = dtheta (short by r): every
red_giant / deep_hot_jupiter run on sp with rt_rad_force = true carried a transverse radiative
force too large by r. Check whether the dhj/rg productions used rt_rad_force before trusting
their thin-layer horizontal dynamics. Not yet ported to rt-integration/polar-average-perf.
See [[he4-presn-global-plan]].
