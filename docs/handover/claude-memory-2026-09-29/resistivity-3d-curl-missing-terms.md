---
name: resistivity-3d-curl-missing-terms
description: 3D RESISTIVITY BUG (Cartesian + spherical polar, NOT cs) since 1123736f 2026-06-20 -- CurrentDensity gated the x2-derivative terms of J1 and J3 on the EXCLUSIVE Mesh::two_d flag, so every 3D run dropped d(B3)/dtheta and d(B1)/dtheta from curl B. FIXED 2026-09-06 (multi_d) in current_density.hpp; found by sp_test iprob=11. Every sp resistive dhj run since June is wrong
metadata:
  type: project
---

**The bug (src/diffusion/current_density.hpp).** The original AthenaK CurrentDensity gated
the x2 terms on `b.x1f.extent_int(2) > 1` ("at least 2D"). Commit 1123736f (2026-06-20,
"Generalize the Ohmic diffusion module") replaced that with `pmesh->two_d`, and
d52f3610 copied it into `geom.two_d`. Mesh::two_d is EXCLUSIVE (mesh.cpp: `if nx3>1
three_d; else if nx2>1 two_d`), so in 3D the block
`j1 += d(dx3*B3)/dx2 ; j3 -= d(dx1*B1)/dx2` never ran: J_r lost its theta-derivative and
J_phi lost its theta-derivative, on BOTH the Cartesian and the spherical-polar branch. The
cubed sphere is unaffected (its own two-pass curl in resistivity_gnomonic.cpp, gated by
cs_test iprob=11 at 2nd order -- which is why the cs test never caught it).

**How it showed.** sp_test iprob=11 (force-free B = b0(sin az, cos az, 0), exact decay
exp(-eta a^2 t)): the field decayed at 0.47x the exact rate at nx2=16 and SLOWER with
resolution (0.75x at n8, 0.41x at n16); L1(B) did not converge (rate 0.09, interior -0.32).
With `geom.multi_d` the log-decay ratio is 1.35 (n8) -> 1.04 (n16) and L1(B) converges at
2.1 (interior 1.7). No 3D resistive test existed anywhere (inputs/mhd/resistivity.athinput
is 2D, tst/ has none) -- see [[validate-the-instrument]].

**Consequences.** Every 3D spherical-polar resistive run since 2026-06-20 used the wrong
current: sp_mhd_diss (11420831), sp_mhd_pole2 (11428062), the sp resistive dhj runs,
ck_mhd_b3 smoke, the xe runs ([[xe-resistivity-long-runs]]), and the "resistive dt / Ohmic
heating" numbers in [[resistivity-audit]] (which certified this kernel CLEAN -- wrongly; the
audit read the geometry, not the dimension guard). The eos-eta production path goes
through the same CurrentDensity. cs production (cs_prod_mhd_rot, rcmfix) is NOT affected.
The sp-vs-cs MHD comparisons of [[sp-hydro-vs-mhd-comparison]] compare a wrong sp
resistivity with a right cs one. Status at write time: fix built in build_sp_test only,
NOT committed; the HIP dhj binary and the running sp jobs still carry the bug.
