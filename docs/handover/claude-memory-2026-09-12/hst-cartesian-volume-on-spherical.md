---
name: hst-cartesian-volume-on-spherical
description: the .hst history output weights every cell by dx1*dx2*dx3, so on the spherical dhj mesh "mass"/"tot-E" are integral dr dtheta dphi, NOT physical volume integrals
metadata:
  type: project
---

`src/outputs/history.cpp` computes every history variable as
`vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3` (lines 123, 209, 309) --
the **Cartesian** cell volume. On the `deep_hot_jupiter_rt` spherical mesh
(x1 = r, x2 = theta in [0,pi], x3 = phi in [0,2pi]) that is `dr dtheta dphi` with no
`r^2 sin(theta)` factor.

**Why:** AthenaK's history output was written for Cartesian meshes and was never given a
coordinate volume; the spherical/cubed-sphere work on this branch did not update it.

**How to apply:** never quote `.hst` mass/energy as a conservation check on this pgen.
Verified 2026-08-24: recomputing `sum(rho * dx1*dx2*dx3)` from the `.bin` dumps reproduces
the hst `mass` column to 7e-8 (float32 noise), while the physical
`sum(rho * (r2^3-r1^3)/3 * (cos th1 - cos th2) * dphi)` gives a drift **8x larger**
(+2.10e-3 vs +2.66e-4 over t = 4.32e5). The two disagree because mass redistributes
strongly in radius, and the unweighted integral scores the tenuous outer shells the same
as the dense bottom ones. Do the integral from `.bin` instead -- recipe and numbers in
[[dhj-conservation-check]].
