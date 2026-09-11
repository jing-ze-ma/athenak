---
name: sp-pole-edge-current-sign
description: sp RESISTIVE POLE BUG #2, FOUND 2026-09-06 late -- the Stokes loop for J_r on the POLE EDGE ran its far-side phi-segment through the polar ghost row, whose B_phi is sign-FLIPPED by the exchange, so the loop never closed for a curl-free field and J_r ~ 3 B sin(phi)/(r dtheta) GREW with resolution (polar E1 error 0.9 -> 2.1 eta*alpha*b0 from nx2 8 -> 16); its Ohmic heating is the m=2 polar pressure pattern. Fixed (sign per pole + 3-face parabola through the pole), rate 1.7
metadata:
  type: project
---

**Instrument:** sp_test iprob=11 finalizer now compares the code's resistive EMF eta*J on
every x1/x2/x3 edge with the exact eta*alpha*B (curl B = alpha B for the force-free field),
in units eta*alpha*b0, interior vs polar rows, north/south maxima. Before the fix: E1
(r-edges) interior 4.3e-2 -> 1.3e-2 (2nd order) but polar L1 0.93 -> 2.08, Linf 3.9 -> 8.8,
max ON the pole edge (j-js = 0). E2/E3 fine (first order in the polar rows, converging).

**Mechanism (current_density.hpp, sp branch, 3D term of J1):** the dual loop of the pole
edge passes through the pole to the far side; the code's far-side phi-segment is
-dx3(j-1) B3(j-1) with B3 from the polar GHOST row, and the exchange stores the ghost's
B_phi with its sign flipped (phihat reverses through the axis). For B = B xhat the true
phi-segments cancel the theta-segments; with the wrong sign they add: circulation
(4/3) r dtheta dphi B sin(phi) over the cap area r^2 (4/9) dtheta^2 dphi = 3 B sin(phi)/
(r dtheta). Predicted 6.4 eta*alpha*b0 at nx2 = 16, measured Linf 8.8. Its Ohmic heating
~ sin^2(phi) is the m = 2 polar-row PRESSURE pattern that GREW with resolution (2.5e-3 at
n16, 4.7e-3 at n32, m0 part converging) and drove the non-converging polar-row v_theta,
v_phi of the resistive test (flat at 7e-3 / 1e-2 from n16 to n32 while ideal MHD halved).
The resistive E_r pole AVERAGE hides it in the CT (sin(phi) averages out) but NOT in the
energy flux, which uses the unaveraged efld_resist. Same class as the c5c85e3b/3694cc40
dual-mesh bugs: the polar row's mirrored ghost geometry is not an ordinary row.

**Fix:** (1) flip the sign of the ghost-row segment, per pole (north: the j-1 term;
south: the j term -- flipping the same term at both poles doubled the south error);
(2) the through-pole theta-segment was a one-point trapezoid (first order): the pole face
is by construction the MEAN of its neighbours so a Simpson with them is a no-op; the
parabola through the faces at -dtheta (ghost), +dtheta, +2 dtheta gives
Int = 2a [B0 + c2 (a^2/3 - dtheta^2)], c2 = [B(2dth) - 1.5 B(dth) + 0.5 B(-dth)]/(3 dth^2).
Result: polar E1 L1 0.195 -> 0.094 (sign only, rate 1.0) and 0.091 -> 0.029 with the
parabola (rate 1.7); Linf 0.47 -> 0.15. See [[sp-pole-fixes]], [[sp-pole-edge-area-zero]].
**nx2 = 32 (committed as the two commits after 0f9959ee):** resistive test polar-row L1(v)
6.0e-3 -> 3.65e-3 (rate 0.43 -> 1.11), Linf 2.0e-2 -> 6.5e-3 (-0.1 -> 1.2), polar L1(B)
4.4e-4 -> 2.9e-4 (2.7); global unchanged (v 1.9, p 1.8, B 2.5). Polar E1 L1 at n32 0.019
(max now at the OUTER radial boundary, south pole: the next thing to look at if pursued).
The sign fix alone gives the same velocity/pressure numbers; the parabola improves E1 only.

## CORRECTION (2026-09-06 01:30): the figure-eight was wrong too -- the POLE-EDGE CURRENT
IS ONE SECTOR OF THE CAP (581a1e43). The b61c5d67 loop (through the pole and back) encloses
zero net area: J_r = 0 for a uniform axial current (the new toroidal test iprob=12,
8addac85, showed E1 error exactly 1.0 on the pole edge). The dual face of the axis edge is
one cap sector bounded by the polar cell's phi-FACE at its own latitude dtheta/2 (NOT the
centroid 2dtheta/3: that gave 3/4 of the current -- the midpoint-vs-centroid trap again),
half-meridians from the pole with the parabola mean of B_theta, area r^2(1-cos a)dphi; the
polar E_r average then yields the exact cap current. Toroidal field: pole-edge eta J to
0.3 %, Ohmic budget ratio 0.90 -> 0.98 (n8 -> n16). FF field unchanged (J_r = 0 there).
Operator errors at t=0 (E1, units eta|J|): FF polar 3.8e-2 -> 2.3e-2 (rate 0.7, max ON the
pole edge), interior 1.3e-2 -> 3.7e-3; toroidal polar 4e-3 -> 6e-3 (max at the outer wall
corner), interior 2nd order. The t=0.5 E1 numbers are dominated by the evolved field.
iprob 12 is also the CLOSED conservation gate (B_r = 0 at the walls) and L_z (exact lever
arm) is now in every sp_test user history.
