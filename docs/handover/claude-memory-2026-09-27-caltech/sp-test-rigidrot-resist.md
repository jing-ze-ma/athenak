---
name: sp-test-rigidrot-resist
description: sp_test now has iprob=3 RIGID ROTATION (hydro or MHD, exact steady state) and iprob=11 FORCE-FREE RESISTIVE DECAY (exact time-dependent solution, div B = 0 by Stokes loops), with exact radial user BCs, L1/Linf errors split polar-rows/interior, div B gate, user history; inputs/tests/spherical_polar_{rigidrot,resist}.athinput. Rates measured 2026-09-06 -- the polar rows do NOT converge for a tangential field
metadata:
  type: project
---

**What exists (src/pgen/sp_test.cpp, build_sp_test, 2026-09-06).**
* iprob 8 (default): the old uniform-field-at-rest test, unchanged.
* iprob 3: v_phi = omega r sin(theta), p = p0 + d0 omega^2 R^2/2; with <mhd> a uniform
  b0 zhat is added (v x B curl-free, J = 0: still exactly steady). Twin of cs_test iprob 3.
* iprob 11: B = b0 (sin(alpha z), cos(alpha z), 0), curl B = alpha B, |B| uniform, J x B = 0;
  exact B(t) = B(0) e^{-eta alpha^2 t}, p(t) = p0 + (gamma-1) b0^2 (1 - e^{-2 eta alpha^2 t})/2,
  total energy exactly constant. Faces from Stokes loops of A = B/alpha with 4-pt
  Gauss-Legendre edge integrals -> div B ~1e-20. The face ON the pole is set the way
  bfield_bcs.cpp sets it (half the difference of the adjacent face at phi and phi+pi), so
  the t=0 pressure is exact to 1e-16 -- without that the polar row starts 4 % off.
* `SPTestRadialBC` (ix1/ox1 = user): exact state incl. ghost FACES, time-dependent for 11
  (`problem/bc_time_frac`). `SPTestErrors` (pgen_final_func): L1/Linf of v, p, rho, faces
  of B, polar rows (`problem/err_nband` cells from either pole) vs interior, location of the
  max B error, max |div B| L/b0; appended to <basename>-errs.dat. `SPTestHistory`
  (user_hist): L1-v, L1-p, L1-rho, L1-b -> <basename>.user.hst.
* TRAPS: the polar boundary needs an EVEN number of MeshBlocks around phi (meshblock/nx3 =
  nx3/2); a cmdline override needs the parameter to EXIST in the input; the explicit
  resistive dt at the pole scales as (dtheta dphi)^2 ~ h^4 (7.4e-4 at n16, 4.4e-5 at n32,
  2.7e-6 at n64: n64 infeasible); bin dumps go to ./bin/.

**Rates (L1, shell r=1..2, nx2 = 8/16/32(/64), tlim 1 for iprob 3, 0.5 for 11), after the
two fixes [[sp-pole-edge-area-zero]] and [[resistivity-3d-curl-missing-terms]]:**
  rigid rotation hydro: v 2.57, p 2.1 (interior 2.55); POLAR ROWS v 1.85 -> 1.51.
  force-free IDEAL (eta=1e-12): v 1.32 -> 1.76, p 1.8, B 1.9; interior 1.8; POLAR ROWS
    v 0.10 -> -0.04 (L1 1.7e-2 at EVERY resolution: an O(1) polar-row defect for a
    tangential field), B 0.6-0.8.
  force-free RESISTIVE: B 2.1 (interior 1.7), p 1.6, v 1.4; polar v 0.31; log-decay ratio
    1.35 (n8) -> 1.04 (n16).
  Existing pole flags `sp_cart_polar_momentum`, `polar_quadratic_recon` change NOTHING here
  (rigid rot polar 7.8e-4 -> 7.3e-4 at best); `use_polar_average_eresist` idem at n16.
The cs twins for comparison: cs rigid rotation n16 L1(v) 3.5e-4 vs sp 4.1e-4.

**div B and the pole EMFs (iprob 11, n16, tlim 2, fixed binary).** div B is 4e-14 at t=0 but
5.7e-4 (in b0/L) at t=2 with the defaults: the resistive E_r at the pole edge is added AFTER
CornerE's azimuthal average, so the k-dependent e1 on the one physical pole edge breaks the
polar-row divergence. `mesh/use_polar_average_eresist = true` restores 6e-13 with the errors
unchanged -> it should be the default (or unconditional). `mhd/polar_emf_diss` (default ON)
COSTS accuracy on this smooth problem: L1(B) 4.2e-3 -> 1.8e-3 and polar-row L1(B) 1.7e-2 ->
2.3e-3, polar v 2.4e-2 -> 1.6e-2 with it OFF. Its stabilising role is [[sp-polar-field-blowup]].

## 2026-09-06 02:00: sp_test now has iprob 12 (TOROIDAL, closed conservation gate, 8addac85)
and iprob 13 (BLAST on the sphere, pole vs equator, 96a18028); L_z (exact lever arm) in
every user history. Toroidal with REFLECTING walls: mass, energy, L_z to round-off at
nx2 8/16/32 (the closed MHD gate PASSES); with resistivity the conducting wall admits only
40 % of the needed Poynting influx (ratio 0.38/0.41), so the Ohmic budget is read from the
open-boundary run: ratio 0.90 (n8) -> 0.98 (n16). Rates: global v/p/B 2.0-2.3, polar v
1.1 -> 2.0, polar B 0.35 -> 0.62 (an axisymmetric poleward drift of B_phi in the polar cell,
first-order hoop-stress balance). Blast n16: pole vs equator agree to 0.5-2 % L1 of the
profile range; the centre bin is a crude comparison at that resolution.
The MHD blast/toroidal/FF radial user BC holds the ambient (open) state; `reflect` with
B_r != 0 at the wall GAINS mass 1 %/0.5 tu with HLLD (HLLE 100x less) and energy 1.6 %
(the wall-face EMF dissipates the 2 B_t jump): reflect is the wrong wall for a threading
field, not a conservation bug.
