---
name: sp-pole-fixes
description: THE sp POLE, 2026-09-06 -- the non-converging polar-row error of a tangential field had TWO causes, both FIXED and committed (a8cfb83e face-basis rotation of x3-face states; de667f32 third-difference polar_emf_diss) plus b84a0502 (eresist pole average default ON): polar-row L1(B) rate 0.76 -> 2.1, v_r 0 -> 2.3. The checkerboard gate on the production dhj run is PENDING (bench/sp_diss3)
metadata:
  type: project
---

**Symptom (sp_test iprob=11, ideal control eta=1e-12, tlim 0.5):** the polar-row RADIAL
velocity grew linearly at the SAME rate at nx2 = 16 and 32 (2.5e-2 at t = 0.5), i.e. a
resolution-independent spurious radial force ~0.6 (B^2/2)/r in the single polar row,
while v_theta, v_phi and the interior converged. The one-cycle force is only 0.14 in those
units and m = 2 in phi; it grows to 0.6 within ~0.05 time units as the polar-row FIELD is
distorted (m = 1 amplitude of bcc_r in row 0 rose 0.065 -> 0.113 by t = 0.5).

**Cause 1 -- `polar_emf_diss` (default ON since 3147de4f).** Its Rusanov term on e2 used
FIRST differences of the face fields along phi. In the polar row dx_phi = r sin(theta)
dphi ~ dtheta dphi while a smooth m = 1 pattern has kappa ~ 1/(r sin theta), so the
damping rate a dphi/(2 r sin theta) is O(1) at every resolution: the stabiliser was
eating the smooth field. With it OFF the polar v_r error dropped 10x. Fix (de667f32):
THIRD differences, -D3/4, identical to D1 on the (-1)^k checkerboard it exists to damp
(by construction) and O(kappa^3 dx^3) on smooth modes.

**Cause 2 -- x3-face reconstruction location (a8cfb83e).** The x3 sweep reconstructs
along phi at the cell's VOLUME centroid theta = x2v, but the x3 face (area r dr dtheta)
is centred at the theta MIDPOINT: polar row 2 dtheta/3 vs dtheta/2. The r-components of
v and B scale as sin(theta) there, so the flux T_phi,r through the polar cell's
phi-faces was 4/3 too large (a discrete balance of the polar cell shows the six face
terms cancel only with the midpoint value). Rotating (v_r, v_th) and (B_r, B_th) of the
reconstructed states into the face's own basis (rhat' = rhat cos d + thhat sin d, hydro
and MHD sweeps) took the diss-off polar v_r error from 6.7e-3 -> 3.4e-3 (rate 1.0) to
2.6e-3 -> 5.2e-4 (rate 2.3). Bit-identical on other grids and on axisymmetric states.
Note the earlier `sp_cart_polar_momentum` (34882bd7) does NOTHING on this evolved test
(one-cycle force 0.14 -> 0.39, worse) although it removed the one-step x-field residual
90x; it stays default OFF.

**Result with both, default flags (FF ideal, nx2 16 -> 32):** polar-row L1(B) 1.7e-2 ->
2.3e-3 and 1.0e-2 -> 5.3e-4 (rate 0.76 -> 2.1); polar v_r 2.5e-2/2.6e-2 -> 2.5e-3/5.2e-4;
global L1(v) rate 1.9, L1(B) 2.05. Rigid rotation bit-identical. Resistive n8->16 B 2.2.
The residual polar-row v_theta/v_phi rate ~1.4 is a fixed-2-cell band, not an order.

**FINAL RATES on b84a0502 (defaults), nx2 8/16/32/64:** FF ideal global v 1.4/1.9/2.1,
p 1.5/1.8/2.0, B 1.9/2.1/2.0; polar-row B 2.1/2.1/0.7, v 0.9/1.3/1.9. FF resistive (8/16/32)
v 1.4/1.9, p 1.6/1.8, B 2.2/2.4; polar v 0.8/0.5 (a fixed 2-cell band, not an order). Rigid
rotation v 2.6/2.6/2.4. Other tests: cs bit-identical; sp_lowbeta monopole identical to 5
digits; sp_lowbeta tangential STILL grows (KE 1.7e-3 vs 1.8e-3 -- not a pole effect); pole
blast 4th digit. One-step sp_test iprob=8 transverse-field polar-row force: r 1.86 -> 0.03,
but PHI stays 0.27 (32x64) -> 0.19 (64x128), rate 0.5 -- the next polar residual to chase.

**GATE SUBMITTED: bench/sp_diss3 job 11474015 (2026-09-06 20:15, HIP binary md5 8f239293
at b84a0502, 1 h apu1).** Pass = 1-ME ~7e31 through rot 1.1 (sp_diss); fail = the
sp_dhj_ctl jump to 9.7e33 at rot 0.8.
**PENDING GATE (old text):** the third-difference form has not yet been shown to stop the dhj polar
checkerboard blow-up ([[sp-polar-field-blowup]]: 56x jump in 1-ME at rot 0.7-0.8 without
diss; sp_diss with the Rusanov form was clean to rot 1.13). bench/sp_diss3 is staged
(sp_diss input, from scratch, 1 h apu1 slot); needs the HIP binary at >= de667f32. Gate:
1-ME at rot 0.8-1.1 tracks sp_diss (~1e32-ish, no jump). Also [[sp-test-rigidrot-resist]].

## Second pass (2026-09-06 evening): can the polar-row TANGENTIAL error be 2nd order?
Per-row, per-component (FF ideal, nx2 16->32): v_r 2.3 after the fixes, v_theta 1.6,
v_phi 1.4. Diagnosis by hand: the polar cell's phi-momentum is a cancellation of three
O(1/dtheta) terms (outer theta-face flux, phi-face flux difference, cot source) each with
O(h^2) face-centre quadrature -> O(h) residual; predicted 0.26 at 32x64 vs 0.27 measured.
* `sp_cart_polar_momentum` REPAIRED (eb8e3cb1: phi-face triads at the face midpoint,
  consistent with the flux basis) and RE-MEASURED: theta 1.7, phi 1.3 -- NO better, level
  slightly worse. It removes the cot source but keeps the same face quadrature. Stays OFF.
* FOURTH-ORDER FACE AVERAGES (F + delta2_t F/24 on all x2/x3 faces) TRIED and REJECTED:
  they break the exact point-value flux / geometric-source balance -- interior rigid
  rotation error x25, through-pole test rate 0.8. A consistent version needs the source
  rebuilt from the same face averages (not attempted).
* FULL THETA SHIFT of the x3-face states (863e8337, all variables, centred theta
  derivative of cell values; supersedes the rotation): fixes an O(1) polar-cell phi-force
  for a scalar gradient across the pole. New test sp_test iprob=3 `rot_axis = 0` (rigid
  rotation about xhat, flow THROUGH the poles, 8e55c4f0): polar-row L1(v) 9.1e-4 -> 3.3e-4
  at nx2=64 and rate 1.0 -> 1.5, Linf plateau (1e-2) removed; global v 2.3/2.3/2.2. FF
  ideal j0 with the shift: r 8.1e-4, th 2.1e-3, ph 4.1e-3 at nx2=32 (rotation-only: 5.2e-4,
  2.5e-3, 3.9e-3) -- a wash there, decisive on the through-pole flow.
VERDICT: the polar-row tangential error is ~1.5 order and stays so; 2nd order needs a
fully consistent higher-order quadrature of the polar cell (fluxes AND source). Not done.
**GATE UPDATE: sp_diss3 (11474015) at rot 0.80 has 1-ME 7.65e31 (control jumped to 9.7e33
there): PASSING; read through rot 1.1.**

**GATE READ (2026-09-06 21:45): sp_diss3 reached rot 0.92 in its 1 h slot: 1-ME 7.5e31 (0.70),
7.7e31 (0.80), 8.0e31 (0.92), tot-ME 9.2e32 -- tracking sp_diss (Rusanov form), NO trace of
the control's 56x jump at rot 0.7-0.8. The third-difference dissipation PASSES the
checkerboard gate through rot 0.92; a 1 h continuation was submitted to reach 1.1+.**
**GATE PASSED (2026-09-06 22:05): sp_diss3 continuation (11476345) reached rot 1.20 with
1-ME 7.3e31, tot-ME 1.0e33 -- the control had 1.45e34 at rot 1.1. The third-difference
polar_emf_diss holds the checkerboard on the production dhj configuration. CLOSED.**

## Third pass: whole-mesh Cartesian momentum (`sp_cart_all_momentum`, commit after dadf1095)
Cartesian update on EVERY sp cell (source masked everywhere), default OFF. Rates, nx2 8/16/32/64:
  rigid z: global v 2.6/2.5/2.3, polar v 2.3/1.8/1.9 (2.0e-5 at n64 vs default 9.2e-5)
  rigid x (through pole): global 2.3/2.3/2.3, polar 1.7/1.7/1.5 (2.1e-4 vs default 3.3e-4)
  FF ideal MHD: global v 1.5/2.0/2.1, B 1.9/2.0/2.0; polar v 0.9/1.3/1.7 (1.15e-3 at n64 vs
    default 9.5e-4: NO gain for MHD); per-row j0 rates n32->64: r 2.0, th 2.2, ph 1.7 --
    the low-res MHD rates were PRE-ASYMPTOTIC (field scale 1/alpha = 0.64 rad, 3 cells at n16).
  Face averages (sp_face_avg) on top: WORSE everywhere (rigid x polar 6.1e-4 -> 8.4e-4 at
  n32); first attempt read UNINITIALISED hydro fluxes (hydro sweeps skip the extra layer
  unless FOFC) -- fixed by extending the ranges; the honest result is still "worse". OFF.
VERDICT: no MHD-specific pole defect remains; tangential polar rows sit at 1.5-1.9 with
the Cartesian update, ~1.5 without. Whole-mesh Cartesian is a candidate default for sp
hydro (needs the dhj hydrostatic/production gate); not for MHD (no gain).
