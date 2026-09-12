---
name: stellar-tide-flag
description: problem/stellar_tide adds the host star's Omega^2(2x,-y,-z) tidal term; implemented and validated, UNCOMMITTED, deliberately left OFF
metadata:
  type: project
---

Written 2026-08-26 (ninth session). **Committed as c4aa730f** on `polar-average-perf`, on top of c37ebe75. Not yet pushed.
Supersedes the "missing" note in [[stellar-tide-at-domain-top]].

**What.** `<problem>/stellar_tide` (default false). The corotating frame already applies
Coriolis and the planet's own spin centrifugal `oor = Omega^2 r sin(theta)`; the full Hill
acceleration is `Omega^2 (3x, 0, -z)`, so the missing residual is exactly
`a_tide = Omega^2 (2x, -y, -z)` -- verified to 3.6e-15 numerically. Spherical components in
`TideAccR/TideAccT/TideAccP` (grid_stretch-style file-scope helpers in the pgen):
`a_r = Om^2 r (3mu^2-1)`, `a_th = 3 Om^2 r sinth costh cos^2ph`,
`a_ph = -3 Om^2 r sinth sinph cosph`, with `mu = sin(theta)cos(phi)` -- the SAME `mu0` the RT
uses for the stellar beam (`phi = x3v - pi`, substellar at lon 0).

Applied at 4 places: 3 momentum components + the work term in the corotating block, and the
RT `tau = kappa p/g` closure via a new `EffGravAt` (3 call sites, clamped at 0.1*g).
Deliberately NOT in the 1D hydrostatic IC integrators (no angle; background must stay
spherical) nor the MHD-only Maxwell ghost clamp. Not in `phicc0` either -- it is a source
term like `oor`, which is fine because **well-balancing is OFF in this campaign**
(`wellbalance_static/dynamic`, `wb_x1..x3` all 0 in the restart file), so gravity itself is
a plain rho*g source.

**Validated (jobs 11059855/6/7, 0.25 rot, production grid):**
* flag OFF is **bitwise identical** to c37ebe75 -- 5 vars x 32 meshblocks x 6 dumps
* ON: terminator isobar **-4.14 cells** (hydrostatic prediction -3.24), substellar **+8.11**
  (predicted +9.30). Right sign and magnitude both ways.
* deep-interior rms v_r decays in both arms but the ON arm sits **4x higher** and had not
  settled at 0.25 rot -- the startup transient from building the IC without the tide.

**DECIDED: leave it OFF for the campaign.** The tidal forcing is EXACTLY mirror-symmetric
about the star-planet axis (a_r, a_th identical at both limbs; a_ph equal and opposite), so
it contributes **zero leading/trailing asymmetry at first order** -- and limb asymmetry is
the observable. Its direct effect is a ~1% uniform radius shift, largely degenerate with the
retrieval reference radius, and ~3% the size of the constant-g error [[grav-point-mass-flag]]
fixed. It matters if the science ever moves to the outer atmosphere: at the domain top
`3 Om^2 r/g` = 13% and the top is at 51% of r_L1 = 4.28 Rp.
