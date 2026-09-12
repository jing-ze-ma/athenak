---
name: sp-pole-bisect-culprit-863e8337
description: sp polar MHD blow-up BISECTED to ONE commit, 863e8337 (x3-face theta SHIFT of every variable); round 3 read 2026-09-10; a8cfb83e (rotation only) and de667f32 are CLEAN; eb8e3cb1 bit-inert
metadata:
  type: project
---

**The sp polar-row radial-field blow-up (sp_mhd_prod dies at rot 10.3) is caused by commit
863e8337 "Spherical polar: shift the x3-face states to the face midpoint, all variables".**
Read 2026-09-10 from bench/sp_pole_bisect (gate: 1-ME, hst column 11, at rot 4.7;
P = 3.05e5 s from omega 2.06e-5; ideal_r input, ideal MHD, no resistivity).

| commit | arm | 1-ME rot 4.7 |
|---|---|---|
| 2a64e2c7 | curl2a64_r (ref) | 8.8e31 clean |
| 54af8380, a8cfb83e (rotation), de667f32 (polar_emf_diss) | r3_* | 7.8-8.0e31 CLEAN |
| eb8e3cb1 | r3_eb8e3cb1 | history BIT-IDENTICAL to de667f32 for 2.8 rot: INERT (sp_cart_polar_momentum default off); cancelled |
| 863e8337 | x3shift_r | 1.61e34 DIRTY, bit-identical to 18f5dd21 |

b84a0502 is resistive-only, inert on this input. No untested gap: interval = {863e8337}.

What 863e8337 does (src/mhd/mhd_fluxes.cpp ~l.786, hydro_fluxes.cpp ~l.658): after x3
reconstruction, add `fac*(w0(j+1)-w0(j-1))` to wl/wr for EVERY variable and to bl/br for
bcc(0..2), fac = (face-midpoint theta - x2v)/(x2v(j+1)-x2v(j-1)). UNLIMITED, applied at
all j, and reaches across the axis into the polar ghost. It replaced the exact (v_r,v_theta)
/(B_r,B_theta) basis ROTATION of a8cfb83e, which is clean. The x2-face states and the polar
corner EMF are NOT shifted, so the x3-face EMF and x2-face EMF at the polar edge became
mutually inconsistent at O(1) in the polar row -- top suspect for the growing B_r
checkerboard (same family as [[sp-polar-field-blowup]]).

**Why:** the change that improved the through-pole rigid-rotation test ([[sp-pole-fixes]],
polar-row L1(v) rate 1.0 -> 1.5) is the one that destabilises the production pole. A test
gain is not a production gate.
**How to apply:** production sp MHD must run with the ROTATION form. The fix in flight
(2026-09-10) is a flag `<mesh>/polar_x3_shift = rotate|shift|off`, default rotate, plus a
3-arm A/B on HEAD (shift / rotate / scalars-shift+vector-rotate) in bench/sp_pole_bisect/hd_*.
Every sp MHD dhj run between 09-05 and this fix (sp_mhd_prod, the 09-07 defaults) ran the
dirty shift. Supersedes [[sp-pole-bisect-round2-null]].

## A/B ON HEAD, RESULT 2026-09-10 ~13:30 (bench/sp_pole_bisect/hd_*_apu1, hd_scalar; binary c5740ec4 = 95f2b07c + flag)
All three reached tlim rot 5.0 cleanly. Gate 1-ME (hst col 11) at rot 4.70:
- shift  6.39e34  DIRTY (800x clean, and 3-KE 3x lower: a materially different solution)
- rotate 7.24e31  CLEAN
- scalar 5.35e31  CLEAN (scalars shifted, vectors rotated)
=> the VECTOR (v, B) x3-face shift is the destabiliser; the scalar shift is harmless. Merge condition met:
polar-x3-shift-flag merged into polar-average-perf, default `mesh/polar_x3_shift=rotate`; sp MHD production relaunched
from scratch on apu (see [[sp-mhd-prod3-run]]). The apudev 15-min chains agreed (shift 7.35e34, rotate 5.79e31 at rot 5).
