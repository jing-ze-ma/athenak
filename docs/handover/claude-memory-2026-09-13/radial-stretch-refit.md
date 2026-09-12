---
name: radial-stretch-refit
description: The production f_stretch_r_c* coefficients are NOT equalising cells per scale height (3.7x spread); a refit against H_rho fixes it at EVERY nx1 and buys a 1.46x dt on the cubed sphere
metadata:
  type: project
---

Measured 2026-09-01, on the dhj MHD test at nx1 = 64/128/234 with the ck production
physics.  See [[radial-grid-stretch]] for the map itself and [[dhj-cubed-sphere-port]] for
the runs.

## The finding

The production coefficients (`0.080646, -3.645655, 5.697513, -5.725033`) leave a **3.7x
spread in cells per density scale height** at EVERY resolution -- 2.93 at the worst radius
against a median of 7.93 at nx1 = 128 -- and a dr ratio of 12.3 across the domain.  The
worst-resolved radius is the TOP of the domain (r/ap = 2.16), which is exactly where the
1e-6 bar level lives.

**The map is resolution-independent**: u(xi) fixes the DISTRIBUTION, so "refit for
nx1 = 128" is really "refit, period" -- the same coefficients are right at 64, 128 and 234.

## The refit that works, and the one that does not

Fit `u_target(xi) = r(xi*N_H)` in the basis xi^k (1-xi), k = 1..4, from a measured H(r).
N_H = 20.55 across this domain, so a perfect map gives nx1/20.55 cells per H everywhere.

  - **against H_rho = -(dln rho/dr)^-1: `-0.068392, -2.191487, 2.464818, -1.366698`.**
    Spread 3.7x -> **1.34x**, min cells/H 2.93 -> 5.50, dr ratio 12.3 -> 4.07.
    CONVERGED: refitting from its own IC moves the coefficients by < 0.009.
  - against the PRESSURE scale height (the objective the header documents): NO BETTER than
    the current one (KE 6.24e32 vs 5.94e32).  H_p/H_rho runs from 0.76 to 4.84 across this
    domain -- the top is far from isothermal -- and equalising H_p leaves the density
    structure unresolved where it matters.

## Measured, 40 cycles, KE and mass drift at matched physical time t ~ 400 s

```
  run                    KE          dmass/m0      dt
  cs 128 current      5.9363e+32    -8.530e-04    19.86
  cs 128 fit H_rho    2.9470e+32    -5.462e-04    28.97   <-- 2.0x KE, 1.6x mass, 1.46x dt
  cs 128 fit H_p      6.2433e+32    -9.224e-04    18.29
  sp 128 current      1.9211e+32    -6.009e-04    19.30
  sp 128 fit H_rho    1.0585e+32    -3.380e-04    19.25   <-- 1.8x KE, 1.8x mass, dt flat
  cs 234 current      5.6260e+31    -2.310e-04    10.86   (the reference)
```

**The dt gain is cubed-sphere only**, and it is a second-order effect of the gentler grid:
the smallest cell grows when the dr ratio falls 12.3 -> 4.1.  Spherical polar sees no dt
change because its dt is bound by the POLAR cells, not the radial ones.

At nx1 = 128 the refit gets the initial eint to 1.7074 (cs) / 1.7080 (sp) against the
nx1 = 234 reference 1.7297; the two grids agree with each other to four digits.  nx1 = 128
still does not reach nx1 = 234 quality (KE 5x higher, consistent with 2nd order).

## Applied (checked 2026-09-06: the refit coefficients ARE in every production input, sp and cs, at nx1 = 128; the note below is stale)

Applying this to the PRODUCTION nx1 = 234 input.  Predicted there: min cells/H 5.33 ->
10.02, dr ratio 12.45 -> 4.09.  It would change production answers, so ask first.
