---
name: dhj-isobar-vs-shell
description: The 1e-6 bar level must be measured on the ISOBAR, not a constant-r shell; the floors never touch it, but the dayside isobar is outside the domain
metadata:
  type: project
---

> **NUMBERS REVISED:** the first pass used a buggy EOS inversion
> ([[eos-inversion-nan-trap]]). The CONCLUSION held and strengthened -- 68.3 % of terminator
> columns had the isobar outside the retired domain, not 43 %. Pressures/temperatures below
> from that first pass are superseded.

Measured 2026-08-25 (eighth session) on the [[ck-hydro-long-run]] evolved state via
[[eos-table-dump]]. **Supersedes the floor conclusions in [[dhj-floors-for-1e-6-bar]] and
my own earlier numbers in this campaign.**

**The mistake to not repeat: a constant-radius shell is NOT the 1e-6 bar level.** At
r = 1.267e10 the horizontal MEAN pressure is 1e-6 bar, but the dayside there is 3.8e-6 bar
and the nightside 3.8e-9 bar -- a three-decade contrast the mean hides completely. Measuring
"floored cells at the 1e-6 bar level" on that shell gave 21 % (42 % nightside). **On the
actual p = 1e-6 bar isobar it is 0.0 %**, in every arm including the dfloor = 5e-13
baseline. Density on the isobar is 2.3e-11 median, 5.0e-12 minimum -- a **10x margin above
dfloor**, and consistent with the ideal-gas estimate p*mu*m_H/(kT) = 2.8e-11, which
independently validates the EOS inversion.

**So the floors do not contaminate the observable at all**, and relaxing `pfloor` buys
nothing for the science: it would only let the nightside near-vacuum ABOVE the observing
level fall further, a region that is not observed and is numerically fragile.

**The real blocker is x1max, on the DAYSIDE.** The outermost cell holds 3.8e-6 bar (median,
up to 7.7e-6) on the dayside against 3.8e-9 bar on the nightside, so **86.7 % of dayside
columns never reach 1e-6 bar inside the domain** (43.6 % of the planet; 0.6 % of the night
side). Extrapolating with the local scale height (2.96e8 cm at the dayside top) the domain
needs **+4.4e8 cm median, +5.9e8 at the 90th percentile -> x1max ~ 1.36e10** to capture the
isobar everywhere. That is a 4.4 % extension, much less than the 1.42e10 that
[[dhj-floors-for-1e-6-bar]] measured as harmful -- and that harmful measurement was
nightside-dominated (a stagnant floor-supported halo), which a lower dfloor now mitigates.
The two conclusions are about different hemispheres and are not in conflict.

**x1max is a grid change**, so it cannot be restarted into; testing it needs from-scratch
runs of several rotations, ideally as a clean 2x2 with dfloor.
