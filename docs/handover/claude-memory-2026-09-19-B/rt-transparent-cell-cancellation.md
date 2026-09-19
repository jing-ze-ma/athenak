---
name: rt-transparent-cell-cancellation
description: "The RT source in an optically thin cell is CATASTROPHIC CANCELLATION: Ft and Fb agree to 1.7e-11 relative and their difference is round-off, which is the same size as the cell's whole internal energy. Fixed by a finite opacity floor, not by the semi-implicit form."
metadata:
  type: project
---

Found 2026-09-09 while adding an ambient medium above the red giant. Measured with
`problem/rt_apply_debug` (now wired into red_giant.cpp), one column, background cells:

```
i=300  T=400  e=6.3301e-12  Fb=5.671403e+09  Ft=5.671403e+09
       divF=5.7825e-12  Qs=0  Em=1.1890e-15  de/e=5.0000e-01   <- clipped
```

`Fb` and `Ft` are identical to every printed digit -- correct, transparent gas transmits
the flux unchanged. But `src = -(Ft-Fb)/dx` turns that into noise: `Ft-Fb ~ -0.098` out of
`5.67e9` is **1.7e-11 relative**, and over `dx = 1.7e10` it gives `divF = 5.78e-12`, the
same order as the cell's ENTIRE internal energy `e = 6.33e-12`. So the radiative "source"
in a transparent cell is floating-point round-off worth 50 % of its energy per step, and
`rt_de_max` clips every one of them.

## What did NOT fix it

- **The semi-implicit linearization.** Replacing it with a form that relaxes to the true
  equilibrium (see below) left the clip count **bit-identical at 18816**. The old
  asymptote `e(A-E)/4E` really does diverge as `E -> 0`, so the fix is right on its own
  terms, but it is not what drives this.
- **Zeroing the opacity below the table window** (log10 kappa = -30). That makes it WORSE:
  `Em ~ kappa` drops to 1.19e-15, three orders BELOW the noise, so both the old and the
  new update degenerate to plain explicit and get clipped identically.

## What does

A **finite** opacity floor, `problem/opac_floor = 1e-5 cm^2/g`, for the nodes below the
table's valid logR window. `Em ~ kappa`, so 1e-30 -> 1e-5 lifts the emission by 1e25, to
~1e10, i.e. 21 orders ABOVE the noise. Absorption and emission then carry the same kappa,
the cell has a real radiative equilibrium, and the noise cannot compete. It stays
physically transparent: at 1e-5 and 1e-22 g/cm^3, tau across the whole ambient shell is
~2e-15.

**Why:** an optically thin cell has no business taking its heating from a flux DIFFERENCE;
the difference is all round-off. **How to apply:** never floor an opacity to zero in a
scheme that deposits `-div F`. Give it a small finite value so the local emission sets the
scale instead. The general lesson is [[red-giant-analysis-opacity-trap]]'s sibling: check
what the *difference* of two large fluxes is worth in units of the cell's energy.

## The semi-implicit change that went in anyway

`two_stream_rt.hpp` `rt_apply` now relaxes toward the true fixed point rather than
linearizing about the current state:
`deq = e((A/E)^(1/4) - 1)`, `x = src*bdt/deq`, `de = deq*(-expm1(-x))`. Exact in BOTH
limits (`src*bdt` as dt->0, `deq` as dt->inf), monotone, cannot cross equilibrium.
`problem/rt_semi_lin` recovers the old form. The user's own commented-out Newton-Raphson
block in the picket-fence path (`two_stream_rt.hpp` ~1888) is the BETTER endpoint: exact
backward Euler on `de/dt = A - E(T)`, general-EOS branch included, no `e ~ T` assumption.
Its `T < 0` branch silently discards the source, which wants a bounded explicit fallback.

## RESOLVED 2026-09-09 (three pieces, all verified by agents)

1. **Direct source** now lives in all three sweep kernels (see
   [[two-stream-rt-three-kernels-trap]]); equals -(Ft-Fb)/dx to 1e-15 in the star; eint
   bit-identical over a 3e5 s regression with the switch off vs on.
2. **The REAL root cause of the fast heating was the opacity LOOKUP, not the arithmetic.**
   The stellar table's density axis stopped at 1e-14; every consumer (pgen KappaTab,
   conduction module, grey two-stream RosselandTable) CLAMPS rho to the table edge, so a
   cell at 2e-22 (ambient) or 1e-18 (thin atmosphere) was evaluated at rho = 1e-14 where
   the node's logR is INSIDE the valid window and the value is real molecular opacity
   (1e-3..1e-2), not the 1e-5 floor -- dtau 1e2-1e3x too large.  The startup guard used
   the TRUE logR and called those cells transparent; the lookups never saw the density.
   Fixed data-side in `ReadOpacityTable`: 200 rows prepended down to log10 rho = -24 at
   log10(opac_floor); nD/ld0 are OUT params so all consumers follow.  Background now holds
   410.8 K at t = 3000 s (11 K = compressional heating of the infall; was 1837 K), star
   bitwise unchanged wherever rho > 1e-14, R0 1-cycle regression identical to the digit.
3. With <hydro>/wb_rmax ([[red-giant-wb-kills-ambient-medium]]) the ambient medium
   free-falls at g and is otherwise inert -- the configuration the user asked for.
