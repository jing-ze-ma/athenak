---
name: no-accuracy-sacrifice
description: User rule 09-25 -- speed-ups must not sacrifice accuracy/order; first-order shortcuts (multi-rate radiation, coarse tensor cadence, looser tolerances) are not acceptable as defaults or recommendations
metadata:
  type: feedback
---
User 09-25: "we don't want to sacrifice accuracy" (after m1-fast4 showed multi-rate radiation k=4 gives rad/hydro ~1.0
but is first order).
**Why:** the implicit VET was just made 2nd order in space and time; speed must not undo it.
**How to apply:** only accuracy-neutral levers (bitwise / round-off: kernel spills, fusion, lookups, launches, exact
reuse); anything that lowers order or raises error beyond the solver's own tolerance noise stays opt-in and
unrecommended (e.g. implicit_mr_every, vet_sc_every > 2, lin_tol > 1e-10, ck cadence levers only after an A/B shows
no physical difference). Related: [[m1-precond-mg]], [[ck-nquad2-production]].
Exception approved by the user 09-25: ck_impl_tol 1e-7 for 10x-metallicity dhj (rsec 20 + maxit 12).

**User 09-26: the accuracy bar is the ROUND-OFF spread, not the time-discretization error.** A speed-up passes only if its deviation, in the region that counts, is no larger than the reference's own round-off spread (1 vs 2 GPUs, or a 1e-14 kick). Hence M1 implicit_tol stays at 1e-8 (1e-7 would give -12 % on the box, at 2e-9 = 60x round-off: rejected).
