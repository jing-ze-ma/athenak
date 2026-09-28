---
name: gate-on-wasp121b-not-prod4
description: user 09-28: WASP-121b is the production candidate -- dhj gates, repro checks and timings use the WASP-121b inputs/restarts, not prod4 (prod4 only where an MHD dhj state is needed, and say so)
metadata:
  type: feedback
---
All dhj checks (bitwise gates, determinism, timing, restart gates) use the WASP-121b setup: w121prod 1x/10x inputs
and restarts (/viper/ptmp2/jinma/w121prod_0927/w{1x,10x}/rst, read-only; later the fresh 09-28 runs), with the
production key set (T4, closed wall, ...).
**Why:** user 09-28 "I don't see why you checked prod4. wasp121b is the latest candidate for production" -- the
ck non-determinism was first reported from a prod4 restart that ran an old non-T4 Jacobian path, which misled.
**How to apply:** write WASP-121b into every worker brief; use prod4 only when an MHD dhj state is required
(WASP-121b has no MHD run yet) and label it as such.
