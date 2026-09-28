---
name: speed-work-20pct-rule
description: user decision (Caltech 09-27, docs/handover SESSION-2026-09-27-caltech.md sect. 6): speed work only if it gains > 20 % of the WHOLE run; ck speed work closed; micro-levers declined
metadata:
  type: feedback
---
Speed work is worth doing only when it gains more than ~20 % of the whole run's cost. ck/hydro micro-levers (split
ConToPrim tasks, PackAndSendCC vector length, lin_sum fold, single EOS solve) were declined; ck speed work closed.
**Why:** user decision recorded by the Caltech session 09-27 (11c9a5be).
**How to apply:** do not propose or pursue optimisations below ~20 % of the whole run; for switches that change
results for small gains, decide by accuracy (noise test) and make them opt-in. Combine with
[[no-accuracy-sacrifice]] and [[measure-cost-on-gpu]]. Timing runs must use ck_impl_verbose = false (+1.3 ms/cycle).
