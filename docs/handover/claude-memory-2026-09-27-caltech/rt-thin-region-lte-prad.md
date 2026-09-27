---
name: rt-thin-region-lte-prad
description: 2026-09-13 CHECKED - eos_radiation=true keeps LTE aT^4 (energy) + aT^4/3 (pressure) in the EOS at ALL optical depths, including above the two-stream handover where the two-stream owns the real radiation field; the B-star FeCZ box top (tau 1e-2) has Prad/Pgas 8.6 (input comment said 1.1), radiation energy 17x gas thermal energy, sound speed inflated ~2.5-3x; force (grad Prad / rho g) only 0.10-0.15 there and within 1.8x of kappa F/(c g) -> stratification ~right, stiffness + heat capacity wrong. FeCZ (tau 580-3300, Prad/Pgas 0.4-0.5, diffusive) unaffected. Same in red giant + dhj above their handovers.
metadata:
  type: project
---
Numbers from bench/bstar_fecz/smoke_rt2/column_used.txt with the GS98 X0.7 Z0.014 Rosseland table, g0 1.4454e4, F 4.21e13:
tau 0.01: Prad/Pgas 8.6, gradPrad/(rho g) 0.10, kappaF/(cg) 0.05 | tau 0.1: 0.85, 0.15, 0.12 | tau 2/3: 0.40, 0.23, 0.21 |
tau 100-300: 0.33, 0.25, 0.24 | FeCZ tau 1e3-3e3: 0.40-0.54, 0.35-0.39, 0.35-0.38. Suspect (unproven) for the smoke_rt2 dt fall
1.49->0.32 s and the top-cell floor/fofc/tclamp counters. Options (user decision): cut the box top at tau 1-3 instead of 1e-2,
or a separate radiation energy above the handover (weighting the table term by the blend w is not exact: the two-stream would
have to own Prad + its force). See [[fecz-box-projects]], [[rt-handover-semi-implicit-bug]].
