# NOTE 2026-11-15 Orion -> viper: cc_idort_head follow-up (hcyl x2) + clean-cc progress

## hcyl_FIN3_g8x32 x2, arm CC1, dc67dba7b (OMP 28 threads, wall 9698 s)
score 0.1681 (x0.5 0.0468, x1 0.2088): BE on 122/146 CLEAN lines, max nreal 101184, 0 dt_retry, 0 CLEAN WARNING.

## Progress on Orion (local branches off clean-cc dc67dba7b; integration branch clean-cc-int, not pushed yet)
- clean_closure = extrap1: ONE formal solution per step, order 2.00/2.00 in time (cell, tau 10); identical results to extrap
  on every IDORT deck where stages are BE (expected), 2-3.5x faster on nltek/hcyl.
- ApplyClosureLimits no longer clips the clean solution: the clip removed 1.3-2.8e-2 of the radiation x-momentum and created
  1-6e-3 of the energy per BE step on cb/cbg/ws66 (no gas compensation) and made K1 inconsistent with U^n.
- Consistency (HOLO) term clean_consistency: face flux F_scheme(U) + [F_ray,face - F_scheme(J, H)] with ray face fluxes taken from
  the short-characteristic sweep itself (x1: intensity at the face crossing; x2: exact transverse-shift flux; open x1 faces:
  entering + outgoing rays), weighted by the rays' C1. vs CC1: cb x0.5/x1 1.686/1.875 -> 0.460/0.953, cbg 1.624/1.631 ->
  1.308/1.423, ws66 1.345/1.447 -> 1.167/1.334, hcyl x0.5/x1 0.0468/0.209 -> 0.0191/0.0137, shd/shdg x0.25 0.264/0.256 ->
  0.0225/0.0332, shd/shdg x0.5 now finish (0.0157/0.0302; CC1 BE non-convergence), nltek1 x1 9.115e-5 (= CC1), rw unchanged,
  energy/momentum round-off; every BE step ends realizable. Stage-1 solves on the beam decks are still mostly rejected
  (frozen rays' D as a linear stress at CFL 256; diagnosis on Orion).
- Policy (user): cell layout only (staggered frozen); all tests at the production tolerance with larger rw amplitudes.
