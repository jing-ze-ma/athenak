---
name: rmdx-analysis-0928
description: 09-28 Rm_dx analysis of prod4 (eos eta unresolved at grid scale everywhere) -- user judged it too pessimistic; two follow-ups DEFERRED (drag time vs cap, effective numerical eta)
metadata:
  type: project
---
/viper/ptmp2/jinma/rmdx_0928/RESULTS.md: prod4 rot 88, median Rm_dx = |v| dr/eta >= 10^1.45 everywhere; closest to
resolved = nightside 2-10 mbar and base 60-500 bar; Ohmic heating 1.5e-4 of stellar input. WASP-121b 1x estimate
(hydro profile + prod4 |B|): Rm_dx <= 1 at 0.03-1 bar (cold, metals condensed) but eta at the max_eta = 5e12 cap in
33-66 % of those cells; dx^2/eta only 5 dt there.
**User 09-28: "that seems very pessimistic".** Agreed flaws: grid-scale Rm_dx is the wrong test (high Rm on the
dayside is physics; what matters is numerical dissipation at energy-containing scales); bulk |v| overstates (shear
Rm_dv 10^0.7-10^4); "the cap decides" unproven (depends on drag time eta/v_A^2 vs advection).
**DEFERRED follow-ups (user: "not now but keep in mind"):** (1) drag time eta/v_A^2 with capped and uncapped eta vs
advection/rotation time; (2) effective numerical eta per layer and scale from the magnetic energy budget between
dumps (induction work - dE_mag/dt - explicit Ohmic); later a resolution comparison of large-scale quantities in the
WASP-121b MHD phase. Raise these when setting up the WASP-121b MHD phase. Related: [[next-prod-ck-c2]].
