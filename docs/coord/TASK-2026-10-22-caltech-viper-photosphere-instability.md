# TASK 2026-10-22 caltech -> viper (USER: TOP PRIORITY): fix the photospheric grid-scale instability (seam/panel imprint) in the global decks

USER: finding the fix is now the top priority (it affects RSG E and HG global; both productions would restart on the fix).
Caltech's queue is slow for many 1-GPU arms -> please run arms on viper in parallel. Split to avoid duplicate work:
**viper = RSG E deck; Caltech = HG deck** (Caltech arms listed below so you don't rerun them).

## What Caltech found (HG global, f55be371, cs N=128, nx1 704)
- The T_eff panel imprint is real (E and F at all r > 1.003 R_ph show face-centre minus edge = -4 % at 2 d; it is exactly 0 at 1 d).
- Source: shells 0.990-0.996 R_ph at the He kappa_R cliff (dlnk/dlnT ~ 15, about 1 radial cell thick). From 1 to 2 d these shells
  turn two-phase: ~50 % "cold" cells, T sd/mean 0.8, near-checkerboard, dT and drho anticorrelated (isobaric), tiny velocities.
  The cold fraction depends on cell size (0.21 at face centres vs 0.49 near corners) -> L 0.98 vs 1.01 -> panel imprint.
- N=32 1-GPU arms (photospheric E sd/mean at 1.0 -> 1.1 d; e-fold ~0.3 d):
  - NO effect: vet_gd_rotate_dt 0 and 2.2e4; implicit_face_opac_n off; implicit_opac_newton off; vet_col_source relaxed;
    cfl 0.15 / 0.6; time_scheme be; vet_gd_replace off; lat_offdiag none; fofc off; am_transport off; cs2_clip gas; cs2_vimp off.
  - HLLE: grows 5x faster.
  - nx1 1408 (2x radial): sd 1e-5 and smooth, ~1000x lower; panel L contrast 2.7e-7 vs 2.7e-5.
- Running/queued at Caltech (HG): nx1 1408 long; nx1 704 and 800 with photosphere-packed spacing; nx1 1056; time2_vet_col = rebuild.

## RSG E (Caltech prodF, f55be371, N=128, from t=0)
T_eff face/edge contrast is +0.02 % up to 39 d, then +0.14 % at 49 d, +0.7 % at 59 d, +1.35 % at 69 d.
T_eff sd grows from 0.4 % to 5.5 % at 79 d, and R_ph sd from 1e-5 to 1.3e-3. The sign is opposite to HG (face centres hotter).
Onset ~49 d = when dt time2_reject retries began. The phot band (0.9-1.2 R_ph) shows cell-scale noise in a thin layer at 0.99-1.0 R_ph.
Your Raven rsg_g25 N=128 to 32 d was quiet (731b45c5): it does not test the post-49 d state.

## Asked of viper (RSG E deck, f55be371 build from 9c2ca92c)
1. Find the H/H- kappa_R cliff at the RSG E photosphere: dlnk/dlnT, and its width in radial cells at nx1 of the prod deck.
2. Reproduce the onset cheaply: N=32 (or 48), fresh, to ~60 d (or until T sd > 1 %). If the onset is far too slow, say so and
   propose a shorter route (e.g. restart from a perturbed state).
3. Compare grids at the same N: base nx1; 2x nx1; cells packed near R_ph (local x2-4 over ~0.97-1.03 R_ph) at the same nx1, if
   the x1 grid supports it (if not, say what is missing; Caltech may add an opt-in key).
4. Report: photospheric T sd/mean; panel contrast (face centre |x|,|y|<0.3 vs edge >0.8); dt and retry count; s/cycle ->
   cost per sim day vs base. Recommend the cheapest grid that suppresses the instability.
Rules: opt-in keys only, default bitwise; incremental builds; NOTE results back to bsg-files-1009.
