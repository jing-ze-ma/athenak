# NOTE from the Raven coordinator: FIN3 is the production method (all sites), final configuration

USER: "make that the production and let other machines know".

## FIN3 = the production method
- Code: fin3-1012 (= rotavg-rst-1012 @ 1fdbe56e plus one default change, hash to follow in an update of this NOTE).
- What it contains:
  - FIN2 keys;
  - the vet_gd twin state from the mean transmission (vet_gd_twin_state = trans);
  - closure averaging over the last 4 rotated direction sets (vet_rotate_avg = 4);
  - vet_interp = auto;
  - the closed-shell cut fix (vet_col_lat_cut_thin = 1);
  - vet_sc rotation on the shared clock;
  - transverse vacuum BCs;
  - bitwise restarts: the conservative2 slope, the averaging history and the Picard window are stored in the rst;
  - link-end outputs keep the output schedule.
- CHANGE since the Caltech TASK of earlier tonight: `implicit_predictor_frozen = false`.
  - The frozen predictor (an item of the sc/gd unification) cost Cartesian time-convergence order (an iteration error at tol
    1e-12), so it was dropped.
  - It is inert on spherical production decks (being verified bitwise).
  - Every production deck NAMES it (false). The Caltech deck docs/handover/prodswitch-1013/bsgP_seed_fin3.athinput.in is updated
    (d22793aa).
- Raven (bsgP, bsgP2n, AG Car E1) and Viper (bsgP, AG Car E1) production run FIN3; nothing has started there yet.

## Per site
- **Caltech** (TASK-2026-10-13-caltech-bsg-fin3.md):
  - use the UPDATED deck (predictor_frozen = false);
  - if you have not built yet, build the final hash given in the update of this NOTE (within ~1 h);
  - if you already built 1fdbe56e, keep it: with the updated deck the results are identical.
- **DeltaAI**: NOT affected. Your BSG run stays on the old method. Only TASK-2026-10-13-deltaai-bsg-d3rst.md applies (restart-fixed
  old-method binary d3rst-1013 @ 9ceb9cdc, swapped between links after the scratch restart check).

## UPDATE: final production hash
**fin3-1012 @ 281dd686** (pushed) = rotavg-rst-1012 1fdbe56e + `implicit_predictor_frozen` default false. Raven/Viper production
is being re-pointed to it after a bitwise check against 1fdbe56e with the same decks.
Caltech: build 281dd686 (same flags as in your TASK) with the updated deck. If you already built 1fdbe56e, it gives identical
results with the updated deck; rebuilding is preferred, so that all sites run the same hash.
