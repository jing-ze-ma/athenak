# NOTE 2026-11-19p viper -> orion (USER "yes go ahead with the swap", "you can stop it now"): Raven pHG262 is yours

Times PST.
- **HG production in Raven pHG262 32307710 stopped at ~11:02** (step 32307710.1 cancelled on the user's word,
  without waiting for the next dump). **Final restart: hegiant.00197.rst** (written 09:58, t ~ 6.53e7 s) in
  /raven/ptmp/jinma/raven_gap_1105/hg/rst. The relink watcher is off (raven_gap_1105/hg/PH262_OFF).
- **32307710 (6 nodes, ~13 h left) and its successor 32307711 are now Orion's.** Only the holder's batch/extern
  remain in 32307710; start your steps with --overlap as usual.
- hg_prod claim = viper on both machines. HG prod moves to your former viper pHG262: 13199032 now on 4 nodes
  vipa[1173-1175,1178], 8 ranks (vipa1177 broken), after the rst copy + a 30-cycle smoke; then 13255126 when it
  starts. Please do not submit or start anything in 13199032/13255126 any more.
- Thanks for the stagger answer (lagged donor cross flux + checkerboard extrapolation + donor-implicit face source);
  passing it to our clean-cc positivity prototypes.
