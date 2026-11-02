# NOTE 2026-11-02 viper -> Caltech (USER): status (halofc MI300A gate, Strang decision) + r1m deck: ck_impl_maxit 16 is stale

1. halofc-1102: merged locally onto fin51 3957c9562 (TU split) as halofcmerge2-test-1102 (head 7f0e3dd9). MI300A gate
   job ran (viper 12574246: r1m cs MHD 1 + 2 ranks, 20+20 vs 40, W121 rot-320 restart, timing); verdicts being
   analysed now. We merge into fin51-r-fix if bitwise and NOTE you the hash. (An earlier cs-hydro smoke FATALed at
   deep_hot_jupiter_rt.cpp:1129 on our hydro deck copy; being checked; likely the deck, not halofc.)
2. Strang: the user has not decided Lie vs Strang yet; we recommended waiting for your rkg_strang_fuse + ck_strang
   arms (please build ck_strang on 3957c9562). r1m stays held for those results.
3. r1m deck: problem/ck_impl_maxit = 16 is a stale leftover. History: the old code default was 8; the 1x/3x fresh-start
   decks raised it to 16; on 09-30 the code default became 24 (merged 9bf6db23: a call stopped at the cap returns an
   unchecked iterate), and the 09-30 A/B (3649979/3649980) showed no dt effect between the key sets. The explicit 16 in
   the deck overrides the default 24. For NEW runs (r1m launch, hydro fresh start) please drop the line (-> 24). For a
   run already going, changing it at a restart would change results only where the cap is hit; leave it unless you
   see non-converged ck calls at the cap in the event log.
