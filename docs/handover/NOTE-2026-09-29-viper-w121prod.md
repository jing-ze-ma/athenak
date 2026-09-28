# NOTE 2026-09-29 (viper -> Caltech): WASP-121b fresh start is ALSO queued on viper (user 09-29)

User 09-29: "it's not so clear if caltech can actually start the production runs faster. let's just queue the 1x
and 10x production here also, each on 2 gpus." This supersedes "viper must not launch it" in
NOTE-2026-09-28-caltech-w121prod.md: both machines run the same fresh start; the user decides later which to keep
(or keeps both as a cross-machine check).

- viper run dir `/viper/ptmp2/jinma/w121prod_0929/{w1x,w10x}`, binary = rt-integration cda4da33 (build_inc_viper
  dhj_gpu, ROCm 6.3, md5 b049e28e).
- Inputs = w121prod_0927 + `rsolver = lhllc`, 1x `ck_impl_maxit 16`, 10x `ck_impl_maxit 24` + `ck_impl_dtmax 0.25`,
  `ck_impl_xstep 8`, `flux_hst_floor = true`; current IC, closed wall, top sponge on / bottom off; tlim 300 rot.
  Same settings as the restarted Caltech runs (NOTE-2026-09-29-dhj-rsolver.md).
- Jobs (apu, 1 node x 2 GPUs each): 1x 12018343 (+12018344 afterany), 10x 12018345 (+12018346 afterany); queue
  estimate ~19 h at submission (apu saturated).
- Smoke (apudev 12018342, 50 cycles both arms): clean — rc 0, no FATAL/NaN, 0 NOT-CONVERGED, Etot_bot = Lrad_bot,
  Mdot_bot = 0, Efloor/Mfloor/Efloor_rt columns present.
