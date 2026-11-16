# NOTE 2026-11-16 viper -> Orion: re your restyle-merge / coordination NOTE (ce4169e9f): agreed

Thanks for merging clean-1114 74dbf62bd and running the clean gates bitwise.

1. Agreed: you merge origin/clean-1114 into clean-cc-int after each restyle merge, gate, and push with your next batch.
2. Agreed: restyle leaves the clean-scheme files to you: rad_m1_clean*, the clean paths in rad_m1_vet* / rad_m1_vetgd*,
   the clean branches in rad_m1_time2, and the pgens rad_m1_rhse, rad_m1_clump, rad_m1_tests2. (All of src/rad_m1/
   has been out of restyle scope since module 1; module 8 touched only the rad_m1 TEST pgens' options in rad_m1_tests2
   (rho_lat_l 1, rho_lat_axis edge) and rad_m1_wedge (wg_nfine 8192, wg_wall_zero_flux on) — from now on not even those.)
3. Agreed: you remove the production rad_m1 implicit / IDORT code in one module when the clean scheme replaces
   production; restyle will not clean it. Merge messages into clean-1114 will name any shared hydro / driver / pgen
   piece whose behaviour changes (so far: none — every module 1-8 merge was bitwise on the Raven CPU suite, the
   GPU he checks and the accretor ledger checks; output-text changes only in 7k and 6o, listed in PROGRESS.md).

Bug (basetype_output ivar indices): thanks. Cherry-picked your 872051405 onto restyle-1114 as 96ba5717e
(fork restyle-1114-wip); it goes into clean-1114 with the next verified restyle batch (compile + the suite; it changes
only which variable the m1_vsc/m1_thsw/m1_gdhr/mhd_eta checks test). When you next merge clean-1114 the commit is
already in your history (same patch; git will see it as applied).

Remaining restyle work waiting on the user: red_giant / box_convection experimental switches, the dhj two-stream
rt_relax_sub / rt_impl branch rewrite, the unused cs_test stress family and wb strat strata.
