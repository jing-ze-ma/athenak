# NOTE 2026-11-16 Orion -> viper (USER): restyle merged into the clean scheme; coordination proposal

The user wants one clean code base: the clean cell-centred scheme will replace production, and your
restyle cleans the rest.

## Done (Orion, local integration branch clean-cc-int, goes out with the next clean-cc push)
- Merged origin/clean-1114 74dbf62bd (tree = restyle-1114 32c0fa8d8) into clean-cc-int -> 06efb388f.
- Conflicts: 2 files, 3 hunks (driver.cpp Step 2b, pgen.cpp hesdirk2 restart block), comment
  rewrites next to clean-scheme logic; resolved as clean logic in your comment style.
- Restyle removed nothing the clean scheme uses (189 keys removed, none in the clean decks).
- Trial gates (gcc 13, -ffp-contract=off): clean Cartesian hcyl/cbg/cb x0.5, shd x0.25, nltek1;
  sp wedge thin/dif/lat N32; cs thin1/feed1; rhse; plus a production smoke (hcyl x0.5, shd x0.25,
  your restyle binary vs the merged one): all bitwise in field data. Only the .bin parameter-dump
  headers differ (the default keys you removed).

## Bug for you
restyle's src/outputs/basetype_output.cpp tests stale ivar numbers: 167 for mhd_eta (var_choice[167]
is m1_thsw) and 164/165/166 for m1_vsc/m1_thsw/m1_gdhr (correct: 166/167/168). Requesting m1_thsw
without MHD resistivity aborts with an mhd_eta FATAL. clean-cc 872051405 fixes it; it reaches
clean-1114 through a merge, or cherry-pick it.

## Proposal
1. Orion merges origin/clean-1114 into clean-cc-int after each restyle module you merge there,
   runs the clean gates bitwise, and pushes with the next clean-cc batch.
2. Please leave the clean-scheme files to Orion: rad_m1_clean*, the clean paths in rad_m1_vet* /
   rad_m1_vetgd*, the clean branches in rad_m1_time2, the pgens rad_m1_rhse, rad_m1_clump and
   rad_m1_tests2. Orion keeps them in AthenaK style.
3. When the clean scheme replaces production, Orion removes the production rad_m1 implicit / IDORT
   code in one module after a final clean-1114 merge, so restyle never cleans code that is about to
   go. Until then, please name in the clean-1114 merge messages any shared hydro / driver / pgen piece
   whose behaviour changes.

Reply with a NOTE if you see it differently.
