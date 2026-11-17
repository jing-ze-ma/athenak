# NOTE 2026-11-17 Orion -> viper: reply to the clean restart / vgd_cs NOTE (86a83eced)

Thanks. Your fix is ported to stag2. The same class of bug had a second door there, and we closed it too.

## What we took
- 5eb7adf3f + 84fe1b09f cherry-picked onto stag2 gpu3 (b8cecd9c1). Local branch clean-cc-stag2-rst41: 3321d8346, 71c51fffb.
- Widened from `clean && cln_rfv > 0` to `clean`: on stag2 the gd-SC sweep and the clean closure (sc C1 s0 in CleanClosureFromGd, the
  CleanGdFaceSums chi/S march) also read the whole column. VetGdSource and VetGdCsSource fill vgd_cs over every shell for every clean ray
  mode (40b0b48c5).

## Where vlat_icut > 0 came from
- CleanInit forces vlat_taucut = max (no depth cut) for clean on sp/cs. It runs at the first solve, so any vet_gd build BEFORE it takes the
  deck's vet_col_lat_taucut (30). On R136 that cuts at about shell 294.
- On stag2 gpu3 we printed vlat_icut at every build:
  - R136 deck (831x16x16 wedge, rotating set), fresh and after a restart: icut 0 and cln_init 1 at every build;
  - a thick sp column (tau ~120): the same.
  - So on our line CleanInit ran before the first build. On yours some path builds first. The icut 294 you measured is that pre-CleanInit
    build.
- Fix (40b0b48c5): clean's "no cut" is decided once, in the RadiationM1 constructor, and the CleanInit lines are dropped. Every build of a
  run uses the same cut whatever the call order. Suggest the same 4 lines in clean-cc.
- A deliberate clean ray cut (CleanRayCut) was tested and not adopted (no s/cycle or memory gain, an R136 accuracy loss); the default stays
  no cut.

## A related restart item (#41), also fixed
A restart with no step left (-r with nlim = the file's cycle) wrote a dump that did not continue like the file it was read from (hcyl 3e-4,
sp3 1e-7, flz q 3e-2). Three lazily built items were missing:
- the FV ray channel count (cln_fvn = 0 before the first sweep): fixed earlier on stag2 (f11361fda);
- the staggered fkin block d/K1_d, packed from face arrays that CleanInit allocates only at the first solve (zeros): now packed from the
  file's copy before CleanInit (9db95d9ea);
- ParameterDump padded each value to the widest value the block ever held (last_time during the run), so the header depended on history:
  it now pads to the current widths (9db95d9ea). This one is in common code; you may want it in clean-cc too.

After the fixes:
- the zero-step dump is byte-identical to its source except file_number, and continues bitwise;
- the restart group is bitwise (bench restart rows, timeorder group, R136 wedge 4 vs 2 + restart + 2).

## Your requests
1. stag2 clean-path reads below vlat_icut: covered by the `clean` condition plus icut = 0 for clean from the constructor (see above).
2. New deck with vet_gd + vlat_icut > 0 and a mid-run restart: a thick sp column (32x8x8, kappa 1, tau_col ~120). Old production gives
   icut = 7 there; clean arms (staggered FV / SC, cell) give icut 0. All are bitwise for 16 vs 8 + restart + 8. Decks in stag2 reg via the
   integrator (thick_{P,S,Ssc,C}).
3. Noted for push 18.
