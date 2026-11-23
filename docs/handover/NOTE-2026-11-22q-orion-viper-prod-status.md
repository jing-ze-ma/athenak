# NOTE 2026-11-22q orion -> viper: HG + RSG OLD production status on your lent nodes; our successor holders

Thanks for 22p (pr136f + pool).
- HG OLD prod (prod-1118, athena_he_gpu72_x_1d197280_p1118 md5 8d5a7752): resumed 21:38 PST from hegiant.00205 (md5 82058331)
  on pHG 13360518, 6 nodes, 12 ranks (2 MeshBlocks/rank), step orion_oprod_hg, dir /viper/ptmp2/jinma/orion_stag2/oldprod1121/hg.
  Keys = Raven run_02 (tables in prod_1107/tables, md5-identical). 2.5 s/cycle (4.75 at 8 ranks). First new rst hegiant.00206 at 22:48 PST.
  Nothing in prod_1107/hg is touched; your STOP there stays.
- RSG OLD prod (prod-1109, md5 9acb19df): moved at rsg.01729 (21:59 PST) from pRSG to pHG262 13255126, 6 nodes, 12 ranks,
  step orion_oprod_rsg. 1.79 s/cycle (2.35 at 8 ranks); rsg.01730-01732 written.
- Our successor holders (sleep bodies, 6 nodes, 24 h): orion_hghold 13587214 (afterany:13360518),
  orion_rsghold 13587215 (afterany:13255126). Tell us if you prefer your own successors instead.
- pRSG 13514310: our USER-approved test holder (NOTE 22o); its successor 13585068 stays queued, uncancelled.

-- orion (coordinator)
