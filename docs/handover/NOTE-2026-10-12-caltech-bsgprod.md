# NOTE 2026-10-12 from Caltech: bsgprod (REAL BSG true reproduction, final method)

Re: TASK-2026-10-12-caltech-bsgprod.md (a450ee38 + deck fix 53995cf9).

## QUEUED (10-12 PDT) -- not started
- Dir /resnick/groups/carnegie_poc/jingze/bsgprod_1012; deck = bsgP_seed.athinput.in with only @FILES@ filled (IC 246a666e, rosseland 450bc0c1, planck e18ca805 checked).
- Jobs: build 4401747 (d3bf42d7 -> athena_gpu_he_star_m1_d3bf42d7_nofma) -> smoke 4401748 (25 cyc, ndiag 1, 2 H200; fails on FATAL/NaN/NON-CONV/dt outside 80-95 s) -> prod links 4401749 (afterok smoke) -> 4401750 -> 4401751 (12 h, 1 node x 2 H200, 2 ranks x 8 MB).
- H200 est. start ~10-14 20:23 PDT, so Raven/Viper will likely start first: tell us (NOTE line "STARTED" on your side) and we cancel ours. STARTED section follows here if link 1 runs first.
