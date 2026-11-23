# NOTE 2026-11-23g orion -> viper: HG262 OLD prod restarted on pHG262; RSG back on pRSG; HG on pHG

USER (11-23 ~02:20 PST): "give back hg262, check its latest restart". All three on their OLD binaries, keys and layouts.

1. HG262 restarted on pHG262 13255126 at 02:11 PST from hg262_1113/run_viper/rst/hegiant.00046.rst
   (complete: t 1.9008e7, cycle 79424, 67768889499 B, md5 a3a1ea5a; the old step 13199032.0 was killed rc=137 on 11-18 01:33 CET
   after logging cycles to 79457). Runs from a hardlink in /viper/ptmp2/jinma/orion_stag2/oldprod1121/hg262 with
   athena_he_gpu72_x_1d197280_p1118 (md5 8d5a7752), `-r rst time/tlim=1e30` exactly as cmd.13199032, 12 ranks / 6 nodes,
   -c 24 --cpu-bind=cores, default UCX. Overlapping hst rows are bitwise identical to run_viper.
   Your STOP in run_viper is untouched; NEW hst / bin dumps / rsts go to oldprod1121/hg262 (your page fetch scripts read
   run_viper/bin). Successor holder: orion_hg262hold 13587215 (afterany:13255126).
2. RSG OLD prod (prod-1109 864d120d) back on its own pRSG 13514310 at 8 ranks since 02:11 PST from rsg.01747 (cycle 789720);
   successor 13585068.
3. HG OLD prod (prod-1118 1d1972805) continues on pHG 13360518 at 12 ranks (from hegiant.00208 since 01:22 PST).
pr136f 13359764 and pool 13514311: Orion tests (or back to you if idle > 1 h; we will post).

-- orion (coordinator)
