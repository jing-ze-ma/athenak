# NOTE 2026-11-24w caltech -> orion: layout batch (24v) for t9l62 (02-17) and M35ZAMS (01-15); arms back to back, one node each (COORD: no idle node; RSG has 7333072)

- Fresh arms use deck copies with the production keys WRITTEN IN (fofc = true, fofc_retest = true, clean_stag_mg_levels = 4): a fresh start FATALs on absent command-line keys. Files: prod_rc3_1124/layopt/decks/.
- Peak GPU memory from memwrap's 2-s nvidia-smi log (mem.log), not a separate mem.csv.
- **Z0**: ARGS: -i prod_rc3_1124/layopt/decks/IA123_m35zams_lay.athinput time/nlim=250 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30
  2026-11-24 01:39:38 rc=0 wall=131 s last: elapsed=1.204216e+02 cycle=250 time=8.701267e+02 dt=3.480506e+00; s/cyc 0.3994 (c50-250); setup elapsed=2.264009e+01 (+ total wall=131); peak GPU mem 0,39077 MiB 1,39077 MiB 2,39077 MiB 3,39077 MiB 
  BiCGStab it 6.39 / stage line, solves 2.11 / stage line (502 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user 0afded9484f46b23ffbf25ae5432cd15, hydro ba148c66be1f64002ca26cacb521474a; dt at c200 3.48051
- **Z1**: ARGS: -i prod_rc3_1124/layopt/decks/IA123_m35zams_lay.athinput time/nlim=250 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 meshblock/nx2=64 meshblock/nx3=64
  2026-11-24 01:41:39 rc=0 wall=104 s last: elapsed=9.403978e+01 cycle=250 time=8.701267e+02 dt=3.480506e+00; s/cyc 0.3225 (c50-250); setup elapsed=1.553198e+01 (+ total wall=104); peak GPU mem 0,19937 MiB 1,19937 MiB 2,19937 MiB 3,19937 MiB 
  BiCGStab it 6.39 / stage line, solves 2.11 / stage line (502 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user 0afded9484f46b23ffbf25ae5432cd15, hydro ba148c66be1f64002ca26cacb521474a; dt at c200 3.48051
- **Z2**: ARGS: -i prod_rc3_1124/layopt/decks/IA123_m35zams_lay.athinput time/nlim=250 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 meshblock/nx2=128
  2026-11-24 01:43:38 rc=0 wall=107 s last: elapsed=9.735580e+01 cycle=250 time=8.701267e+02 dt=3.480506e+00; s/cyc 0.3329 (c50-250); setup elapsed=1.622143e+01 (+ total wall=107); peak GPU mem 0,22747 MiB 1,22741 MiB 2,22741 MiB 3,22747 MiB 
  BiCGStab it 6.39 / stage line, solves 2.11 / stage line (502 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user 0afded9484f46b23ffbf25ae5432cd15, hydro ba148c66be1f64002ca26cacb521474a; dt at c200 3.48051
- **Z3**: ARGS: -i prod_rc3_1124/layopt/decks/IA123_m35zams_lay.athinput time/nlim=250 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 meshblock/nx2=64 meshblock/nx3=64 mesh/layout_tile=64
  2026-11-24 01:45:36 rc=0 wall=111 s last: elapsed=9.574936e+01 cycle=250 time=8.701267e+02 dt=3.480506e+00; s/cyc 0.3237 (c50-250); setup elapsed=1.610387e+01 (+ total wall=111); peak GPU mem 0,19987 MiB 1,19987 MiB 2,19987 MiB 3,19987 MiB 
  BiCGStab it 6.33 / stage line, solves 2.17 / stage line (502 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user 1e0c57b57e944bf7f534fe471424f9d2, hydro 7e8b70e0f8f80341fba8ebbf18fad8e5; dt at c200 3.48051
- **Z0r**: ARGS: -i prod_rc3_1124/layopt/decks/IA123_m35zams_lay.athinput time/nlim=250 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30
  2026-11-24 01:48:01 rc=0 wall=133 s last: elapsed=1.220681e+02 cycle=250 time=8.701267e+02 dt=3.480506e+00; s/cyc 0.4055 (c50-250); setup elapsed=2.312683e+01 (+ total wall=133); peak GPU mem 0,39077 MiB 1,39077 MiB 2,39077 MiB 3,39077 MiB 
  BiCGStab it 6.39 / stage line, solves 2.11 / stage line (502 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user 0afded9484f46b23ffbf25ae5432cd15, hydro ba148c66be1f64002ca26cacb521474a; dt at c200 3.48051
- **Z4**: ARGS: -i prod_rc3_1124/layopt/decks/IA123_m35zams_lay.athinput time/nlim=250 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30 meshblock/nx2=64
  2026-11-24 01:50:14 rc=0 wall=115 s last: elapsed=1.050396e+02 cycle=250 time=8.701267e+02 dt=3.480506e+00; s/cyc 0.3536 (c50-250); setup elapsed=1.888350e+01 (+ total wall=115); peak GPU mem 0,27763 MiB 1,27753 MiB 2,27753 MiB 3,27763 MiB 
  BiCGStab it 6.39 / stage line, solves 2.11 / stage line (502 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user 0afded9484f46b23ffbf25ae5432cd15, hydro ba148c66be1f64002ca26cacb521474a; dt at c200 3.48051
- **T0**: ARGS: -i prod_rc3_1124/layopt/decks/IA123_t9l62_lay.athinput time/nlim=250 output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30
  2026-11-24 01:50:39 rc=0 wall=792 s last: elapsed=7.817418e+02 cycle=250 time=2.421289e+04 dt=9.684980e+01; s/cyc 2.8459 (c50-250); setup elapsed=1.032968e+02 (+ total wall=792); peak GPU mem 0,106529 MiB 1,106365 MiB 2,106365 MiB 3,106529 MiB 
  BiCGStab it 4.07 / stage line, solves 2.36 / stage line (502 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user b578af5a0b7ae6839f552b4dbc336169, hydro 44d04279b221b4c4d52b07f8281723fc; dt at c200 96.8505
- **ZD0**: ARGS: -r prod_rc3_1124/m35zams/rst/m35zams.00003.rst time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/floorcap_v=3e7 time/nlim=10708 job/final_restart=false output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output
  2026-11-24 01:52:34 rc=0 wall=134 s last: elapsed=1.197647e+02 cycle=10700 time=2.691706e+04 dt=2.478668e+00; s/cyc 0.3645 (c10410-10700); setup  (+ total wall=134); peak GPU mem 0,39019 MiB 1,39019 MiB 2,39019 MiB 3,39019 MiB 
  BiCGStab it 6.42 / stage line, solves 2.80 / stage line (600 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user eecc3f5cea010e44a616fb17fe45f6f6, hydro 7e2f297369157c33669d14e272707f05; dt at c200 2.47867
  DTREP/DTREG (180 lines; first 2): `DTREG rest  cycle=10410 npol=7 cfl*cell_dt=2.385996e+00 row=41 theta_deg=8.935453e+01 phi=3.750000e-03 r=5.561691e+11 dir=1 rho=1.000000e-17 T=2.341378e+13 v=-2.781347e+07 cgas=6.246836e+06 ceff=6.279103e+06 h=2.711491e+08 cfl*cell_dt_gas=2.388256e+00|DTREP cycle=10410 time=26198.2 cfl*dt_hyd=2.38254 cfl*dt_gasonly=2.38441 rank=2 gid=10 kji=22,9,336 dir=1 r=5.5644e+11 rho=1e-17 T=3.22261e+13 pg=0.000322261 G1=1.66667 cgas=7.32873e+06 v=-2.68032e+07 ceff=7.35554e+06 h=2.71281e+08 tauT_h=9.50329e-10 Erad=9270|`
  floorcap_report: 325 lines; last: `CLEAN FALLBACK cycle 10707: pos_gas 0 (from E 0.000e+00) e_floor 0 (from gas 0.000e+00, column 0.000e+00) ffct_left 0 auto_be 0 floorcap v 1518978 T 1089670 c2p_dfloor 10 (mass 2.265e-07) floor_pay 7.233e+23 (unpaid 0.00`
- **ZD1**: ARGS: -r prod_rc3_1124/m35zams/rst/m35zams.00003.rst time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/floorcap_v=3e7 time/nlim=10708 job/final_restart=false output2/dt=1e30 output3/dt=1e30 output5/dt=1e30 output
  2026-11-24 01:55:34 rc=0 wall=148 s last: elapsed=1.283138e+02 cycle=10700 time=2.720859e+04 dt=3.480506e+00; s/cyc 0.3915 (c10410-10700); setup  (+ total wall=148); peak GPU mem 0,39123 MiB 1,39123 MiB 2,39123 MiB 3,39123 MiB 
  BiCGStab it 8.14 / stage line, solves 2.94 / stage line (600 lines); CW 0, FATAL 0, fallbacks 0, DT_RETRY 0
  md5 rows 0-200: user 41ec86b7bee60c44bb4ffb70474e3e5e, hydro bd43fe29174f005e44c44346a1011933; dt at c200 3.48051
  DTREP/DTREG (180 lines; first 2): `DTREP cycle=10410 time=26199.2 cfl*dt_hyd=3.48051 cfl*dt_gasonly=3.87656 rank=2 gid=10 kji=24,20,3 dir=1 r=4.73055e+11 rho=0.000183423 T=1.91129e+14 pg=3.50575e+10 G1=1.66667 cgas=1.78479e+07 v=0.323046 ceff=1.98789e+07 h=2.30629e+08 tauT_h=25959.2 Erad=3.1621|DTREG polar cycle=10410 npol=7 cfl*cell_dt=3.48051 row=2 theta_deg=88.2357 phi=0.04275 r=4.73055e+11 dir=1 rho=0.000183423 T=1.91129e+14 v=0.323571 cgas=1.78479e+07 ceff=1.98789e+07 h=2.30629e+08 cfl*cell_dt_gas=3.87656|`
  floorcap_report: 325 lines; last: `CLEAN FALLBACK cycle 10707: pos_gas 0 (from E 0.000e+00) e_floor 0 (from gas 0.000e+00, column 0.000e+00) ffct_left 0 auto_be 0 floorcap v 1621016 T 1032036 c2p_dfloor 6 (mass 1.220e-07) floor_pay 0.000e+00 (unpaid 0.000`
- ZD2 skipped (24zb: start Z3 now).

### m35zams: FRESH restart with the Z3 layout (orion 24zb, USER): 4 blocks of 384x64x64, layout_tile 64, on hpc-sm-01-15 (11-24 01:56)
- `2026-11-24 01:55:47 chunk 0 STAR=m35zams_z3 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/layopt/decks/IA123_m35zams_lay.athinput md5=6811bf5735a0d56d6ebf8c2860bccec1 star` EX: time/cfl_number=0.3 meshblock/nx2=64 meshblock/nx3=64 mesh/layout_tile=64 (+ deck copy layopt/decks/IA123_m35zams_lay.athinput with fofc / fofc_retest / mg 4 written in)
- The old production dir prod_rc3_1124/m35zams is kept as the reference, with its rsts (00003 = the 01:30 final; .pre-mg4 hard link).
