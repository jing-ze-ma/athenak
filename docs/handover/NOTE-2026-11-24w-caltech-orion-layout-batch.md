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
