# NOTE 2026-11-23zz6 caltech -> orion: t9l62 arm W (23zz4) + dt-limiter arms T0/T0d/T1/T1d/T2d/LAY/LAY0 (23zz5)

- W: binary caltech-t9wk-floort-2b1d35094 -> wt_orion1121/bin/athena_t9wkfl_2b1d35094 (md5 in the gate file caltech_results/H200CHK/gate_t9wkfl_2b1d35094.txt).
  Keys: IA121_t9l62_restart.args + clean_wb_floor_diag=8 + clean_wb_mcap + clean_stag_fcone_emax + hydro/floor_from_ekin=true, nlim 8050, from rst_keep/t9l62.00012.
- T/LAY arms: athena_m35fast_13679216a, env ATHENA_DT_REPORT=5, IA121_t9l62_restart.args + the arm key, rst/bin dt 1e30, job/final_restart=false.
  rst and nlim (absolute): T0/T0d prod_fresh_1121b/t9l62/rst_p1_seg/t9l62.00037.rst (00031 no longer exists; 00037 = chunk-5 final rst, cycle 30085, t ~1.574e6), nlim 30115;
  T1/T1d rst_keep/00010 (cycle 4514, t 4.31e5), nlim 4544; T2d rst_keep/00012 (cycle 6457), nlim 6487; LAY/LAY0 rst_keep/00010, nlim 4614.
- One arm per node (t9l62 needs ~104 GB/GPU on 4 H200, so two per node do not fit). Test stop rule 23zk inside each run script (t9q_1123/run_T.sh).
- W/SB continuation with ATHENA_DT_REPORT=10: skipped. No rst is kept, so there is no state to continue from.
- Files: caltech_results/IA121_T9L62/q_<arm>/.

## IA11 production t9l62 output write times (23zz5 "report back")
- rst (88.1 GB): `<restart> rst/t9l62.00010.rst: all (file closed) written at 92.3 s`; 00011 at 85.9 s.
  Per-part times for 00011: cell/face vars 58.4 s, hesdirk2 slope 83.3 s, w0 84.1 s, rad_signal_speed 85.5 s, face flux 85.9 s.
- The bins print no timing line. Estimated from the file mtimes of one dump (11:17, output index 12/3): m1_face (1.07 GB) closed at 26.81 s, m1_vet (1.61 GB) +1.77 s,
  m1 (1.88 GB) +2.08 s, hydro_w (1.34 GB) +1.54 s. So each bin takes ~1.5-2 s; m1_face's own start is unknown.
