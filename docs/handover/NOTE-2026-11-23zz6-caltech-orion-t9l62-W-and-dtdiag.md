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
- (11-23 21:03) status: WB [cycle=7550, t 4.85395e+05 dt 12.2, CW 28, fb 16, DT_RETRY 12, pos_floor 1]; SB [cycle=7000, t 4.79148e+05 dt 9.31, CW 3, fb 2, DT_RETRY 2, pos_floor 0];
- (11-23 21:33) status: WB [cycle=8000, t 4.89076e+05 dt 3.16, CW 48, fb 28, DT_RETRY 21, pos_floor 1]; SB [cycle=7650, t 4.86040e+05 dt 11.9, CW 3, fb 2, DT_RETRY 2, pos_floor 0]; W [cycle=6750, t 4.75700e+05 dt 4.76, CW 38, fb 20, DT_RETRY 21, pos_floor 0];
- (11-23 21:38) status:

## T0 (11-23 21:38)
- 2026-11-23 21:37:34 rc=0 last: elapsed=7.001133e+01 cycle=30100 time=1.579165e+06 dt=1.173417e+02 
- FATAL 0 
- CW 0, stage fallbacks 0, DT_RETRY 0, pos_floor created 0; s/cycle n/a; dt min 95.29 at t 1.580408e+06; skip_de per cycle: max abs 1.238e+17 at c30105; 0 cycles > 1e41
- first DTREP: `DTREP cycle=30085 time=1.57756e+06 cfl*dt_hyd=107.368 cfl*dt_gasonly=187.657 rank=2 gid=5 kji=34,3,146 dir=1 r=4.5059e+12 rho=2.89929e-11 T=1.67374e+13 pg=485.265 G1=1.66667 cgas=5.28164e+06 v=-3.96447e+07 ceff=1.32938e+08 h=3.08832e+10 tauT_h=0.595814 Erad=2.56419e+06 Ddd=0.333333 cdt_kP=5083.75 frz=-1 sc=1 cell_dt=178.947 cell_dt_gas=687.417 theta_deg=90.1758 phi=1.76408`
