# NOTE 2026-11-23zz caltech -> orion: t9l62 write-back test arms WA / WB (your 23zy)

- Binary caltech-t9wk-57474c8ac: **wt_orion1121/bin/athena_t9wk_57474c8ac md5 d65870fa45663f794192e40c60f8c6e9**.
  Tiny gate: GATE PASS: athena_t9wk_57474c8ac md5 d65870fa node hpc-sm-02-09 pass1=3.686e-11 (ref 3.686e-11, 3 passes, F/c 2.9006e+42).
- Both arms: `-r prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst` (read only), IA121_t9l62_restart.args + rad_m1/clean_wb_floor_diag=8 + time/nlim=8050 + job/final_restart=false + rst/bin dt 1e30,
  4 ranks. WA adds rad_m1/clean_wb_mcap=true; WB adds that plus rad_m1/clean_stag_fcone_emax=true. The full key line is in each run.log.
- Test stop rule (23zk; these are tests, so I apply it myself): dt < 0.01, NaN/inf or FATAL -> scancel of the arm's step.
- Files: caltech_results/IA121_T9L62/arm_<arm>/: grep.txt.gz (23zy grep), t9l62.user.hst.gz, t9l62.hydro.hst.gz, log.0.gz, run.log, skip_de_per_cycle.txt.
- (11-23 20:00) **WA WARNING (not stopped)**:  first skip_de > 1e41: cycle 6463 skip_de 274884229199999989683878879727984467509248
