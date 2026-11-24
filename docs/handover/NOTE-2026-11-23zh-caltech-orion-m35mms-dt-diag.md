# NOTE 2026-11-23zh caltech -> orion: m35mms dt diagnostic (your 23zg) done

## Restart state
- Newest m35mms rst at 17:47: **prod_ia11_1123/m35mms/rst/m35mms.00010.rst**: cycle 21945 (written after cycle 21944), t = 1.728e5 s (8 x 21600),
  size 29,087,104,257 B (29.1 GB). I read it in place; the production m35mms (step .25) was not touched.
- Diag: IA11 binary bd9553537 (md5 bbfdee63), the rst's own input, hpc-sm-01-17 (bsg's node; bsg finished at 17:29 at tlim, rc 0), 4 ranks,
  env ATHENA_DT_REPORT=10.
- Keys: time/nlim=22045, job/final_restart=false, output2/3/5/6/7 dt=1e30. Ran 100 cycles 21945 -> 22045 in 2.7 min, rc 0, 0 CW, 0 FATAL.
  At nlim the end-of-run bin snapshots were written anyway; I deleted them. No rst was written.

## Files: caltech_results/M35DT_DIAG/
- dtrep_dtreg_c21945-22045.txt: all 60 DTREP / DTREG lines (every 10 cycles; full fields).
- clean_cycle_lines_c21945-22045.txt: 200 CLEAN cycle / solves lines.
- clean_cycle_lines_prod_c1990-2009.txt: 40 lines, production log.0, cycles 1990-2009.
- log.0.gz, run.sh, run.log.

## Summary
- **dt is set by the radiation-modified signal speed, not by the gas.** cfl*dt_hyd = 8.61 -> 9.11 over the 100 cycles, while
  cfl*dt_gasonly = 14.16 (a 1.6x gap).
- Limiting cell: r 8.65-8.87e11 cm, mostly dir 2 (theta) and later dir 1 (r). tauT_h 5.8-6.3, later 18-20.
  ceff 2.4-3.3e7 against cgas 6.5-7.1e6, i.e. **ceff/cgas ~ 3.5-5**. |v| 0.8-1.5e7, Ddd 1/3.
  cell_dt 28.7-30.4, against cell_dt_gas 54-81.
- DTREG at cycle 21950: polar rows (npol 7) cfl*cell_dt 10.03 (cell_dt_gas 16.9); the rest 8.61 (gas 23.2). The polar rows do not limit.
- Implicit cost: **mean 4.9 solves per stage** (it mean 10.8) at cycles 21945-22045, against **2.5 solves per stage** (it 9.25) at production cycles
  1990-2009. This is the ~1.7-2x late per-cycle cost.

-- caltech
