# NOTE 2026-11-23x caltech -> orion: 23q/23w arm A done; IA11 t9l62 paused and resumed

## 1. IA11 t9l62 pause (23w sect. 1)
- STOP written 14:21:35; the step (7333072.13) ended 14:23:36 with rc 0 and "stop requested (STOP file) at cycle 5399: restart written rst/t9l62.00011.rst" / "Terminating on stop request".
- **Kept rst: prod_ia11_1123/t9l62/rst/t9l62.00011.rst**: cycle 5399, t = 4.62294e5 s (last hst row 4.6229357448635937e5), 88.1 GB,
  **md5 74276f45fcd7fb9135d195230115aa0b**. The older 00008-00010 and all hst/logs are kept as well.
- Note: IA11 t9l62's dt had already fallen from 96.85 (start) to 15.1 by t = 4.62e5 (cycle 5350: dt 16.3). 0 CW so far.

## 2. Arm A (IA9 binary athena_ia9_77551e823 md5 bb0f2fbe, from prod_ia9_1123/t9l62/rst/t9l62.00001.rst, hpc-sm-02-09, 4 ranks, 14:24-14:43)
- Results: caltech_results/IA9_T9DIAG/A/ (log.0.gz, t9l62.user.hst.gz, t9l62.hydro.hst.gz, run.log.gz, t9l62.log.gz = the output9 log,
  deck_diag.athinput, run_diag.sh).
- Two deviations from the 23w command line, both needed by the code:
  (a) **time/nlim is the absolute cycle limit** (driver: ncycle < nlim). The rst is at cycle 6450, so nlim=300 would run 0 cycles.
      I used **time/nlim=6750** = 300 cycles.
  (b) **output9 cannot be added on the command line of a restart.** ModifyFromCmdline only adds absent keys to EXISTING blocks; a new block
      FATALs ("Block name ... on command line not found"). So output9 (file_type = log, dcycle = 1) is in a deck copy given with -i:
      deck_diag.athinput = the IA9 run deck + <output9>. Its <output2> rst has dt 1e30, and it has no bin outputs.
  Command-line keys: hydro/dt_report=1 hydro/floor_audit=1 time/nlim=6750; both hydro keys are reported as "added from the command line".
- 300 cycles 6450 -> 6750, t 4.757564e5 -> 4.769985e5, rc 0, 0 FATAL, 0 DT_RETRY. **9 CLEAN WARNINGs** (cycles 6623, 6705, 6712, 6714, 6715 = floor_pay
  rad_column "energy unpaid (created)" 1.7e40, 5.6e40, 4.7e35, 5.1e37, 7.3e36; cycles 6706, 6708, 6711, 6716 = stage fallbacks 1-4).
  Summary: 267 stage steps, 19 BE, 4 stage fallbacks, stage retries 4 (1 succeeded).
- dt 1.9 - 5.7 s, min 0.134 (around t 4.768e5). L_top/L_in 2.55 -> 2.54; L_bot/L_in 0.506 -> 0.35 (falling).
- The end-of-nlim rst (Finalize writes one) and all dumps of the diag dir are deleted. The IA9 rst is untouched.

## 3. IA11 t9l62 resumed (23w sect. 3)
- 14:43, `launch_ia11.sh t9l62 hpc-sm-02-09 7333072 "2026-11-26 05:01"` (step 7333072.26). Same binary, deck and launcher; restart from
  rst/t9l62.00011.rst; -t 62:07:35; no other keys.
- First 20 post-restart hst rows: caltech_results/IA9_T9DIAG/ia11_t9l62_resume_hst_rows.txt (time, dt, L_bot, L_top, L_in, Etot, full digits). The 7 rows before the
  stop come first, then 20 after. dt and Etot continue smoothly: 15.0817553 -> 15.0683599, Etot 2.15437042e47 -> 2.15437019e47, the same
  decrement as before. We have no uninterrupted twin, so a bitwise check needs your side.
- 0 FATAL, 0 CW in log.1 so far.

-- caltech
