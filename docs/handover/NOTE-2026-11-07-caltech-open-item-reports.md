# NOTE 2026-11-07 Caltech -> viper/Raven: open-item reports pushed (reply to TASK-2026-11-07-caltech-push-open-item-reports.md)

All under `docs/handover/caltech-to-raven-1107/reports/`. Unchanged copies unless marked (new). $J = /resnick/groups/carnegie_poc/jingze. No dumps/rst/binaries; largest file 375 KB; total 1.2 MB.

1. `seamE_1106/`: AUDIT.md, REPORT.md, RESTART_PROD.txt, RESTART_PROD2.txt, runs/{r1,l1,r7,r8}/arm.info.
2. `r136open_1106/`: REPORT.md, r136open.png, an.py, wbforce.py, force.py, tcmp.py, plot.py, r136noise_1106_an1.py (loader exec'd by wbforce/force), arms/*.arm.info (B BL C M ML T0-T6);
   (new) he_wb_src_linear_rmin.patch (diff of the test source vs dc2dc20c), (new) r136a1_runB_Etot_drop_timeline.txt (grep/sed/awk excerpts of runB/prod/run.log + user.hst, commands printed per block), (new) README.md.
   Not found: saved stdout of wbforce.py/force.py (numbers only in REPORT.md).
3. `rsgdiag_1107/` REPORT.md + ana_{base,ecap2,off3,on3,prednone,tsold}.txt; `cstan_1106/REPORT.md`;
   `rsg_fblocal_floors/`: (new) README.md (arms, counter column, cycle windows; note: 2.6x is raw totals over 3647-3947 vs 3000 cycles, ~2.0-2.1x per cycle), rsg_pass1.sh / rsg_pass2.sh (= cstan_1106/real/rsg.sh, rsg2.sh), arm.info + event-counter files (rsg.log) of rsg_on1/off1 (pass 1, 1173/1174) and rsg_on3/off3 (pass 2, fb_local, 3044/3031).
4. `hg_wall_layer/`: hg_relax_1106 fig1-4 + look_ts.png, results1-5.txt, an1-4.py, lib.py, plots.py, collect.py, (new) README.md pointing at the lines (L_bot/L_in results2.txt l.1/31; odd-even mvr results3.txt WALL block).
   Not found: any keys tried on the HG wall layer (none were run; it is open/untested). HG prodG deck = ../hg_runG.cfg, ../hg_last_cmd.txt.
5. `$J/r136a1_1104/runB/cont/cont.sbatch` now names `rad_m1/time2_adm_poslim=true` (KEY var on the srun command line and in the launch.log line); backup cont.sbatch.pre_poslim; copies of both in `r136a1_runB_cont/`. Not submitted, nothing started.
   CAVEAT: held job 6264757 (r136B2) was submitted from the OLD script and sbatch snapshots it, so if released it runs WITHOUT the key on its command line (the rst params still carry poslim=true). To make the job name the key it must be cancelled and resubmitted from the fixed file: USER decision; runB/prod/STOP and the hold were not touched.
