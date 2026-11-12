# TASK 2026-11-12 viper -> Orion: IDORT paper linear-wave ladders at amplitude 1e-6 (conv_rw_P001 / conv_rw_P10)

Same builds, runner, settings and keys as TASK-2026-11-12-orion-rw-amp1e6.md (runner `docs/handover/scripts/rw_amp1e6`:
drive.py = run_radwave.py with dt_retry forced off, athwrap.sh; PROD = 7d233ac30, FINAL = 84cd20888 = TIME keys + physflux_zold fc).
The viper amp-1e-6 grid lacks the six regimes of the IDORT paper figures (Ma+ Sect. 4.3.2), so the IDORT TIME page shows
conv_rw_P001 / conv_rw_P10 as pending. The viper apu queue is jammed, so these run ONLY on Orion (viper jobs cancelled, no duplicates).

| file (docs/handover/scripts/rw_amp1e6/) | cases | binary |
|---|---|---|
| orion_paper_FINAL.txt | P_rad/P_gas 3.33e-3 and 3.33 x tau 0.01 / 1 / 100; x1 N 16-256 and xy N 16-128; nt = 2N; amp 1e-6; tol 1e-12 / lin 1e-13 / cnorm 1e-13 / tol_rel 0; FINAL keys (physflux_zold fc, hr_cmem_rst) | 84cd20888 |
| orion_paper_PROD.txt | the same 54 cases, PROD19 keys | 7d233ac30 |

108 cases. The slow ones are the last 4 lines of each file (tau 0.01 x1 N 256 and xy N 128); start them first or on separate cores
(the grid's thin N 256 cases took ~1.5 h per case on one viper core). An N 16 smoke per arm passed on viper (both "ok").

```
cd <rw_amp1e6 dir>
export OMP_NUM_THREADS=1 RW_PKF=$PWD/prodkeys_all.txt
RW_EXE=<wt_rwa6_fin>/build_none/src/athena  RW_RUNS=$PWD/runs_paper_FINAL python3 drive.py run orion_paper_FINAL.txt --np 16 > runs_paper_FINAL.out 2>&1 &
RW_EXE=<wt_rwa6_prod>/build_none/src/athena RW_RUNS=$PWD/runs_paper_PROD  python3 drive.py run orion_paper_PROD.txt  --np 16 > runs_paper_PROD.out 2>&1 &
```
(same RW_EXE / RW_PKF convention as the main TASK.)

Pass: every case ends "ok" (a FATAL is a result; keep dt_retry off). Return in a NOTE on this branch: the results.json of each
case (small) under docs/handover/scripts/rw_amp1e6/results_paper/{FINAL,PROD}/<case dir>/results.json, plus the two .out logs.
Viper's rwconv.py then rebuilds conv_rw_P001 / conv_rw_P10 (FINAL vs PROD, x1 + xy, amplitude 1e-6) for the IDORT TIME page.
