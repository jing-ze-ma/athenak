# TASK 2026-11-12 viper -> Orion: rw-wbsolve, gate implicit_cs2_wbsolve = true on the radwave suite (CPU)

## Why (viper rc18_1110/RW_P33_UNSTABLE_1112.md)
The conservative2 write-back under time2_vstage finds the gas kick with 4 fixed-point passes, gain g = om (1 + d11) E/(rho c chat)
(rad_m1_impl_wb_kernel.hpp `for (int q = 0; q < 4 && !wcf; ++q)`). After 4 passes the kick is (1 - g^4) x exact:
- your rw-drag eps = -3.06e-4 (P 3.33, g 0.133) and -1.79e-5 (P 1.667) are -g^4 (-3.16e-4, -1.98e-5);
- g > 1 (P > 25 at c 10) reverses the kick: your P 33.3 PROD FATAL / FINAL divergence (grid-scale mode, x26 per cycle, dt-independent).
The closed form `rad_m1/implicit_cs2_wbsolve = true` (455a5265, in BOTH your builds 7d233ac30 and 84cd20888) removes both on viper.
It is opt-in and not bitwise: the change is predicted O(g^4) per cell. This TASK checks that across the 18-cell suite.

## Build
Unchanged: PROD 7d233ac30 and FINAL 84cd20888 (the rw-amp1e6 builds, Release, MPI, -ffp-contract=off). No new build.

## Lists (docs/handover/scripts/rw_wbs, siblings of rw_amp1e6 as before)
gen.py wrote all_PROD / all_PRODW / all_FINAL / all_FINALW.txt, 180 cases each: 18 cells (P 0.0333/0.333/3.33 x tau 0.01/0.1/1/10/100/1000)
x N 16..256 (nt 2N) x edd (ny=4, ref=edd) and vetf (ref=vetqs_d), amp 1e-6, tol 1e-12 / lin 1e-13 / tol_rel 0 / cnorm 1e-13.
PROD = arm=PROD19; FINAL = arm=TIME + physflux_zold=fc + hr_cmem_rst; the W lists add rm_implicit_cs2_wbsolve=true (otherwise line-identical).
viper is running the subset without vetf tau 0.01 N >= 128 and tau 0.1 N 256. Those heavy cases (36 per list) are the part only you can do.
Please run ALL 4 lists anyway, so the comparison uses one machine.
```
cd <run area>/rw_wbs
export OMP_NUM_THREADS=1 RW_PKF=$PWD/../rw_amp1e6/prodkeys_all.txt
RW_EXE=<wt_rwa6_prod>/build_none/src/athena RW_RUNS=$PWD/runs_PROD   python3 ../rw_amp1e6/drive.py run all_PROD.txt   --np 30 > all_PROD.out 2>&1
RW_EXE=<wt_rwa6_prod>/build_none/src/athena RW_RUNS=$PWD/runs_PRODW  python3 ../rw_amp1e6/drive.py run all_PRODW.txt  --np 30 > all_PRODW.out 2>&1
RW_EXE=<wt_rwa6_fin>/build_none/src/athena  RW_RUNS=$PWD/runs_FINAL  python3 ../rw_amp1e6/drive.py run all_FINAL.txt  --np 30 > all_FINAL.out 2>&1
RW_EXE=<wt_rwa6_fin>/build_none/src/athena  RW_RUNS=$PWD/runs_FINALW python3 ../rw_amp1e6/drive.py run all_FINALW.txt --np 30 > all_FINALW.out 2>&1
```
Run the 4 lists concurrently if cores allow (the heavy vetf tau 0.01 N256 case is ~5000 s per core).
Pre-check: `python3 ../rw_amp1e6/drive.py input all_FINALW.txt | grep -E "^implicit_cs2_wbsolve|dt_retry =|^time2_vstage|^implicit_gas_src" | sort | uniq -c`
must show wbsolve = true, dt_retry = false, vstage true, conservative2 (and all_FINAL.txt: no wbsolve line).
A FATAL is a result: report it, do not add keys.

## Return (NOTE-2026-11-12-orion-viper-rw-wbsolve.md on this branch)
```
python3 cmp_wbs.py runs_PROD  all_PROD.txt  runs_PRODW  all_PRODW.txt  > cmp_PROD.md
python3 cmp_wbs.py runs_FINAL all_FINAL.txt runs_FINALW all_FINALW.txt > cmp_FINAL.md
```
Paste both whole: per case dev rho, dev E (max |W - B| at t = P in initial-amplitude units), |dw|/|w|, dev/g^4, plus the
per-(closure, P) worst line. Add the ok/FATAL count per list and the Picard/NC tail for any case where W and B differ by > 1e-3.
Wall per list. NOTE only, no run trees.
