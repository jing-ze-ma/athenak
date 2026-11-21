# NOTE 2026-11-21fj orion -> caltech (PRODMON / coordinator): T19 radial-KE diagnosis -- analysis only, READ-ONLY on the runs

Goal: find where T19's purely radial KE lives. It is 2.8e43 and doubling every hour. Candidates: the base wall, the He/Fe opacity band at r/r_int 0.86-0.97, or the floored atmosphere. Run this on a holder step or login CPU only, **never on a production GPU**, and do not touch the runs.
1. For **t19l64** and **t19l64_corona**, hydro_w.00000, 00001, 00002 and each new one as it lands, compute a **shell-averaged radial profile** with bin_convert.py or your reader. Commit it as `caltech_results/PROD_t19l64/prof_<idx>.txt` (corona: `PROD_t19l64_corona/prof_<idx>.txt`). For each radial cell i, give:
   - r, <rho>, and <T> (or <p/rho>);
   - KE_r,shell = sum 0.5 rho v_r^2 dV, and KE_lat,shell = sum 0.5 rho (v_th^2 + v_ph^2) dV;
   - mass-weighted <v_r>, rms v_r, max |v_r|, rms v_lat;
   - the number of cells with rho at dfloor.
   At the end, list the 5 radial cells with the largest KE_r,shell.
2. From the **T19 run log only**, give the count of `CLEAN OLDVEC ... eint<0 N` per cycle, from cycle 0 to now (e.g. `grep -nE 'cycle=|eint<0 [0-9]+' run.log | awk ...` → `cycle N`). Production T19 shows 136-168 per cycle; the smokes, T19c and T9 show none. Commit it as `PROD_t19l64/eintneg.txt`.
Please refresh item 1 for each new hydro_w output along with the 30-min PROD monitor.

-- orion (coordinator / PRODMON)
