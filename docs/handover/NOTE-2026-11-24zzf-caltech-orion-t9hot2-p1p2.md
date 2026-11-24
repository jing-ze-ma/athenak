# NOTE 2026-11-24zzf caltech -> orion: T9HOT2 bisect B / F / K / W (24zzc) and T9HOT P1 / P2 (24zy, trimmed by 24zze)

- Binary athena_t9hot_e8eaa95af (md5 fbc14236..., gate PASS 3.686e-11). B/F/K/W from rst_keep/t9l62.00013.rst (cycle 5121), 300 cycles, hydro_w + vet bin every 50 cycles (single precision), stats + B sub-boxes per dump (dumps deleted after).
  P1 / P2 from t9l62's final rst (final_1124). All keys named: time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/floorcap_v=3e7 rad_vet/floorcap_rho=1e-14 rad_vet/floorcap_report=25 hydro/dt_report=10 hydro/rad_signal_speed_phase=true + the arm key.
- Stop rule: first CW / FATAL -> scancel (B/F/K/W per 24zzc). The hydro_fofc count variable is not an output variable of this deck/binary, so it is not posted.
- Files: caltech_results/T9HOT2/<arm>/ (stats.txt, dtrep.txt, hst, log.0.gz, run.log; B: box_c<cycle>.npz).

## W (11-24 10:49)
- wellbalance_dynamic=false nlim=5421 ; 2026-11-24 10:49:41 rc=0 wall=51 s 
- CW 0, FATAL 4, created inf 0, DT_RETRY 0; s/cyc n/a
- dt every 10 cycles: 
- DTREP every 50 cycles:
```
```
- hst tail (t dt L_bot L_top L_in Etot M_tot Lz_tot):
```
```
- **Arm W FATAL at startup (10:49), as designed by the pgen, not a run failure:** `he_star_vet: needs <hydro>/etotgrav = true, wellbalance_dynamic = true and wb_x1 = true (the gravity source reads the x1 well-balanced cache)` (he_star_vet.cpp:656). So hydro/wellbalance_dynamic=false cannot be bisected on this binary. **W is not rerun**; tell me if you want another WB-off form. Its node 02-17 goes to T9HOT P1 (from t9l62.00020).
