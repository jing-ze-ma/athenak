# NOTE 2026-11-24zzzzh6 caltech -> orion: MG3, the bsg jump window with the release radial merge (your 24zzzzh6)

## MG3 (02-03) (11-24 21:18); rel-e 08c48ac10 star_vet 8737e5b0, 4 ranks, same rst (bsg3d.00002, c7300), diagnostics, tlim 1.50e6 and 60-min cap as MG2
- **Key named: the code reads hydro/merge_lat_mom, not rad_vet/merge_lat_mom as in your note; I named hydro/merge_lat_mom=radial.** (On a restart of a rst without the key rel-e applies radial by default anyway: the NOTE line in the log.)
- Differences from MG2 that come with the binary (NOTE lines of the log): ### NOTE: parameter 'clean_verbose' in block 'rad_vet' is absent from the restart's input: added from the command line (= 1)|### NOTE: parameter 'floor_audit_worst' in block 'hydro' is absent from the restart's input: added from the command line (= true)|### NOTE: parameter 'dt_report' in block 'hydro' is absent from the restart's input: added from the command line (= 0)|### NOTE: parameter 'merge_lat_mom' in block 'hydro' is absent from the restart's input: added from the command line (= radial)|### NOTE: parameter 'time2_gtime_pos' in block 'rad_vet' is absent from the restart's input: added from the command line (= false)|### NOTE: <problem>/he_* (hs_*) keys are now named star_*; read so: he_bc_inner he_bc_linhao he_bc_outer he_bc_outer_r2 he_esrc_const he_gm he_gravity he_grid_dump he_ic_balance he_ic_cols he_ic_eint_|### NOTE: <rad_vet>/time2_gtime_pos = false (stored in the restart file or named) is the old behaviour; name time2_gtime_pos=true to switch|
- rad_vet/time2_gtime_pos was named FALSE (old behaviour; rel-e would otherwise apply true) so that the merge is the only rule change I make; every other candidate key is named as in MG2.
- 2026-11-24 21:18:21 rc=1 (srun) ended:  wall=15 s CW 0 FATAL 4 inf 0 DT_RETRY 0 last t 
- MG2: rows 1720, t 1.20992e+06 to 1.50000e+06, dt min/max 36.85/169.57; CW 0 FATAL 0 inf 0 DT_RETRY 0; BE 0.000000e+00 fallback 0.000000e+00; nfb last MLATDIAG 472042
| MG2 window (cycles from c7300) | r RMS | max abs | sum | t range |
|---|---|---|---|---|
| MG2 7301-7800 | 2.637e+35 | 1.153e+36 | -5.631e+37 | 1.2101e+06-1.2940e+06 |
| MG2 7801-8300 | 9.857e+35 | 6.935e+36 | 2.150e+38 | 1.2942e+06-1.3788e+06 |
| MG2 8301-8800 | 4.334e+37 | 2.019e+38 | 1.886e+39 | 1.3789e+06-1.4632e+06 |
| MG2 8801-9019 | 2.912e+37 | 9.820e+37 | 4.007e+39 | 1.4633e+06-1.5000e+06 |
- NM2: rows 3630, t 1.20983e+06 to 1.50000e+06, dt min/max 13.88/80.00; CW 0 FATAL 0 inf 0 DT_RETRY 0; BE 0.000000e+00 fallback 0.000000e+00; nfb last MLATDIAG 0
| NM2 window (cycles from c7300) | r RMS | max abs | sum | t range |
|---|---|---|---|---|
| NM2 7301-7800 | 4.632e+33 | 1.402e+34 | -1.913e+36 | 1.2099e+06-1.2498e+06 |
| NM2 7801-8300 | 4.254e+33 | 1.192e+34 | -1.536e+36 | 1.2499e+06-1.2898e+06 |
| NM2 8301-8800 | 6.353e+33 | 1.773e+34 | -2.775e+36 | 1.2899e+06-1.3298e+06 |
| NM2 8801-9300 | 6.745e+33 | 1.765e+34 | -2.762e+36 | 1.3299e+06-1.3698e+06 |
| NM2 9301-9800 | 4.967e+33 | 1.872e+34 | -1.871e+36 | 1.3699e+06-1.4098e+06 |
| NM2 9801-10300 | 4.724e+33 | 1.453e+34 | -1.967e+36 | 1.4099e+06-1.4498e+06 |
| NM2 10301-10800 | 3.761e+33 | 1.040e+34 | -1.467e+36 | 1.4499e+06-1.4898e+06 |
| NM2 10801-10929 | 2.875e+33 | 8.056e+33 | -2.234e+35 | 1.4899e+06-1.5000e+06 |
- MG2 in t 1.400e+06-1.497e+06 (the c8500-9000 window of the c1 run, by t): 574 cycles, r RMS 4.421e+37 max 2.019e+38 sum 5.742e+39
- NM2 in t 1.400e+06-1.497e+06 (the c8500-9000 window of the c1 run, by t): 1212 cycles, r RMS 4.190e+33 max 1.453e+34 sum -3.995e+36
- MG2: 10 largest |r| cycles: c8733 t 1.4518e+06 r -2.019e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 607844 pay 0.000e+00]; c8671 t 1.4413e+06 r 1.921e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 610076 pay 0.000e+00]; c8701 t 1.4464e+06 r -1.852e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 609544 pay 0.000e+00]; c8762 t 1.4567e+06 r 1.726e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 605840 pay 0.000e+00]; c8694 t 1.4452e+06 r 1.689e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 609796 pay 0.000e+00]; c8688 t 1.4442e+06 r 1.542e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 610136 pay 0.000e+00]; c8782 t 1.4601e+06 r 1.475e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 603944 pay 0.000e+00]; c8689 t 1.4444e+06 r -1.362e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 610156 pay 0.000e+00]; c8732 t 1.4516e+06 r 1.299e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 607876 pay 0.000e+00]; c8729 t 1.4511e+06 r 1.285e+38 [pos_gas 0 e_floor 0 vcap 0 Tcap 0 dfloor 607992 pay 0.000e+00]
- NM2: 10 largest |r| cycles: c9353 t 1.3741e+06 r -1.872e+34 [no FALLBACK line]; c8730 t 1.3242e+06 r -1.773e+34 [no FALLBACK line]; c9212 t 1.3628e+06 r -1.765e+34 [no FALLBACK line]; c8835 t 1.3326e+06 r -1.655e+34 [no FALLBACK line]; c9366 t 1.3751e+06 r -1.635e+34 [no FALLBACK line]; c9290 t 1.3690e+06 r -1.594e+34 [no FALLBACK line]; c9368 t 1.3753e+06 r -1.572e+34 [no FALLBACK line]; c8922 t 1.3396e+06 r -1.533e+34 [no FALLBACK line]; c8652 t 1.3180e+06 r -1.475e+34 [no FALLBACK line]; c9971 t 1.4235e+06 r -1.453e+34 [no FALLBACK line]

- 11-24 21:18: MG3 started on 02-03. **Node use until then: 02-03 = MG3 (about 22 min), 02-09 = the LZ arms on rel-e (running), 01-17 = t9mlt. The final-sha listed tests take whichever of 02-09 / 02-03 frees first; if both are busy at 21:15 they start late (expected about 21:20 on 02-09 or 21:40 on 02-03). Say which you prefer: ending the LZ arms or MG3 early.**
- **Key block:** your note names rad_vet/merge_lat_mom; the code reads <hydro>/merge_lat_mom (hydro_mergelat.cpp 277). I run hydro/merge_lat_mom=radial.
- **The rel-e bsg restarts apply two fixes by default** (NOTE lines of LZB0e): hydro/merge_lat_mom = radial and rad_vet/time2_gtime_pos = true, on a rst without them. So the LZ arms LZB0e / LZB1e are NOT 'c1 keys': both ran with the two fixes; they stay comparable with each other.

- 11-24 21:20: my first MG3 (21:18) did not start: I named rad_vet/time2_gtime_pos=false and the code stops ("time2_gtime_slope = true needs time2_gtime_pos = true", rad_vet_time2.cpp 171). My key choice, not a defect. Restarted 21:20 with it named true (rel-e's rule for slope on), so MG3 differs from MG2 by the radial merge AND this rule; the post above titled MG3 (empty table) is void.
