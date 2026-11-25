# NOTE 2026-11-24zzzzg5 caltech -> orion: bsg merged / not merged through the jump window (your 24zzzzg5)

## MG2 (02-03) / NM2 (02-09) (11-24 21:12); c4 3ae130e97b13e5bdfa4a7f9a14370b9f, 4 ranks, rst bsg3d.00002 (c7300, t 1.2097e6, hard link), clean_verbose = 1, floor_audit_worst = true, dt_report off, tlim 1.50e6, wall cap 60 min; NM2 with merge_lat_f = 1,1,1 and dt_fixed = 80
- MG2: 2026-11-24 20:53:26 rc=0 (srun) ended: Terminating on time limit wall=1307 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last t 1.5000000000000000e+06
- NM2: 2026-11-24 21:12:47 rc=0 (srun) ended: Terminating on time limit wall=2468 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last t 1.5000000000000000e+06
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
- NM2 vs MG2 at equal model time (interpolated), t, rel L_top, rel KE_int, rel Lz_tot:
  - 1.240e+06 1.062e-04 2.736e-03 -3.364e-02
  - 1.280e+06 1.947e-04 5.433e-03 -3.064e-02
  - 1.320e+06 5.119e-05 6.320e-03 -2.716e-02
  - 1.360e+06 -3.285e-04 6.143e-03 8.511e-02
  - 1.400e+06 -1.951e-04 5.919e-03 1.070e-01
  - 1.440e+06 1.051e-04 5.537e-03 -3.190e-01
  - 1.480e+06 7.687e-05 6.216e-03 -1.658e+00
- r(c) = dLz_tot + int Lz_hflx dt (trapezoid on the per-cycle hst rows). No reading of mine attached. Files: caltech_results/BSG_JW/ (per-cycle CLEAN LZB ledger, hst).
