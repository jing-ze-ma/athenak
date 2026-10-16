# NOTE Caltech -> raven: BSG atmosphere-band diagnostic on bsgli_1014 (Linhao-IC BSG, 57.4 d), analysis only

Reply to TASK-2026-10-16-caltech-bsgatmos.md (04d6ee6a). Scripts from docs/handover/bsgatmos-1016/ (MD5SUMS OK), run unchanged.
Files: docs/handover/bsgatmos-1016/caltech_bsgli.txt, caltech_bsgli_meta.txt, caltech_bsgli.json.

## Run directory
RUNDIR = /resnick/groups/carnegie_poc/jingze/bsgli_1014/pagerun (symlinks, built by pagework_li/merge_run.sh): dumps 00000-00032 from
prod/ (1 node x 2 H200; ran to t 1.61e6 s, wide seed taken from its rst near t 1.385e6) and 00033-00116 from wide4/run (1 node x 4 H200,
restarted from the prod rst, finished t = tlim = 4.96e6 s = 57.41 d, DONE). 117 hydro_w + 117 m1 dumps, t 0.00-57.41 d. Covers the whole run.

## R_ph: override needed
The TASK rule on the FIRST m1 dump (t = 0) finds no radius: the Linhao IC sets E uniform (1.22e5) and F_r = 0 above the profile edge
(rmax 52.49 Rsun), so E/(4 F_L/c) never drops below 276 (minimum at 52.6 Rsun). The same rule on later dumps gives
R_ph = 54.02 (0.5 d), 55.49 (1 d), 54.18 (2 d), 59.53 (5 d), 54.02 (20 d), 51.99 (39.5 d), 51.52 Rsun (57.4 d).
We used the first dump where the rule resolves (dump 1, t 0.5 d, L 4.59e38): RPH_RSUN = 54.02345092838868.

## Results (caltech_bsgli.txt)
```
== caltech_bsgli: 117 dumps 0.00-57.41 d, cadence 0.500 d, R_ph 54.023 Rsun (given)
  band 0.99-1.10: max|<v_r>| 135.86 km/s (t 6.5 d), median over t>=1 d 27.10, max rms 122.92, median rms 16.02, max drho 9.47e+05, median drho 8.08
     signed band-mean v_r: min -122.55 max 37.64 km/s; fraction of dumps with mean < -1 km/s: 0.62
  band 1.01-1.10: max|<v_r>| 132.07 km/s (t 7.0 d), median over t>=1 d 26.35, max rms 118.42, median rms 16.07, max drho 9.47e+05, median drho 8.08
     signed band-mean v_r: min -118.16 max 39.03 km/s; fraction of dumps with mean < -1 km/s: 0.61
  band 0.50-0.90: max|<v_r>| 28.00 km/s (t 2.5 d), median over t>=1 d 2.02, max rms 7.14, median rms 0.72, max drho 1.17, median drho 0.426
     signed band-mean v_r: min -4.53 max 3.31 km/s; fraction of dumps with mean < -1 km/s: 0.20
  FFT: cadence 0.50 d cannot resolve ~1.3 d periods (Nyquist 1.00 d)
```
caltech_bsgli_meta.txt (dump 0):
```
pagerun: R_ph 54.02 Rsun, dr(R_ph) 0.1563 Rsun, H_P 1673.4930 Rsun, cells per H 10709.31, nx1 256, r 20.0-80.0 Rsun
```
Caveat: R_ph 54.02 lies above the IC profile edge (52.49 Rsun), where dump 0 holds floor gas (dfloor 1e-20) with the uniform E, so
the atmosphere drho (relative to dump 0) and the meta H_P / cells per H are against a floor reference, not a stellar atmosphere.
caltech_bsgli.json: per-dump t, vmax, vrms, dsw for the three bands (period_d empty).

## Run metadata
- Binary: athena_gpu_he_star_m1_8e018272_nofma (md5 4a92728d016609effe6d4586fb2e55bf), HEAD 8e0182723da2 (origin branch
  fin4-linhao-bc-8e018272 / fin4-linhao-bc); same binary for both parts.
- Deck: files/bsgLI_bc_caltech.athinput (md5 d77f9a10): Raven prod/LI/bsgLI_full_fin3_bc.athinput (md5 d37e861a) with Caltech paths;
  run.cfg KEYS empty in both parts. Grid 256^3, r 20-80 Rsun, theta pi/3-2pi/3, phi 0-pi/3; rk2, cfl 0.3, HLLC + plm + fofc,
  point-mass GM 2.6542488e27, M1 hesdirk2 + vet_gd, inner flux BC L 5.581e38, outer marshak; tlim 4.96e6 s.
- Dump cadence: hydro_w, m1, m1_vet, m1_face every 4.32e4 s (0.5 d); hst 1000 s; rst 1 d.
- IC: he_ic_linhao = true, he_ic_linhao_file = files/linhao_profile_full.txt (md5 41edc225), amp 0.05, rho_thr 1e-11, rng 0,
  gconst 6.67408e-8; he_ic_file ic_ma2026_f0.txt (md5 246a666e) for the rad part (he_ic_rad = file), he_ic_balance false,
  he_bc_inner zerocontact.
- Sponge: no he_sponge_* key in the deck -> code defaults: he_sponge_rate 0 (off), he_sponge_r0 = r_int (3.680e12 cm = 52.9 Rsun),
  mode/dmax unset (dmax off). Log: "sponge rate = 0 above r = 3.67986e+12".
- Floors: dfloor 1e-20 (dfloor_keep_velocity true), pfloor 1e-12, tfloor_kelvin 5e3, vceil 1e8; M1 e_floor 1e-4, implicit_pos_floor true.

## Wall time
atm_runs.py 105 s + meta.py 1 s on one expansion CPU node (-c 4, MaxRSS 0.95 GB), reading 117 hydro_w + 1 m1 dump (~39 GB). Plus a
first pass (96 s) that stopped on the missing R_ph, and a login-node check of the R_ph rule on 8 m1 dumps (~1 min).
