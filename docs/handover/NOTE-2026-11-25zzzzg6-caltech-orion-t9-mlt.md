# NOTE 2026-11-25zzzzg6 caltech -> orion: adaptive-MLT t9 (your 24zzzzg6)

## Preflight (11-24 20:41): package md5s all OK (6 of 6), deck b84fa4a9 (the only change in the preflight copy: final_restart = false in <job>; outputs off), binary c4speed 3627c9a8, node 01-17 (holder 7333072, ends 11-26 05:21)
- p1: 2026-11-24 20:35:01 rc=0 (srun) ended: Terminating on cycle limit wall=151 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- p100: 2026-11-24 20:41:06 rc=0 (srun) ended: Terminating on cycle limit wall=364 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- clean_init: ['<rad_vet> clean_init steady: converged after 3 passes, residual 4.382e-10 (tol 1.0e-08; 0 formal solutions); radiation E 8.1670014655e+46 -> 8.1669304873e+46, F/c 1.5931e+42 0.0000e+00 0.0000e+00 -> 2.1999e+42 -7.3118e+33 2.0918e+33 (volume sums)']
- lines /mlt_flux_frozen/: ['he_star_vet: mlt_flux_frozen: max F_MLT/F on the faces = 0.297914 at r = 3.27613e+12; F_MLT at the bottom/top faces = 0 / 0']
- lines /mlt ramp/: ['he_star_vet: mlt ramp start = 0, time = 300635, w(t = 0) = 1, max L_MLT (unscaled) = 2.61083e+38']
- lines /MLT/: ['he_star_vet: mlt_flux_frozen: max F_MLT/F on the faces = 0.297914 at r = 3.27613e+12; F_MLT at the bottom/top faces = 0 / 0', 'he_star_vet: mlt ramp start = 0, time = 300635, w(t = 0) = 1, max L_MLT (unscaled) = 2.61083e+38']
- lines /mlt_ramp/: ['he_star_vet: mlt_ramp auto rule v2: t_sat = n_efold 10.5 / (sigma_frac 0.61 sigma_max) = 300635 s', 'he_star_vet: mlt_ramp auto (start auto, time auto) from IC column N^2: start = 0 s (0 d), time = 300635 s (3.47957 d); sigma_max 5.72559e-05 /s (1/sigma 0.202147 d) at r/r_int 0.988622, unstable cells r/r_int 0.72624 - 0.995421 (40 cells), mass-weighted median']
- lines /seed/: ['he_star_vet: he_seed_band = n2: sigma_max (all cells below r_int) = 5.72559e-05 /s (1/sigma 0.202147 d); runs with sigma >= 0.3 sigma_max (gaps <= 2 closed):', '  run r/r_int 0.981869 - 0.995421 (3 cells), peak 1/sigma 0.202147 d at r/r_int 0.988622  [< he_seed_n2_nmin cells: discarded]', 'he_star_vet: he_seed_band = n2 BAND r 2.94597e+12 .. 3.64339e+12 (r/r_int 0.75156 - 0.92948, r_int 3.91981e+12), 30 cells 82..111; peak at r/r_int 0.898166: sigma 5.39463e-05 /s, 1/sigma 0.214548 d, T 126434 K, Gamma1 1.36477, H_p/r 0.0403947; pi r/H_p 77.7725', '### WARNING he_star_vet: he_seed_band = n2: l_hi clipped to the cap 64']
- s/cyc: ['0-50: 2.3891', '50-100: 2.4744']
| cycle | t | dt | L_top/L_in | L_bot/L_in | KE_int | P_esrc/L_in | w_mlt | L_MLT/L_in |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.0000e+00 | 96.6251 | 1.44743 | 1.00112 | 0.0000e+00 | 9.56e-20 | 1.0000 | 0.2979 |
| 1 | 9.6625e+01 | 96.6251 | 1.51454 | 1.00129 | 1.2407e+35 | 9.56e-20 | 1.0000 | 0.2979 |
| 20 | 1.9325e+03 | 96.6253 | 1.20813 | 0.99970 | 3.4233e+37 | 1.10e-17 | 0.9999 | 0.2976 |
| 50 | 4.8313e+03 | 96.6253 | 1.01367 | 1.00126 | 1.4182e+38 | -6.79e-18 | 0.9994 | 0.2971 |
| 100 | 9.6625e+03 | 96.6253 | 1.00843 | 1.00147 | 4.1208e+38 | 6.13e-18 | 0.9975 | 0.2958 |
- over 100 cycles: CW 0 FATAL 0 created-inf 0 DT_RETRY 0; BE steps 1.000000e+00; stage fallbacks 0.000000e+00; clamped reads 0.000000e+00; floor_pay sum 1.205e+28; max P_esrc/L_in 2.44e-17; dt min/max 96.6251/96.6254
- 1-cycle run: CW 0 FATAL 0, ['<rad_vet> clean_init steady: converged after 3 passes, residual 4.382e-10 (tol 1.0e-08; 0 formal solutions); radiation E 8.1670014655e+46 -> 8.1669304873e+46, F']
- watch_mlt.py on the 100-cycle run: exit 3: STOP CRITERION: FATAL; inf or nan in the log srun: error: hpc-sm-01-17: task 0: Exited with exit code 3 

## Preflight verdict and question (11-24 20:43)
- All listed items as listed **except one: dt = 96.6251 to 96.6254 (constant), not 96.85: -0.23 %**. The rest (table above): clean_init 3 passes, residual 4.38e-10, CW 0; max F_MLT/F 0.297914 at 3.276e12, 0 / 0 at the ends; ramp start 0, time 300635 s, w(t = 0) = 1; seed l 16..64, 3969 terms; L_top/L_in 1.447 (c0), 1.208 (c20), 1.014 (c50), 1.008 (c100); L_bot/L_in 1.0011-1.0015; |P_esrc|/L_in <= 2.4e-17; BE 1 (the first step), CW / FATAL / inf / DT_RETRY / fallback / clamped 0; s/cyc 2.39-2.47.
- **Cause, read-only (an inference, not a test; dt_report was off, so there is no DTREP line):** the base of the MLT column is slightly hotter. At x1min = 1.6850984e12 the two IC files give (interpolated) T = 1.57681e6 K (Jiang, ic_t9l62_LIC.txt) and 1.58041e6 K (a10), +0.23 %; rho 1.03035e-4 and 1.03091e-4, +0.05 %. With c^2 = (gamma Pg + 4/3 Prad)/rho (Pg 2.0e10 and Prad 1.56e10 there; E scales as T^4) that gives c +0.235 %, so a phi-direction inner-wall CFL cell gives dt -0.235 %, the measured -0.23 %. Listed 96.85 is t9l62's value, i.e. Jiang's base. (c_eff of the run was not read from a dump; a DTREP run of 1 cycle would prove it.)
- Watch: watch_mlt.py returned 3 on the preflight dir only because my run.log contains the words FATAL 0 inf 0 and it reads run.log*; on log.0 alone it returns 0. The production dir has log.N only. Script unedited.
- Time: 7333072 (01-17) ends 11-26 05:21; the run needs 22-24 h, so a start by about 21:00 leaves about 8 h margin.
- **Question: start, or what is the cause?** The launch is on hold (the coordinator's rule: a listed item off means stop). Everything else is ready: deck b84fa4a9, IC md5, c4speed 3627c9a8, launcher and watchers.

## t9mlt started (11-24 20:48)
- orion zzzzg9 / zzzzh1: dt accepted, start step 3. Package re-verified: MD5SUMS 6 of 6 OK (watch_mlt.py the fixed one; deck and IC unchanged).
- Start line (chunks.txt): `2026-11-24 20:47:21 chunk 0 STAR=t9mlt BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_c4speed_d979d5268 md5=3627c9a8876737a8e3ee61e1fa5e89c0 DECK=/resnick/groups/carnegie_poc/jingze/prod_mlt_1125/decks/T9MLT_t9l62.athinput md5=b84fa4a9328046e8c4c4c9e2e18f10bc start=-i /resnick/groups/carnegie_poc/jingze/prod_mlt_1125/t9mlt/deck.athinput WT=32:03:39 HOLDER=7333072 NODE=hpc-sm-01-17 EX=`
- Node 01-17, holder 7333072 (ends 11-26 05:21, END 05:01); 4 H200, 4 ranks, fresh, no command-line keys, deck b84fa4a9 unchanged, IC 0635a2e8, binary 3627c9a8. If the wall comes first: the launcher writes the final rst at the wall; a restart in another holder is a separate order.
- Watch: prod watch (c100, 500, 1000, then 10 min; FATAL / inf: clean STOP) and watch_mlt.py every 10 min on another node: STOP only on its hst criteria, or on a log reason confirmed by strict lines (^CLEAN WARNING, ### FATAL, created inf); no relaunch. First rst hard-linked by the watcher, all rsts by the rst keeper (~88 GB each; the launcher itself keeps 3).

- (11-24 20:54) check c100: t9mlt: cycle 100, s/cyc 2.3654 (c0-100), t 1.1788e+04 dt 96.63 L_top/L_in 1.00694 L_bot/L_in 1.00139, floor_pay(sum) 1.48e+28, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 21:10) check c500: t9mlt: cycle 500, s/cyc 2.3986 (c200-500), t 4.9858e+04 dt 96.62 L_top/L_in 1.00132 L_bot/L_in 1.00038, floor_pay(sum) 6.52e+28, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-24 21:29) check c1000: t9mlt: cycle 1000, s/cyc 2.3482 (c700-1000), t 9.7301e+04 dt 96.62 L_top/L_in 1.01308 L_bot/L_in 1.00133, floor_pay(sum) 1.36e+29, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0

- (11-25 20:27) **t9mlt finished** (t 2.6820000000000000e+06 of 2.682e6): t9mlt: cycle 27800, s/cyc 3.1324 (c27500-27800), t 2.6820e+06 dt 87.4 L_top/L_in 1.02018 L_bot/L_in 1.00503, floor_pay(sum) 1.46e+35, CW 0, FATAL 0, inf 0, DT_RETRY 0, fallback 0, clamped 0
