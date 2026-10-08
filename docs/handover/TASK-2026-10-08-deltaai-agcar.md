# TASK for DeltaAI: AG Car 3-D shake-down runs A and B (general EOS), smoke then run (copy of the viper set)

## 1. What
AG Carinae (LBV) envelope, he_star_m1 sp wedge (theta pi/3-2pi/3 x phi 0-pi/3), 480 x 128 x 128, implicit M1 + vet_gd,
general EOS (table, Saha H/He, X 0.36 Y 0.62), TOPS X0.36 Z0.02 opacities, MLT++ scaffold ramped off, SEED RULE.
- A = cool/violent state: log L 5.95, Teff 9 kK, R_ph 388.3 Rsun, base 130.74 Rsun (T 5e5 K), MLT++ alpha* 0.970755;
  tlim 8.064e6 s (4 Fe-peak turnovers).
- B = hot state: log L 6.17, Teff 20 kK, R_ph 101.3 Rsun, base 35.8 Rsun; tlim 4.624e5 s (4 turnovers).
Write-up on viper: /viper/ptmp2/jinma/lbv_1008/agcar/geos/IC_AGCAR_GEOS.md (user has it).

## 2. Code
Branch **he-ic-eint-from-t** on the fork, code commit **2fbc6aa1** (= rt-integration 82b7c212 + key
`problem/he_ic_eint_from_t`, IC eint from the column T with the code's own EOS; the commit adding this note is docs
only). Build with build_inc_deltaai.sh, PROBLEM=**he_star_m1**, GPU, stack as in TASK-2026-10-07-deltaai-hegiant.md
section 2 (if the He giant he_star_m1 build dir exists, an incremental rebuild at 2fbc6aa1 is enough). Record md5.
viper binary: athena_he_gpu72_2fbc6aa1_eft md5 4c66621bc838501fbd948f2fed83866a (ROCm 7.2, MI300A).

## 3. Files
`docs/handover/agcar-files-1008/` (ICs, tables, inputs with @AGCAR_DIR@, MD5SUMS). Copy the directory to work
storage (e.g. /work/nvme/bivj/jma20/agcar_1008/files), then `bash SETUP.sh <that dir>`: checks md5s and writes
agcar_shake{A,B}_ge.athinput with the paths filled in. Keep the directory in place for the whole run (the inputs and
restarts point at it).

## 4. Layout
4 MeshBlocks 480 x 64 x 64 -> **1 node x 4 GH200, 4 ranks** (viper: 1 node x 2 MI300A, 2 blocks per rank).

## 5. Smoke (60 s scale), then compare with viper
`<bin> -i agcar_shake{A,B}_ge.athinput -d <rundir> time/nlim=10`. viper reference in `agcar-files-1008/smoke_ref_viper/`
(jobs 12130897 A / 12130898 B, inputs WITH rad_m1/implicit_opac_newton_slope_max = 3 (user 10-08, opacity cliff at H/He recombination), rc 0, 0 NON-CONVERGED, 0 FATAL):

| | A | B |
|---|---|---|
| t after 10 cycles | 1.2525591814642336e+04 s | 2.4186988199012844e+03 s |
| dt at cycle 10 | 1.2525587305066963e+03 | 2.4186976857007841e+02 |
| mass (hst col 3) | 6.1192879615035610e+32 | 1.7891809043883109e+31 |
| tot-E (col 7) | 1.6552578848539149e+47 | 1.0633773260127629e+46 |
| IC T check "IC column T(rho,eint)/T_col" | 7.1e-14 | 6.4e-14 |
| "he_ic_balance cells" max | 3.24e-7 | 1.08e-7 |
| Picard mean / max | 13.5 / 16 | 6.9 / 24 |
| zone-cycles/cpu_second (2 MI300A) | 3.09e6 | 3.18e6 |

Pass: same rc / no FATAL / no NON-CONVERGED, the two IC T-check lines identical in magnitude, hst columns 1-7 within
~1e-6 relative (FMA, not bitwise; see TASK-2026-10-07 section 2), Picard within +-1-2 passes.

## 6. Runs (only if both smokes pass)
Run A and B to their tlim with restarts (every ~1 h wall) as chained jobs, 1 node x 4 GH200 each. Outputs as in the
inputs. Estimate the wall time from the smoke speed first and put it in your NOTE. If the queue allows only one
at a time, B first (shorter, 4.6e5 s).
The viper copies are NOT queued yet; the user decides whether to keep a viper backup. Do not touch other runs.

## 7. Report
NOTE-2026-10-08-deltaai-agcar.md on this branch: binary md5, smoke table against section 5, job ids, wall estimate.
