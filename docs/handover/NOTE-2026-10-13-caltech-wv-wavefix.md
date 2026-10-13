# NOTE Caltech -> raven (WV): wavefix-1013 IDORT suite — DONE (10-13 15:45 PDT)

## Binaries (plain CPU, same flags as every Caltech CPU build: gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI on,
Release, PROBLEM built_in_pgens, NO Athena_M1_VGD_NMAX)
- OLD 9972a56d: `athena_cpu_omp_9972a56d` md5 **cf843e67b88d9d01606d02f4036bf13a**
- NEW wavefix-1013 @ **3972a05d** (TASK UPDATE 5c99c030): `athena_cpu_omp_3972a05d` md5 **84f8b7e1be477923f15a9813fbf9ede2**
- (41c8b210 build md5 0661dd47... superseded; its 4 BIT dirs were deleted and rerun on 3972a05d. No result uses it.)

## Run
UPDATED run_wavefix.sh step (wv_one.sh): BIT O / a (fixes 1+3) / b (+ implicit_hr_q_face=true) / z (keys off), 10 cycles,
same threads; FULL a / FULL b for their DIFF dirs. Caltech fixes: rstonce.py (one rst block), both shockref paths, gtime,
BIT verdict = FAIL if any arm rc != 0 (a FATAL run still leaves the t = 0 rst, which compares "BITWISE").
Long dirs (cart cloud / shd / rs n1024) BIT at 8 threads per arm and FULL at 32 threads; cart/ddf + cart/ddiff in the
expansion job 4456442 (4 threads). Debug jobs 4461226-30, 4468250-53 (6 nodes each).

## bit_all (359 dirs): arm (a) 207 BITWISE + 150 DIFF + 2 FAIL; arm (b) 105 BITWISE + 252 DIFF + 2 FAIL;
## arm (z) 357 BITWISE + 2 FAIL — NO (z) DIFF.
- (a) DIFF dirs (150) are all sp: rsw*/rswt*/rswv*/rswtv* (incl. _T, _cart, _r1000, _h20, _q), rlat*/rlatt* (incl. om5, _T,
  h40), marsh_r1* / marsh_r1000_u / marsh_cart, relax1(_early)_T, relax2_early_T. No cart / ss / lad / p2 dir differs in (a).
- The 2 FAIL dirs (every arm, OLD too, FULL a/b too): p2/pinc_ha20_sr30_n40/HFfZR0G2_66bdd7_64 and _d747a5_64:
  "FATAL <rad_m1>/vet_gd: too many directions" (GL 20 x 24 = 480 > default M1_VGD_NMAX 432). Rebuilding OLD+NEW with
  -D Athena_M1_VGD_NMAX=480 (or more) would run them; not done (say if wanted).

## FULL (full_lines.txt: arm dir rc time nt Picard-mean reach clamped): 406 runs (FULL_a 152, FULL_b 254)
**rc != 0 (besides the 2 p2 FAIL dirs x 2 arms): 23 FULL_a runs, all arm (a) on the rlat family, all FATAL non-finite
hydro state -> "dt = -3e-01 ... non-finite or non-positive time step" within cycles 13-137 (first "hesdirk2 stage 1 NOT
ADMISSIBLE ... redone with backward Euler", then rho = nan):**
  sp/rlat/HFfZ_1841da_{32,64,128}, sp/rlat_h40/HFfZ_1841da_64, sp/rlat_om5/HFfZ_1841da_{32,64,128},
  sp/rlat_T/HFfZ_1841da_{1,2,4,8}, sp/rlat_T/HFfZR0_72ca83_{1,2,4,8}, sp/rlat_om5_T/HFfZ_1841da_{1,2,4,8},
  sp/rlat_om5_T/HFfZR0_72ca83_{1,2,4,8}.
The same dirs run to the end with rc = 0 in arm (b) (+ implicit_hr_q_face), e.g. sp/rlat/HFfZ_1841da_32 8.3 s,
sp/rlat_T/HFfZ_1841da_1 42.8 s. rlatt* (FULL_a) all rc = 0. Every other FULL_a / FULL_b run rc = 0.

## Files (docs/handover/caltech_wavefix/)
- wavefix_runs_light.tgz (19 MB, md5 file) = the TASK tar command + `--exclude='*.bin'` (with bin: ~1 GB, mostly the
  cart/pipe dirs; over the push limit; kept at Caltech /resnick/groups/carnegie_poc/jingze/wavefix_1013/caltech_wavefix/runs).
- bit_all.txt, full_lines.txt.
