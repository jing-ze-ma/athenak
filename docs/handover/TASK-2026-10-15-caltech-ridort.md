# TASK raven -> Caltech: FULL IDORT for option R and ROL2U (relative Picard stop + opacity lag) — production code fin51-r-fix @ b16b66db (cross-site check, TEST ONLY)

From Raven worker RID (ridort_1014) via the coordinator. USER (10-15): send it. ROL2U is already in production (agcE1r, hgF4w) and Raven's full suite PASSED (292 dirs x REF/R/ROL2U, all
families; worst ROL2U rkw_kr_cart_T 8.3e-2 of its own ladder step; 0 FATAL / NC). This is the cross-site copy of that check on the
exact production source. TEST ONLY: no production job, deck or restart touched. Normal priority; production GPUs keep priority.

R = `rad_m1/implicit_tol_rel = 0.1` with cap `rad_m1/implicit_tol_rel_max = 1e-4`: a Picard solve that starts from the step
predictor stops at tol_eff = max(implicit_tol, min(0.1 r_0, 1e-4)), r_0 = its first change (picard-1014 study: -20 % HG, -28 % late
E1 A100 s/cycle). Default 0 = off.
ROL2U (USER approved as a third arm) = R + `implicit_opac_lag = pred`, `implicit_opac_lag_passes = 2`, `implicit_gas_newton_u4 = true`
(opacities evaluated in passes 0-1 then frozen; gas Newton in T^4 where emission dominates; oplag-1014: ROL2U A100 S3 -30 %, P268 -60 %,
HG -18 %; order 2 on kappa(T) ladders). All keys default off. Raven checks: default-off BITWISE vs FIN5.1 21c4839d (CPU 20 IDORT dirs + A100 E1).

## 1. Binary
- Branch **`fin51-r-fix` @ `b16b66db`** (the source of your production binary md5 2df1a774) = fin51-r 1bd0c2ce + rst header limit
  5ab8c87e + dt_retry restore fixes 63919cdd / 0da548c5. Straight runs are bitwise to 1bd0c2ce, so Raven's results compare 1:1.
- ONE plain CPU build with EXACTLY the flags of your ENTHBE / FIN4-suite builds (gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI on,
  Release, built_in_pgens, no Athena_M1_VGD_NMAX). Fresh build dir (CMakeLists adds driver/dt_retry.cpp). Report md5 + HEAD.

## 2. Run: `ARMS="ref r o" bash run_ridort.sh <fin51-r CPU> <cores per node> list_a.txt list_b.txt list_rk.txt` (bundle docs/handover/caltech_ridort/)
- Decks: your wavefix tree (DECKS=/resnick/groups/carnegie_poc/jingze/wavefix_1013/caltech_wavefix/runs, SHOCK=…/shockref; edit
  if moved). Part-A decks there are byte-equal to the ENTHBE tree. Set PY (python with numpy) and BINCONV=<athenak>/vis/python.
- Per dir (one.sh): `ref` = deck AS IS (FIN5.1 defaults, dt_retry on), `r` = deck + the two R keys, `o` = r + the three OL2U keys; all to
  the END, concurrently (3 x threads cores per dir),
  same threads; finished arms (rc 0) are skipped on a resubmit. list_a = the 214 part-A dirs; list_b = 35 dirs (families of the
  20-dir CPU set: cart cb cbg marsh nltek1 hcyl ddn shd pipe; sp marsh_r1_T dif_r1_T core_L16_u pulse_r1 shd patmT_L16_s);
  list_rk = 43 kappa(T) ladder dirs (OLS rk_{k0,k2,kb,kr}_T, rk_kr_cart_T, rkw_{kr,kb}_T, rkw_kr_cart_T; 1 thread; decks in the
  bundle's decks_rk/: set DECKS2=<bundle>/decks_rk).
- Threads column: 2/4, 16 for cloud + rs n1024 (Raven); use 32 as in ENTHBE if your nodes allow. ddf / ddiff ~35 min at 4 threads.
- Expected on Raven where cb / cbg hit dt_retry: the dtretry= field may differ ref/r. rc != 0 lines: report in full.

## 3. Eval (on Caltech, then push)
- `$PY reval.py $OUT/runs 'sp/*/*' 'lad/*/*' 'cart/*/*' > $OUT/reval.txt` and `ARM=o $PY reval.py ... > $OUT/reval_o.txt` (R-vs-REF L1 per variable vs the family's own ladder step).
  Then `$PY summ.py $OUT/reval.txt > $OUT/summ.txt` (and reval_o -> summ_o.txt) (one line per family: worst D_R/D_lad, Picard REF -> R, verdict).
- Push under docs/handover/caltech_ridort/: rid_cmp_list_a.txt, rid_cmp_list_b.txt, reval.txt, summ.txt, md5.txt and a light tarball
  `cd $OUT && tar czf ridort_light.tgz --exclude='*.rst' --exclude='*.bin' runs` (+ md5; split -b 44M if > 90 MB), plus
  `tar czf ridort_bin.tgz runs/sp/rlat*/*/*/bin runs/cart/eq1/*/*/bin runs/cart/eq2/*/*/bin` (the rows evaluated from bin).
- Report: docs/handover/NOTE-2026-10-15-caltech-ridort.md: md5, job ids, counts, every rc != 0 / FATAL / NON-CONVERGED line.

## 4. Optional, only if GPUs are idle: GPU spot check
- Raven ran IDORT with ROL2U on CPU only. If an H200 is free, build b16b66db for GPU (HOPPER90, nofma recipe of builds_fin51) WITH
  built_in_pgens (production binaries are he_star-only), then on 1 H200:
  `ARMS="ref o" bash run_ridort.sh <gpu build> 8 list_b.txt list_rk.txt` with OUT=<...>/ridort_gpu (one GPU, arms share it; threads
  column only sets OMP threads). Eval the same way (ARM=o reval.py -> summ_o_gpu.txt). Verdict = same per-family D_o/D_lad pass rule.
  Skip this part if it needs production GPUs.

## Coordinator notes
- Raven's results: /raven/ptmp/jinma/ridort_1014 (summ.txt / summ_o.txt); compare family by family, expect the same verdicts
  (CPU numbers within round-off of Raven's; arm ref vs r/o differences, not bitwise to Raven).
