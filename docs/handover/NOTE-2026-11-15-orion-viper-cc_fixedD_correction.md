# NOTE 2026-11-15 Orion -> viper: correction to 024cd1ca5 (cc_fixedD rerun on dc67dba7b)

"0 failures among finished runs" in 024cd1ca5 meant rc = 0 only. The logs carry backward-Euler non-convergence
warnings (CLEAN WARNING ... dt_retry) in 45 staggered and 24 CELL runs (12013 warnings in total), so those rows ran with
dt_retry-modified steps and their orders are not valid:
- cell, tau 1000, time-only ladders (N256 nt 512..2048, amp 1e-4 and 1e-5, clean_tol = clean_tol_nl = 1e-15, i.e. effective
  Newton tolerance 1e-14): 280-370 warnings per run (all 12 rows); nt 4096 rows 1-2.
- cell, P 3.3 tau 10 N256 nt1024 amp 1e-6 with clean_tol 1e-16: 185.
- cell, tau 1000 space rows (tol 1e-14 / nl 1e-13): N64/N128 34-44 warnings each, N256/N512 1.
- staggered: 45 rows (tol_nl raised to 1e-12 as you asked); the 4 N4096 rows ~1230 each and still running.
Cell rows at tau 0.1 and 10 with tol 1e-14 / nl 1e-13 are clean. Suggestion: the tau 1000 ladders need a looser floor
(e.g. tol_nl >= 1e-13 effective) or a scaled tolerance; we will rerun them on our side as part of the clean-cc work.
Run dir: /orion/ptmp/jinma/ccfixedD_1115/pkg_head (logs per run).
