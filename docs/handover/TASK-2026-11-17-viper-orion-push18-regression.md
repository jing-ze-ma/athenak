# TASK 2026-11-17 viper -> Orion (USER): CPU regression of fork `clean-cc-push18` b53e35565 vs push 17 5785ee51d

viper owns fork clean-cc (push 17 = 5785ee51d). Push 18 goes there only after a green CPU regression from
you, the MI300A and A100 gates, and the user's word. Please run your full pre-push regression (reg/
harness, all decks, restarts) of fork `clean-cc-push18` = **b53e35565** against push 17 **5785ee51d**.
Reply as NOTE-2026-11-17-orion-viper-push18-regression.md on bsg-files-1009.

## Content (5785ee51d..b53e35565)
- 5eb7adf3f, the restart fix: VetGdSource / VetGdCsSource fill vgd_cs over every shell when the clean
  FV rays are on (`clean && clean_rays fv*`). Before, the clean FV sweep read a stale source below the
  run-time vlat_icut (from the run's first build), so R136 restarts were not bitwise. Decks whose
  vlat_icut stays 0 are unchanged.
- 84fe1b09f: opt-in restart fingerprint diagnostic (`dbg_rst_fp` extended, `dbg_rst_fp_cycle`,
  `dbg_rst_fp_dump`). Bitwise off/on. b53e35565 is only its style fix.
- 88e5e0486 (speed3 a12b572dd): CleanFvRows, the fv1 sweep in row chunks along x1 where the cells'
  order allows. Key-free, expected bitwise.
- 1e4cfd388 (speed3 4c3b9f372): `clean_fv_lag = extrap` also keeps one pass when miss <=
  `clean_fv_lag_frac` (new key, default 0.5) x the inflow's change since the last build (OR with the old
  clean_fv_lag_tol rule). Run summary counts one-pass / swept-on builds. Opt-in path only.
- clean_yres default stays **false** (the push-17 candidate 243146188's local default is not in this branch).

## What to run
- (A) no keys: expected **bitwise** to 5785ee51d on every deck, including multi-rank, x1line, WB, RT feed decks
  and restarts. Exception: a deck whose run-time vlat_icut > 0 with clean FV rays may move through
  5eb7adf3f. If one moves, name it and the first differing file.
- (B) your usual table (maxrel, CLEAN lines, WB vcs_max, dE/dP/dLz, thin1 L1, WARNING / DT_RETRY / FATAL).
- Opt-in `rad_m1/clean_fv_lag=extrap` on a few fv1 decks (sp, cs, IDORT beam with clean_rays=fv1, rw):
  rc / warnings / retries, maxrel vs (A), the new "clean_fv_lag = extrap: N one-pass builds, M swept on"
  summary line, and dE/dP.
- Restart: (A) and the extrap arm restart bitwise to the uninterrupted run. A forced dt retry
  (dt_retry_force_cycle = 4, fac 1) with extrap is bitwise to the straight run.
- Command-line keys must exist in the deck (fresh runs), so write them into the deck.
- Please time flz_q / cs_closed on CPU vs push 17 (CleanFvRows also changes the host loop order).

Reply in this form: "A bitwise / not (first deck/file); B table; extrap table; restart bitwise / not; timing".
GPU gate records: /viper/ptmp2/jinma/ccpush18/summary.txt (viper only).
