# TASK viper -> caltech (USER): fin51-r-fix = 6fac5c7ff (mhd-fixes-1025); switch the Caltech runs to it

**The user wants the Caltech runs switched to this commit.** fin51-r-fix was fast-forwarded 83d97f8c9 -> 6fac5c7ff.
At each run's next restart, rebuild from 6fac5c7ff and continue. Reply with a NOTE that lists the runs switched, the
binary md5 and the first cycles (s/cycle, no FATAL).

## What changed
1. **Bug fix, hydro too (5731c97d6):** the dt_retry checkpoint did not include the `time2_vet_col = extrap2` history
   (`m1.x2_tt`, `x2_hr`, `x2_q`). A retried cycle read the failed attempt's build. It has been live since 83d97f8c9
   (extrap2 = the sph fresh-start default). It matters only on cycles where a retry fires; those were not exact before,
   and are bitwise equal to the straight run now.
   - Cost: the checkpoint grows by ~21 cell-centred fields per rank. On a GPU deck that is device memory unless
     `<time>/dt_retry_host = true`; set it if device memory is tight.
2. **MHD support** (b15d5cb5d, 20f018a56, b4ef9a1b8): `time2_posfix` (default 4 now for MHD too), `time2_start = slope`,
   `vet_gd_rotate_slope = correct`, and dt_retry all work under `<mhd>` + `<rad_m1>`. b0 is untouched by the radiation
   step (bitwise). dt_retry saves/restores u0, w0, b0 (3 faces), bcc0, the c2p caches and res.eta_b, whole Views incl.
   ghosts, so no B exchange after the restore.
   - **MHD defaults change for fresh starts:** posfix 0 -> 4, dt_retry forced off -> on (as hydro). A restart whose input
     lacks the dt_retry key keeps it off.
3. **Energy budget:** the posfix "budget gap" (-1.16e-5 vs -6.6e-6) was a measurement artifact of the trapezoid rule
   on hst rows. Measured exactly in code, residual/Etot(0) is -1.18e-9 with posfix 4 vs -1.39e-9 with posfix 0 (MHD and
   hydro). New read-only diagnostic: sph_wedge `<problem>/wg_hist_ebnd = true` writes the applied boundary radiation
   energy as user hst column 17 `Ebnd_rad`. Off by default (bitwise).

## Gates (viper; details /viper/ptmp2/jinma/mhdfix1025/RESULTS.md)
- Hydro (HG, wedge, wedge_prodkeys, slab, hsa45) and MHD decks where no retry fires: bitwise unchanged vs 83d97f8c9's
  behaviour.
- Forced retries (cut dt0 -> dt0/2; same-dt at cycles 1/4/5/7, MHD and hydro): bitwise equal to the straight runs.
- Real MHD retry on the hot-spot wedge: clean, max |div B| 6.98e-17.
- cpplint: no new findings.
- GPU (MI300A, viper apu job 12388566, 2 GPUs): the same-dt forced retries (MHD c1/c5, MHD-only poison, hydro c5) are
  bitwise equal to the straight runs; the exact budget and the real MHD retry reproduce the CPU numbers. The cut-vs-cap
  pair differed only because the cap deck hard-coded the CPU dt0/2 (GPU dt0 differs in the last digits); a rerun with the
  GPU dt0/2 is job 12388579.
