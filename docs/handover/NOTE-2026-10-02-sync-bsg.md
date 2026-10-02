# NOTE 2026-10-02 ~05:30 (viper): SYNC for the BSG arm-2 runs on all three sites

Follows NOTE-2026-10-02-viper-update.md (05:00). Purpose: keep the three BSG arm-2 runs (viper, Caltech, DeltaAI)
identical in setup so they can be compared, and agree how we report.

## Common setup (please confirm in your next NOTE)
- Code: rt-integration **30bf6c03** (he_star_m1), no later commits for the BSG run.
- Input: the bundle bsg3d_arm2.athinput with exactly two changes: `<output1> dt = 1000.0` and
  `implicit_eos_cache` removed (commented out). Everything else byte-identical to the bundle (time_scheme be,
  cfl 0.3, tlim 4.96e6 s = 57.4 d).
- Log must show: implicit_src_stable = 1, implicit_face_weight = distance, he_bc_hse_flux = face,
  implicit_det_reduce = 1.
- Outputs as in the bundle (hst every 1000 s; bins every 0.5 d; rst cadence of the bundle). Restart links reset every
  <outputN>/last_time = floor(t/dt)*dt (viper helper bsg_1001/backfill/rst_info.py; t is 232 bytes after <par_end>).

## Viper status (05:30)
- Arm 2 and arm 1 queued as guarded jobs (one shared run dir per arm; 24 h and 6 h jobs; whichever starts first runs,
  the others exit): arm 2 bsg_a2 12059937-42, arm 1 bsg_a1 12059943-49. Not started: apu queue, our priority fell after
  tonight's many gate jobs; expected start this afternoon/evening (CEST).
- Smokes of the exact production scripts/inputs passed (dt 88.18 s, 0 FATAL/nan).

## Reporting (each site)
1. When the production job STARTS: push NOTE-2026-10-0x-<site>-bsg-started.md with job id, start time, nodes x GPUs,
   s/cycle over the first ~200 cycles, and the four log defaults above.
2. Then a short NOTE at 5, 10, 20, 30 d of simulated time (or at each link end): t, cycle, mean s/cycle, Picard mean,
   number of NON-CONVERGED solves and the largest NC residual, NEWTON-FALLBACK count, FOFC firings if fofc_report is
   on, L_top/L_in, any FATAL/STOP. Keep the hst, the bins at 5/10/20/30/40/50/57.4 d, and the final rst.
3. Comparison plan (viper does it): hst columns and the paper-figure reductions (bsg_1001/figs, analysis/ss_steady.py)
   between sites at equal t. Sites will not be bitwise (hardware), so we compare statistics after convection develops.
4. Nobody cancels another site's run; if one run finishes first, the others continue as the cross-check unless the
   user decides otherwise.

## Known behaviour to expect (not a fault)
- Convection starts in the Fe zone (~41 Rsun) around 3-4 d; v_r grows ~1.5x per 0.5 d, reaching v_MLT (~13 km/s) near
  6-7 d; the MLT scaffold ramps down 3.5-6.9 d.
- Before the stall fix, Picard stalled at the inner wall from ~4 d (round-off floor in the gas-eliminated row); with
  30bf6c03 (implicit_src_stable default) it should not. If NC solves at i = 3-4 reappear, report the count and resid.
- An extra hst row with Picard 0 at each wall-limit boundary (end-of-run dump) is normal.
