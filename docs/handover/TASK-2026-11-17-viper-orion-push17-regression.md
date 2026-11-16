# TASK 2026-11-17 viper -> Orion (USER): CPU regression of fork `clean-cc-push17` e5c54306a vs push 16 dd468ab99

viper owns fork clean-cc (push 16 = dd468ab99) and will push 17 only after a green CPU regression from
you (plus the MI300A and A100 gates, NOTE to follow). Please run your full pre-push regression (reg/
harness, all decks, restarts) of fork `clean-cc-push17` = **e5c54306a** against push 16 **dd468ab99** and
reply as NOTE-2026-11-17-orion-viper-push17-regression.md on bsg-files-1009.

## Content (dd468ab99..e5c54306a)
- Races (fork clean-cc-race3-viper a9f481263, always on): hydro/MHD barriers where a flux sweep's face
  states change thread (a482258b5), two-stream partition barrier before the clamp-count sum
  (b5b31bffb), sp hydro barrier after the x2/x3 reconstruction (f6d3340e4), shearing box x2 loops over
  nj cells, not nj+1 (a9f481263).  No result should change where no race fires.
- Memory (clean-cc-mem2-viper 60efe679f), bitwise: cell-layout arrays sized to what the preconditioner
  reads, iw without the implicit appended blocks, no implicit-solve scratch arrays, FV pass flags and
  BiCGStab s in existing slots, packed x1-line pivots, clean_halo_vars = 16 by default (at most 16
  variables per ghost exchange of the clean solve).
  Opt-in: `clean_pc_factor_float`, `clean_fv_float` (not on a 2-D axisymmetric mesh: FATAL),
  `clean_extrap_share`, `time/dt_retry_views = clean`, `clean_line_teams`, `clean_line_chunk` (slower).
- Speed (clean-cc-speed2-viper 35d58f973): bitwise kernels of the cell-layout solve and the FV rays;
  opt-in `clean_fv_lag = off | lag | extrap` (default off) + `clean_fv_lag_tol`; with extrap the restart
  carries the inflow history in second-layer ghosts of the fv intensities (nothing new otherwise).
- yres (clean-cc-yresphys-viper 5dfc5c1a0): `clean_yres = false | true (= post) | post | local`,
  default false.
- Integration: dt_retry checkpoint also restores cln_fvis (256bfc658) and the clean_fv_lag = extrap
  inflow history (e5c54306a); clean_fv_lag works with clean_fv_float; CUDA form of the FV level kernel
  (6b5348e8e, host build unaffected); clean_fv_float restart of a rotating deck read the empty
  float array (b5a614047, CleanFvI::Of).
- **Your CPU fix (NOTE 64e442cf3, P3), ec96434ff:** on the host each x1-line team takes the contiguous
  chunk [t nline/T, (t+1) nline/T) (CleanLineRange, `if constexpr` host only; GPU order unchanged).
  Please time flz_q / cs_closed on Raven CPU vs push 15 / push 16 as in your bisect.  Details:
  /viper/ptmp2/jinma/ccpush17/INTEGRATION.md (copied below in the NOTE if you need it).

## What to run
- (A) push 17 with no keys (defaults = push 16 + clean_halo_vars 16). Expected **bitwise** to
  dd468ab99 on every deck incl. multi-rank and restart.  MI300A: 29/29 gate decks bitwise (1 and 4 ranks).
- (B) the same runs: your usual table (maxrel, CLEAN lines, WB vcs_max, dE/dP/dLz, thin1 L1, WARNING /
  DT_RETRY / FATAL).  Include one x1line deck (`rad_m1/clean_pc=x1line`).
- Opt-in, one arm each on a few fv1 decks (sp, cs, IDORT beam with clean_rays=fv1):
  `rad_m1/clean_fv_lag=extrap`; `rad_m1/clean_fv_float=true rad_m1/clean_pc_factor_float=true
  rad_m1/clean_extrap_share=true time/dt_retry_views=clean`; `rad_m1/clean_yres=local` on rhse_cs6.
  Report rc / warnings / retries and maxrel vs (A).
- Restart: (A) and the extrap arm restart bitwise to the uninterrupted run; a forced dt retry
  (dt_retry_force_cycle = 4, fac 1) with clean_fv_float on bitwise to the straight run.
- Note: command-line keys must exist in the deck (fresh runs); write them into the deck.

Reply "A bitwise / not (first deck/file); B table; opt-in table; restart bitwise / not".
