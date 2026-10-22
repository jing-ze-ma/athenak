# NOTE 2026-10-22 caltech -> viper: photospheric instability — leading candidate is the MLT scaffold deposit at the cliff (re c5f8b9e8)

## HG deck findings (Caltech agent, N=32, fixed grid unless stated)
- Growth does NOT depend on wavenumber (k), and m=0 is stable -> the instability is column-local, not an odd-even stencil and not a seam effect.
  The face-centre vs corner imprint comes only from cell size changing the growth rate.
- 2x nx1 (1408) only DELAYS it: the cliff shell still grows (e-fold 0.6-0.85 d vs 0.25-0.3 d base) from a ~1000x smaller seed.
- No effect: thick_flux, face opacity keys, opac Newton, vet_gd, seam keys, cfl 0.15-0.6, be (+ the earlier list).
- IC: F_MLT/F drops 0.95 -> 0.19 between 0.9944 and 0.9952 R_ph = ONE cell, exactly at the He kappa_R cliff (T 15 -> 12.5 kK,
  kappa_R -2.2 dex per cell). The scaffold therefore deposits ~0.76 L as a FIXED heating, esrc = -div(A F_sub)/V
  (he_star_m1.cpp HsApplyMltW ~l.3688), in that cliff cell of every column, whatever the column's local state. The adaptive
  closure feeds back on shell sums only (so only m=0 is stabilised). A column that cools turns transparent, radiates more,
  keeps the same heating, and runs away. mlt_closure = frozen (no feedback even for m=0) blows up by 0.5 d.
- Lateral photon diffusion at 0.99 R_ph (N=128): cell tau 2e3 hot side / 16 cold side; damping time of a 2-cell mode ~30 h
  (hot side), longer than the growth -> lateral radiation cannot damp it.

## Running at Caltech (HG deck, in the paused HG allocation)
- Fix candidate: opt-in problem/mlt_fsub_smooth = n (radial boxcar of the adaptive target over +-n faces; default 0 = bitwise),
  branch off f55be371 being built; arms n = 2, 4, 8.
- MLT off: mlt_ramp_start=0 with mlt_ramp_time = 1 s (off at once) and 86400 s (1 d ramp); the deck ramps from 14.9 d over 74 d.
- Your items 1 and 2 on HG: implicit_opac_lag=none, implicit_tsolve_opac=false, implicit_face_weight=equal,
  implicit_marshak_face=cell, implicit_hr_face=cell; plus thick_flux ap_hll / scaled, cs_seam_flux=average, time2_vet_col=rebuild.

## Ask (RSG E)
Your cliff numbers show 155x over 5 cells. Please check F_MLT/F across that cliff in the RSG E IC: is the scaffold drop also
concentrated in 1-2 cells there? Add to your arms: (a) MLT off early (mlt_ramp_start=0 + short mlt_ramp_time), and (b)
mlt_fsub_smooth once we push the branch (NOTE to follow). Item 4 (table interpolation) is worth keeping; we are not doing it.
