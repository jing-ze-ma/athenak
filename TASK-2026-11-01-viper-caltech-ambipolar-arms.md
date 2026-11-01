# TASK 2026-11-01 viper -> Caltech (USER): addendum to TASK-...-maxeta-sts-rkghalo: does ambipolar diffusion matter?

Same start rst, deck and duration (1.5 rot) as part A of 6e5e226a6. Compare each arm with its no-AD partner.

## Arms
- **A15 = R15 + `mhd/ambipolar = eos`.** Binary: fin51-r-fix head. With Saha x_e this is an UPPER bound on AD, because Saha overestimates the night eta_A by the same factor F (50-4000 at 1e-6..1e-2 bar).
- **T15 / TA15.** Binary: fork/tracers-defaults-1101 head (this includes the x_e tracer; use the commit with the floor-dilution fix when we NOTE it, otherwise ac98652f2 or later).
  - T15 = R15 + xe_tracer = true, AD off.
  - TA15 = T15 + ambipolar = eos.
  - This pair is the realistic test.
- **Optional:** if part A shows max_eta 1e14 is acceptable, A14 = R14 + ambipolar = eos. The AD cap max_eta_ad defaults to max_eta.

AD parameters: the defaults (ad_mion_amu 30, ad_sigv_in 1.9e-9 cm^3/s). The rho_i comes from the same x_e as eta_O: Saha on fin51, the tracer on the tracer branch.

## Report (each arm vs its partner: A15 vs R15, TA15 vs T15)
- **Cost:** s/cycle, mean dt, RKG stages/cycle. Our viper restart test gave an explicit dt 0.82 -> 0.61 s with AD; RKG barely changed.
- The same wind, T, B_rms, E_B, KE and event diagnostics as part A, at 0.5 / 1 / 1.5 rot.
- If cheap to add: the ratio eta_A/eta_O (horizontal medians, day / night / terminator) at 1e-6..1e-1 bar, from one dump.

## Our prediction (offline, viper)
- AD/advection < 1e-3 on the day side and terminators at all p > 1e-6 bar. It is ~1 only at the night top (1e-6..2e-4 bar), and only with Saha x_e.
- So we expect TA15 = T15 within the noise at p >= 1e-5 bar, and A15 to differ from R15 only at the night top.
- If that holds, ambipolar stays OFF in production (the USER decides).
