# NOTE 2026-11-01 Caltech -> viper (USER): the integ (b) E_B 450x spike is OUR roche_wb_tan probe; fix wbtan-etafix-1101 (d0d1a738), being tested

Follow-up to our integ-gates NOTE (0e39821b): gate (b) on e32a104e and e1ecd210 fails, because E_B jumps 450x in cycle 1 when xe_tracer + ambipolar = eos + roche_wb_tan are all on.

## Cause (in the MHD WB probe, our code)
- WbTanProbe calls MHD::ConToPrim(nullptr, stage 0). With use_rkg_sts, stage 0 satisfies `(!use_rkg_sts || stage == 0)`, so ConToPrim recomputes presist->eta_b (Ohmic) and presist->ad_coef (ambipolar) from the PROBE state.
- The probe state has its passive scalars zeroed, so x_e = 0 (the xe_tracer scalar). That gives a huge eta.
- The probe saved and restored u0, w0, b0, bcc0, wder, wtemp, uflx and wbq0, but NOT the eta arrays. So the first cycle ran with the probe's eta, until TracerSource recomputed it. Hence the jump and then the decay, and why it needs all three keys.
- Otherwise the probe is side-effect-free: B and the face fields are zeroed for the probe and restored; the scalars are restored with u0/w0; no EMF or tracer source is called; the tracer init runs before the probe.

## Fix
- Branch wbtan-etafix-1101 @ d0d1a738 (off fin51-r-fix d3f878be, pushed).
- ProbeT in wbtan_common.hpp saves and restores eta_b and ad_coef around the probe, MHD only. Key off is unchanged.
- **Being tested now:** e1ecd210 + the fix (cherry-picked as a8b0dc5a), integ (b) deck, 3 cycles, E_B vs the failing INT binary. Then the gates (key off bitwise; fresh-start E_B like the non-WB arm). We NOTE the result.

## Request
- Please hold the push of e1ecd210 to fin51-r-fix until this fix is in.
- Either we merge wbtan-etafix-1101 into fin51-r-fix first and you rebase integ on top, or you cherry-pick d0d1a738 into integ-test-1101. Tell us which.
- The restart path of the integration is clean (gates a and a2 PASS).
