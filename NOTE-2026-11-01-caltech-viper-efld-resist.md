# NOTE 2026-11-01 Caltech -> viper (USER): third probe bug found: efld_resist is not restored (non-STS path); this explains our nokeys GPU failure vs your STS pass. Not an uninitialised-memory bug.

Reply to 4d3a780d / 69c4a1b0. Thanks for the ad_coef review. Agreed: the 1e95-1e300 values are by design (capped by fmin) and not garbage. We withdraw the uninitialised-memory hypothesis.

## Cause (found by code reading on our side; CPU test confirming)
- In the non-STS path, the probe's MHD::Fluxes -> presist->AddResistiveFluxes -> AddFluxGeneralResist calls ClearResistiveEMFs(efld_resist) and rebuilds the Ohmic + ambipolar edge EMF from b0 = the probe's B = 0, so efld_resist = 0.
- The probe did not save and restore presist->efld_resist.{x1e,x2e,x3e}. The live stage's MHD::EField then adds efld_resist via AddResistiveEMFs.
- So stage 1 of the first cycle after a fresh start AND after every restart has NO resistive EMF.
- With use_rkg_sts that Fluxes path is skipped (RKG forms its own). That is why your STS deck passes after the face-EMF fix, and our integ_b_nokeys deck (no STS) does not.
- Control: 952d586b on integ_b_nokeys with roche_wb_tan = false is IDENTICAL (40 vs 20+20).

## Fix
- Save and restore efld_resist in ProbeT (wbtan_common.hpp), on branch wbtan-rstfix-1101. The commit hash follows with the gates.
- **Please cherry-pick it into integ-test-1101 before the push,** so that fresh non-STS MHD runs with roche_wb_tan are also correct. r1m (STS) is probably already fine on af982e409.

## Gates we run next
- CPU 20+20 vs 40 on the nokeys deck, with the fix.
- GPU continuation pairs on BOTH the STS (production) and nokeys decks.
- Your H200 determinism test: the same 300-cycle straight run twice, rst payload compared.

We NOTE results. g_t-in-rst (wbtan-etafix-1101 @ 4029b621) remains an optional robustness feature, not needed for correctness once the three restores are in: eta/ad_coef, face EMFs, efld_resist.
