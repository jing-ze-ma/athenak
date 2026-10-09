# NOTE 2026-10-08 viper: thin-atmosphere defect of the central sp flux, and the blend flux (for Caltech, He giant N897)

From: viper. To: Caltech. Docs only. **No action without the user.**

## What we found
The default `rad_m1/implicit_flux = central` x1 face flux on the spherical-polar wedge corrupts the optically
thin atmosphere above the photosphere. We saw it in two independent setups:
- AG Car A/B 3-D shake-downs.
- The finished BSG reproduction (viper `bsg_truerepro_1008/TRUEREPRO.md`, section 6).

Your running He giant N897 (chain 4242926-29, binary f3a66907, central flux) uses the same transport. It should
have the same defect above its R_ph.

Numbers (viper gate `spblend_1008/gate2/GATE2.md`):

| quantity above R_ph | central | blend (`implicit_flux = blend`, `implicit_blend = tau_f`) |
| --- | --- | --- |
| shell L / L0 | collapses to ~0-30 % | 0.83-1.0 |
| cells or faces with inward radial flux | 36-49 % | < 1 % |

- Below the photosphere the interior is identical. In the diffuse limit the blend is bitwise central.
- Cost: +41 %.
- In the BSG the transported field above 1.05 R_ph is saturated beams (|F| ~ cE), about half of them pointing
  inward. The lateral face fluxes are ~10^2 x the radial mean. The radiation force on the atmosphere gas comes
  from these fluxes.

## Code
The blend (opt-in, default central) is **not yet in rt-integration**.
- It is on viper's release candidate rc-1009 / accel-1009, which is being gated now. It goes into rt-integration
  later in one gated step.
- Wait for viper's note announcing that rt-integration sha before rebuilding anything.
- An earlier gate on viper (merge2-1008 = sp-blend + truerepro on rt-integration, superseded) measured the
  following, with the keys off against the previous rt-integration on GPU:
  - AG Car A rst + hst: bitwise.
  - BSG restart: bitwise in run.log.

## WARNING: switching central -> blend ON A RESTART can fail
- Raven smoke 31011465 (AG Car B) restarted from a central run's rst 00010 with
  `rad_m1/implicit_flux=blend rad_m1/implicit_blend=tau_f`.
- Its first Picard solve hit 200 passes, then FATAL `Picard solve DIVERGED: resid=2.99e6 > implicit_resid_fatal=0.01`.
- The worst cells (i 426-434) are in the thin atmosphere: the field that central left there is a bad start for the
  blend solve.
- A second Raven smoke, 31011669, failed the same way: the first Picard solve was FATAL DIVERGED in the thin
  atmosphere.
- A restart from a CENTRAL-flux state records `implicit_recon_lag = picard`. Any blend restart from it must therefore
  also name `rad_m1/implicit_recon_lag=step` on the command line, because the sp blend path defaults to step only on
  fresh starts. The divergence above happens even with that key.
- On viper it also failed twice (both runs without `implicit_recon_lag=step`):
  - AG Car A central rst at cycle 10: DIVERGED resid 3.55 at the first blend solve (job 12134855; the sp-blend
    and the merged binary gave the same result).
  - The finished BSG reproduction rst 00051: DIVERGED resid 1.39e5 at the first cycle (job 12134860).
- Blend fresh starts and blend -> blend restarts are fine.
- A fresh BSG true-reproduction start with blend ran 10 cycles clean (viper apudev 12134859: Picard mean 5.4,
  max 7, 0 NON-CONVERGED).

## Suggestion (the user decides)
1. Do **NOT** switch the running He giant at a link boundary without a smoke from its own restart.
2. If the user wants to try: at a link boundary, rebuild at the rt-integration sha that viper announces (CUDA, your stack). Then smoke 10 cycles from
   the newest rst with `rad_m1/implicit_flux=blend rad_m1/implicit_blend=tau_f rad_m1/implicit_recon_lag=step` on
   the command line (allowed on restarts). Expect it to fail.
   - If it cycles: compare the interior with the central link, then continue with blend.
   - If it diverges like AG Car B: the switch needs a fix first (viper is looking at it). Keep running central
     meanwhile; the interior is not affected.
