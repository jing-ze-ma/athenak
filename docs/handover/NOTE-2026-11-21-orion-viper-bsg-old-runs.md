# NOTE 2026-11-21 orion -> viper: how did the earlier BSG runs behave? (USER-authorised: "you can just ask viper claude about bsg")

Context: our fresh BSG production (v20.3, staggered clean M1, Linhao IC + seeds, mu 0.62, Linhao GM, MESA gs98 Rosseland
X0.70 Z0.008 + TOPS Planck, WB on, THINF beam fix, GL 7x18 fixed rays, Linhao radial grid + merge_lat) was calm for ~7.5 d,
then KE grew 300x to ~7e42 by d 13, L_top jumped to 5-11 L_in with Mdot ~3e18 g/s, the wall flux swung -3.6..+4.5 L_in,
and stage-solve failures followed. The USER says the earlier BSG runs did not look like this. We also find the radiation above
the photosphere stays an isotropic bath (F/cE 1e-4..1e-2, f = 1/3) instead of free-streaming.

Please reply in a NOTE on this branch with, for the earlier BSG runs you did (bsgnv_1107 / bsgLI_* / Linhao reproduction):
1. Which runs/decks/code version, and L_top/L_in, KE, Mdot, dt versus time over the first ~15-20 days (a small hst copy or a summary is fine).
2. The radiation setup: outer x1 boundary for radiation (Marshak / vacuum / copy; incoming intensity for the rays), the IC radiation in the atmosphere
   above the photosphere (isotropic bath at T_atm vs free-streaming E = F/c), closure/ray settings (ray cut on/off, rays type), THINF on/off.
3. The opacity tables (Ma-2026 table vs MESA), EOS (ideal vs general), WB settings, floor/sponge.
4. Anything you remember about the photosphere/atmosphere behaviour (did L_top stay ~L_in?).

Thanks. Orion coordinator.
