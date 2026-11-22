# NOTE 2026-11-21jh Caltech -> orion: USER sees grid-scale noise near the star in t9l62 after the binary switches; measured, it steps up twice

Metric: rms of the 2nd difference of log10 T_gas cell-to-cell along phi and along r, equatorial slice, 0.3-1.5 R_ph,IC (from the page reductions, every 3rd dump):

| segment (t9l62) | t [d] | phi | r |
|---|---|---|---|
| B_j3 (old production) | 7.5-15.5 | 0.0013-0.0028 | 0.0008-0.0013 |
| B_fdef/B_fwbs (harmonic + opac_update + fwb/fcone + flim off + conv 1e-4; cfl 0.9), restarted from 00031 at t 15.61 d | 15.63-16.9 | **0.0033 -> 0.0067** (2-3x, rising) | 0.0013-0.0020 |
| B_p1 (+ gas Newton switch 4/6, VGD_XFE=0; cfl 0.9 then 0.6), from 00035 at ~17.0 d | 17.25-19.25 | **0.018-0.030** (~10x B_j3) | **0.009-0.014** (~10x) |

The first step coincides with the k9/kp keys (15.6 d), the big one with B_p1 (17.0 d). The slices show it as cell-scale texture in T_gas and much larger log rho contrast (+-4 vs +-1) at 0.3-1.5 R_ph. CW also rose on B_p1 (3 stage fallbacks since 17:53 on testH1, dt 123 -> 45). We have no idle Caltech node for an A/B (testH1 runs t9l62 to 00:10). Suggest an A/B on Raven/viper from t9l62.00035: B_p1 with vs without clean_gas_newton_switch, and fwb/fcone vs flim, measuring this metric after ~200 cycles.

-- Caltech
