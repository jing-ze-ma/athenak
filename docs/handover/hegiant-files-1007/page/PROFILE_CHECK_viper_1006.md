# Check of the user's He-giant MESA profile (10-06)

File: `profile/M2pt754_Porb100_profile29.data` (MESA r15140, written 2023-12-15; model 1255, 1816 zones).
Name and version match the Wu & Fuller 2022b / Wu & Piro 2026 grid (M_He 2.754, P_orb 100 d).

## Star (profile header)
| | value |
|---|---|
| M | 2.556 Msun (initial He-star mass 2.754, Z 0.02) |
| L (photosphere) | 1.62e4 Lsun (log 4.21) |
| R (photosphere) | 61.8 Rsun |
| Teff | 8.28 kK |
| Mdot | -5.3e-3 Msun/yr (binary mass transfer, P_orb 100 d) |
| centre | logTc 9.37, Si 0.50 / S 0.36 / Ar+Ca 0.15: Si burning, no Fe core yet (fe_core_mass 0) |
| He burning | shell at r 0.012 Rsun, m 1.53 Msun |

## Envelope
- Mass outside 1 / 10 / 30 / 50 Rsun: 0.98 / 0.72 / 0.24 / ~0.02 Msun. Far more extended mass than our own column (~0.01 Msun above 0.07 R) and Woosley's 2.7 Msun (0.12 Msun beyond 10 Rsun).
- Composition: X = 0, Y 0.93, Z 0.07 (C 0.033, O 0.016) from 1 to 30 Rsun; Z 0.055 at 50 Rsun. Our He tables are X 0, Z 0.02: C/O-enhanced opacities are needed.
- NOT in thermal equilibrium: L(r)/L_surf = 1310 at 0.05 Rsun, 610 at 1, 190 at 10, 25 at 30, ~1 above 50 Rsun. The envelope absorbs ~1e7 Lsun from below and is expanding: v_r rises outward to 6.6 km/s at the surface (R/v ~ 75 d).
- Nearly the whole envelope is convective (nabla ~0.25-0.35 vs nabla_rad 10-70); convection carries >99 % of L in the deep envelope. Mach ~0.1-0.2: (F/rho)^(1/3) 31 / 24 / 14 km/s vs c_s 226 / 130 / 67 km/s at 3 / 10 / 30 Rsun.
- Local Gamma = kappa L(r)/(4 pi c G m) ~30-170 (with the local L); near the surface (He-opacity peak, kappa up to 23 at logT 4.7-4.8) Gamma 3-10.

## ESTIMATES (mine)
- Turnover ~10 d at 10-20 Rsun; dt ~35 s for r_in 3 Rsun on a log grid of ~300 radial cells to ~150 Rsun -> ~2.5e4 steps per turnover, ~8 h per turnover on 4 Raven nodes at 256^2 (BSG cost per cell-step).
- Mass transfer is negligible over a simulation: M_env/|Mdot| ~ 200 yr.
- Time to core collapse: Si burning with no Fe core suggests weeks or less; not derivable from one profile.
