# TASK for Orion: LTE Planck (and Rosseland) means at very low density for the AG Car thin atmosphere (CPU)

## Why (user 10-09)
AG Car A's optically thin atmosphere sits at rho 1e-20..1e-15 g/cm^3, T 3000-10000 K (X 0.36, Y 0.62, Z 0.02, GS98).
Every available table stops above those densities (TOPS clamps below log rho -15.6 at 5800 K; Ferguson 2005 / AESOPUS
stop at log R = -8). Viper built tables "ext2" (TOPS above log T 4.2, Ferguson 2005 GS98 Rosseland + Planck below,
blend 4.0-4.2) with a physical low-density extension below each source's floor (Rosseland = Saha electron scattering
+ absorption x (rho/rho_e)^s; Planck = kP_e (rho/rho_e)^s, s = local slope at the floor clamped to [0,1]). The Rosseland
there is mostly electron scattering and robust (Saha vs TOPS within 6 %), but the PLANCK mean is uncertain by ~1-2 dex
(Ferguson and TOPS already disagree by 1.3 dex at 9000 K). kappa_P sets how strongly the thin gas couples to radiation.

## Files
Branch `rsg-ck-1008` is yours already; this TASK is on rt-integration. The tables and their write-up are on viper:
copy what you need from the user (or ask viper to commit them): /viper/ptmp2/jinma/lbv_1008/agcar/tables_ext/
(TABLES_EXT.md, rosseland_ext2_*, planck_ext2_*, the extension script). Viper will push them on a data branch
`agcar-opac-1009` if you ask in a NOTE.

## Task
1. Compute LTE Planck means (absorption only, no scattering) and Rosseland means (with electron scattering) on log T
   3.45-4.3 x log rho -21..-12 for the AG Car mixture: atomic + ionic bound-bound lines (Kurucz/VALD-type line lists or
   whatever you used for the RSG work), bound-free, free-free, H-, H2/molecules where relevant (they matter below ~4000 K
   even at low rho? check), Rayleigh. Saha/Boltzmann + molecular equilibrium (FastChem as in your RSG work). State line
   broadening at these densities (Doppler-dominated) and the frequency grid needed to converge the Planck mean.
2. Validate where the tables overlap (Ferguson 2005 GS98 at log R -8..-3, TOPS above 5800 K at log rho >= -14): agree to
   ~0.1-0.2 dex or explain.
3. Deliver both means on the same grid as ext2 (format of TABLES_EXT.md / the AthenaK table reader) for the low-density
   corner, plus a comparison with ext2's extension: how wrong was the extension, where.
4. NOTE-2026-10-09-orion-lowrho-planck.md on a branch you push (e.g. `agcar-opac-orion-1009`), with the table files
   (small) and plots. CPU only; no AthenaK changes.
