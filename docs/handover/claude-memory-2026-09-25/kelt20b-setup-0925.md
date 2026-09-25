---
name: kelt20b-setup-0925
description: 09-25 KELT-20b (MASCARA-2b) dhj setup done (/viper/ptmp2/jinma/kelt20_0925), not run; mass 2.0 M_J adopted (unmeasured, 1.5-3.38); open issues: 14 % of stellar flux below 0.26 um folded into band 10, albedo sources disagree
metadata:
  type: project
---
Singh+2024 system (Teff 8980 K, R* 1.60, R_p 1.8042 R_J, a 0.05533 AU, P 3.4741039 d, Teq 2328.6 K; Lund set gives 2256.6).
M_p 2.0 M_J (Chachan+2025 log g 3.1 +- 0.1; RV limit < 3.38; Fossati 3.0-3.5 model-dependent); 1.5 -> 3.38 M_J changes
grav 1347 -> 2775, prod dt 3.86 -> 1.75 s. Met 10x (same tables as WASP-121b; retrievals ~1-30x). SED: ATLAS9 8980 K
-> ckdata10/sw_flux/sw_band_flux_K20_11.txt; 14.4 % of stellar flux < 0.26 um is deposited with band-10 opacities (too
deep; upper inversion under-heated) -> a real limitation of the 11-band grid for A-star hosts. Albedo 0.32 (range
0.1-0.4; Singh+2024 0.36 optical vs Dang+2025 0.1 +- 0.1 Spitzer). Tint 592 K, RCB 1.35 bar, p_ref 1.344e-4 bar (TESS
effective radius; tau 0.56 chord gives 9.3e-5). Keys: ap = x1min 1.211677e10, grav 1725.786, omega 2.093261e-5, x1max
1.478866e10, dfloor 1e-16; ic_k20.txt 9.5e-11..1.1e4 bar. Grids: sparc nx1 78 (dt ~8.6 s), prod nx1 256 (dt ~2.92 s);
one KELT-20b rotation ~4.3-4.5x the cost of a WASP-121b rotation. CPU gates pass; IC T up to 2 % off the table above
1.7e-4 bar (unexplained). run.sub PLANET=k20 (arms under /viper/ptmp2/jinma/sparc_k20); smoke.sub shares sparc_0925/smoke
-> never two smokes at once. NOT submitted; decisions for the user in the morning.
