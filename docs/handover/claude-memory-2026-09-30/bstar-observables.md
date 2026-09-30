---
name: bstar-observables
description: 2026-09-14 OBSERVATIONAL TARGETS for the 15 Msun B-star FeCZ box (log script-L 3.10, above the FeCZ threshold 2.5 but at the low-L end; no MS 10-20 Msun 3-D sim exists). SLFV (Bowman+2020 Table A.2, B0-B1 IV/V): alpha_0 6-130 umag, nu_char 0.7-4.3 /d (median ~1.2), gamma 1.2-2.5, C_w 1.5-9 umag; nu_char = 1/(2 pi tau) -> our 2.8 h turnover = 1.36 /d (on the median). METHOD (Schultz+22): light curve = bolometric F through the top boundary, fit alpha(nu) = alpha_0/(1+(nu/nu_char)^gamma) + C_w, disk correction alpha_0/sqrt(N), N = 2 pi R^2/L_box^2 = 145-850 for our box (photospheric vs FeCZ H_p) -> single-patch dL/L must be 0.7-4 % to match 50 umag. Fit both amplitude-spectrum MCMC and celerite2 SHO GP (nu_char differs up to 2x). Micro/macroturbulence (Nieva & Przybilla 12): xi 3-8, zeta_RT 4-20 km/s, LPV ~1 %; our v_MLT 2.4 km/s = Cantiello+09's detectability line. Spots (Cantiello & Braithwaite 11): dT/T 2-3e-3, size H_p, 40-150 G. Nearest sim Schultz+23 M13TAMS: v_Fe 8.3 km/s, 130 umag after /sqrt(n), nu_char 1.5 /d. RUN IMPLICATIONS: 20 turnovers = 2.3 d resolves the knee with 2-3 elements only -> extend to 40-60 turnovers; the amplitude from this run is a FLOOR (taper damps the flux-setting layer, photosphere under-resolved); check the dt-forcing floor (0.6 % v_MLT) vs the measured rms before quoting C_w.
metadata:
  type: reference
---
Refs: Bowman+19 NatAst (1905.02120), Bowman+20 A&A 640 A36 (2006.03012), Bowman & Dorn-Wallenstein 22 (2211.08347), Shen+24
(2408.11082), Cantiello+21 ApJ 915 112 (2102.05670), Lecoanet+19 (1910.01643), Anders+23 (2306.08023), Schultz+22 ApJL 924 L11
(2110.13944), Schultz+23 (2306.08034), Cantiello+09 (0903.2049), Grassitelli+15 (1507.03988), Simon-Diaz+17 (1608.05508),
Nieva & Przybilla 12 (1203.5787), Cantiello & Braithwaite 11 (1108.2030), Jiang+15 (1509.05417), Jiang+18 (1809.10187).
Extract from the box: box-integrated and per-column F_top(t) (hst every 100 s is ample cadence), T and v at tau 2/3 and
tau 0.01-1 (rms v_r, v_h; Gaussian/RT fits), dT/T map + correlation length. See [[bstar-box-growth-and-production]].
