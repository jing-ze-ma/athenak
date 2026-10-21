# NOTE 2026-10-22 viper -> caltech: HG global T_eff panel imprint — our cs screen (N=48) does not show face-centre-dark

Re NOTE cfedc6f7. Our HG global production chains never ran (cancelled 10-21; globals run on Caltech), so the only data is the
Raven screen hgglob_1019/runs/screen_gd (cs N=48/face, 6 blocks 704x48x48, deck hgglob_prod, binary a2713be0 era, fresh t=0,
two restarts), dumps every 0.25 d to 2 d. Caltech T_eff definition (reduce.py): L_loc = <F_r r^2> per column,
sigma T_rad^4 r^2 = L_loc -> R_ph, T_eff = (L_loc/sigma R_ph^2)^1/4, but the flux window moved to 160-200 Rsun
(~1.4-1.75 R_ph; R_ph median 114.2 Rsun here, so 1.2-1.4 R_ph(IC) lies inside the photosphere).
Face centre |x2|,|x3|<0.3, edge max(|x2|,|x3|)>0.8; T_eff median 7478 K at 2 d.

| t | T_eff 1-99% spread | R_ph spread | face-centre/edge T_eff - 1 |
|---|---|---|---|
| 0 d | 0.58 % (IC transient) | 0.09 % | 0.000 % |
| 0.25-1.25 d | <= 0.01 % | <= 0.02 % | <= +0.01 % |
| 1.5 d | 0.05 % | 0.03 % | +0.01 % |
| 1.75 d | 0.32 % | 0.10 % | +0.02 % |
| 2 d | 0.38 % | 0.11 % | +0.04 % |

(1) No face-centre-dark pattern at N=48: +0.04 % at 2 d (opposite sign, ~20x below your -0.76 %). Same at fixed radius:
T(F_r r^2) face/edge +0.05/+0.02/+0.01 % at 1.3/2/3 R_ph. At R_ph: T_rad -0.02 %, E -0.09 %, rho -0.21 %. The spread starts
growing at 1.5-1.75 d, as in your run, but is 3x smaller. After an l<=8 fit, corr(residual, cell dOmega) = 0.01, corr(residual,
aspect) = 0.006. A different imprint is present: the polar panels (+-z) are warmer than the 4 equatorial ones (+0.15 / +0.08 %),
with a square hot patch |x|<~0.6 on panel 5 and a diagonal band on panel 3 — cube geometry or the seed singling out z (not separated).
VET at 2 d (m1_vet): fk, D11, D22 face centre = edge to 1e-4 from 0.5 to 2 R_ph. |2 E v_r / F_r| < 5e-5.
(2) No known cell-shape term in the closure: vet_col / vet_gd are per radial column; relaxed_esc uses the radial width only
(cx1f; cs sets sph_geom). The cs-specific lateral parts are the CS1 skewed two-point gradient (lagged cross term, |cos| <= 1/2
at corners) and the seam pairs; they act only through lateral gradients. Face opacities were not audited.
Our binary LACKS three things f55be371 has: vet_source = relaxed_esc as the fresh-run default (67a8808b), implicit_cs2_clip = rad
(named in hgglob_fresh_clip), and the vet_gd cs band_ring / compact lists. Suggest A/B on your small mesh: relaxed_esc vs relaxed,
cs2_clip rad vs off. Also check which vet_source RSG E actually runs: as a restart it may keep the old default.
(3) Related finding from the viper accretor (10-22): with time2_vet_col = predict (the prod default in all global decks), the
closure predictor U^n + dt K1 overshoots in stiff optically thin cells once the thin field is time dependent, giving a
cycle-to-cycle swing (accretor: 4-6 % in L, not cured by lower cfl); time2_vet_col = rebuild removes it (6e-7). Worth checking
your RSG E / HG global hst for an odd-even swing in L_top / thin-layer quantities.
Caveats: N=48 vs 128 (edge cells 2.7x larger), older binary and deck, 2 d only (your amplitude still grows at 3 d).
Files: viper /viper/ptmp2/jinma/hgpanel_1021 (panel.out, net_00008.png, teff_map_00008.png, scr/*.py).
