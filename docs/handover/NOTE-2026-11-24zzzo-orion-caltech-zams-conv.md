# NOTE 2026-11-24zzzo orion -> Caltech: M35ZAMS: is the iron-opacity layer convecting? (analysis of existing data only, no run)

The USER sees plots that don't look convective. From the dumps you already have of m35zams (T2 run; its last dump near t ~ 8e4 s, plus an earlier one ~ 2e4 s if it exists), on idle CPU (srun --overlap --gres=none; stop if any arm slows):
1. Against radius (shell averages): rho, T, v_r rms, lateral v rms (v_theta, v_phi separately), Mach rms and max, superadiabaticity nabla - nabla_ad (gas + radiation as the EOS has it) or N^2, kappa_R, L_rad/L_total and the convective (enthalpy + KE) flux / L. Mark where tau (from the top) = 1, 2/3 and the iron-peak T (~1.5-2e5 K).
2. The same with the thin layer (tau < 1) excluded from the rms, so the ray-effect dipole doesn't dominate.
3. KE by radial shell against time from the hst / dumps: growing, flat or decaying, with the e-fold time.
4. A lateral (theta-phi) slice of v_r and T' at the iron peak, and say what the page's plots show (which quantity, which radius, scale limits) so we can tell "too early / too weak" from "no instability".
Post a short table and figures (txt + png) in caltech_results/M35Z_CONV/. Low priority against the batch.

-- orion (coordinator)
