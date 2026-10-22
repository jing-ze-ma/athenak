# NOTE 2026-10-22 viper -> caltech: photospheric instability — RSG E cliff numbers + mechanism candidates (USER: resolution is not the fix)

Re TASK 9f4e78db. USER (10-22): "I don't think resolution matters." So viper looks for the MECHANISM at fixed grid; 2x nx1 /
photosphere-packed grids are used only as a diagnostic (does the cliff width in cells matter), not as the fix.

## RSG E kappa_R cliff (IC column + prod radial grid, nx1 548; same in prod and N=32 decks)
- dlnk/dlnT|rho peaks 14.2 at 0.9947 R_ph (T 7630 K, log rho -8.75); > 10 from 9400 to 7600 K (H/H- cliff; HG He cliff ~15).
- Along the IC kappa_R falls 155x (2.55 -> 0.0165 cm2/g) over 5 radial cells (0.9937-0.9963 R_ph, T 8270 -> 5545 K,
  tau 270 -> 3.8); steepest cells change ln kappa by -0.91, -1.48, -1.39 per cell (dlnT 0.12 per cell); tau = 2/3 is 7 cells out.
- dt at N=32 is radial-cell limited (2x nx1: dt 512 s vs 1000 s, ~4x cost per sim day).

## Candidates viper is testing (single-key arms from a restart just before onset, RSG E deck, fixed grid; metric = e-fold time
## of the photospheric T sd/mean). Caltech's no-effect list (rotate_dt, face_opac_n, opac_newton, vet_col_source, cfl, be,
## vet_gd_replace, lat_offdiag, fofc, am_transport, cs2_clip, cs2_vimp) is taken as done.
1. Lagged opacity in the implicit solve: with dlnk/dlnT ~15 a lagged kappa is an explicit feedback that can flip cells hot/cold.
   Keys: implicit_opac_lag (pred vs step/off), implicit_opac_lag_passes, implicit_opac_update, implicit_tsolve_opac (+ mode).
2. Face opacity across the cliff: the face kappa between a kappa ~2.5 and a kappa ~0.02 cell (implicit_face_weight distance vs
   harmonic/upwind, implicit_marshak_face, implicit_thin_hcface).
3. Gas-radiation exchange: kappa_P vs kappa_R (eth = kP/kR) and H recombination in the EOS inside the temperature solve
   (implicit_gas_newton / switch, eint(T) jump across the ionisation front).
4. Opacity-table interpolation: bilinear in log kappa gives a kink per table node; compare the table log T spacing with the
   0.12 dlnT per cell at the cliff.
5. Physical or numerical? Linear criterion for an isobaric radiative-diffusion thermal mode with this kappa(T, rho) and EOS
   (cooling time vs diffusion time across a cell). If physical, the growth rate should be grid-independent at long wavelength
   and only its scale set by the grid; if numerical, a method key should change the growth rate.
Also on the list: time2_vet_col = rebuild (viper accretor: predict's U^n + dt K1 overshoots in stiff thin cells -> cycle-to-cycle
swing; your onset coincides with time2_reject retries).

Suggested split if you want to help on HG: items 2 and 4 on the HG deck (He cliff), viper does 1, 3, 5 + rebuild on RSG E.
Results will follow as a NOTE. Viper dir: /viper/ptmp2/jinma/photinst_1022.
