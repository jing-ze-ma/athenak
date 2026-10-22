# NOTE 2026-10-22 viper -> caltech: RSG E confirms the MLT-scaffold-deposit lead (re 0099ea86)

## F_MLT/F across the RSG E cliff (IC on the prod deck grid, nx1 548; active-cell index)
| cell | r / R_ph | F_MLT/F in -> out | deposit |
|---|---|---|---|
| 317 | 0.99396-0.99446 | 0.974 -> 0.948 | 0.026 L |
| 318 | 0.99446-0.99497 | 0.948 -> 0.838 | 0.110 L |
| 319 | 0.99497-0.99548 | 0.838 -> 0.155 | 0.683 L |
| 320 | 0.99548-0.99600 | 0.155 -> 0.0005 | 0.154 L |
68 % of L is deposited in ONE cell (319), 95 % in 3 cells; cell 319 is exactly the steepest kappa_R cell (dln kappa -1.48 per cell,
T 7000 K). Same picture as HG.

## Why it runs away there (discrete criterion)
A fixed deposit in a cell whose outflow F = K_face (T_i - T_j)/dr runs away only if dF/dT_i < 0, i.e. w * n * dlnT < -1
(n = dlnK/dlnT|p ~ -9 at the cliff, w the face weight). IC: -0.85 (harmonic face; just stable); 2x nx1 at IC: -0.41. In the evolved
N=32 base run the max |dlnK| between radial neighbours at the cliff is 2.8-3.1 per cell from 1 d on, i.e. past the bound.
The continuum is monotone (pure isobaric radiative-diffusion thermal modes over 0.975-1.01 R_ph are all stable, slowest decay
100-1600 d, continuum and deck grid, harmonic/arithmetic/face-T face K) -> the runaway needs the pinned deposit + an under-resolved
(steepening) cliff; refinement only delays it. Note the cliff also sits at the top of the CZ with a density inversion
(rho 1.67e-9 at 0.992 -> 2.15e-9 at 0.996 R_ph; buoyant time ~8-10 d), but the base run shows no convective signature so far
(corr(v_r, T) at the cliff shell -0.16..+0.16).

## RSG E base N=32 (fresh, f55be371) at 13.4 d: no onset yet
T_eff sd 0.06 %, panel contrast |fc/ed| < 1e-4, cliff-shell Trad sd 2-4 % (corrugation from the seed), corr(T_cliff, F_r 3 cells
up) +0.13..+0.93 (hot columns radiate more so far; the runaway predicts the opposite sign once it starts).

## Item 4 (table interpolation; viper does it on RSG E)
Table node spacing 0.025 dex in log T (dlnT 0.058) vs 0.037-0.065 dex per cliff cell (1.5-2.6 nodes per cell); the cliff is smooth
in the table (node-to-node dlnk/dlnT 12.4-14.2, slope jumps <= 1.1 per node; up to ~5 only at the cold end log T 3.70-3.75).
PCHIP vs bilinear changes dlnk/dlnT by < 0.5 in cells 314-320. Prediction: not the cause; a x4 PCHIP-resampled kappa_R/kappa_P arm
(tab4) is running anyway.

## Running on viper (RSG E, N=32)
base (to ~60 d), MLT off early (mlt_ramp_start 0 with ramp 1 s and 86400 s), tab4, time2_vet_col = rebuild. We will add
mlt_fsub_smooth n = 2/4/8 on RSG E as soon as you push the branch. Dir /viper/ptmp2/jinma/photinst_1022.
