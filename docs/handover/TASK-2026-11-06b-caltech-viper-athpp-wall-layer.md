# TASK 2026-11-06 (2nd) Caltech -> viper (USER): near-wall stratification in your real Athena++ 1-D BSG runs

Numbers from existing outputs in `/viper/ptmp2/jinma/athpp_ke_1106/` (arms A_af1, A_tc4); no new runs needed unless the dumps don't exist.

## Why
Our 3-D AthenaK arms with the v* mimic have an exponential, non-spherical KE jump at 2.5-4 d (e-fold 0.28 d in KE, 64x), which Linhao's 3-D Fig. 2 doesn't show.
It is a buoyancy instability of a super-adiabatic layer only 3 cells thick on the inner wall (r = 20.16-20.78 Rsun, dr 0.3125). The layer is present in our 1-D too, with the same N^2:

| | r=20.16 | 20.47 | 20.78 | 21.09 | 22.03 |
|---|---|---|---|---|---|
| IC, N^2 (1e-10 s^-2) | +9.1 | +8.9 | +8.9 | +8.8 | +9.0 |
| our 1-D 3 d, N^2 | -16.1 | -19.8 | -10.5 | +8.9 | +11.4 |
| our 3-D shell mean 3 d, N^2 | -17.6 | -17.5 | -8.8 | +7.4 | +11.5 |
| our `none` (no v*) 3 d, N^2 | +2.5 | +4.0 | +4.7 | +6.4 | +9.7 |

Our T relative to the IC at 3 d, cells 1-5: +2.8, +1.1, -0.9, -1.7, -1.5 %. Your A_af1 2nd-cell T at 1, 2, 3 d is 1.04, 0.87, 0.44e6 K, against ours 0.88e6 K (384 uniform) and 1.03e6 K (256 grid) at 3 d, so the Athena++ wall may be cold and stable instead. That's only one cell, so we need the profile.

## Please report, for A_af1 and A_tc4, at t = 0, 0.5, 1, 2, 3 d (add 5 and 10 d if available)
For each of the first 12 cells: r, rho, P_gas, E_rad, T, v_r, plus F_r if available. Raw values in double precision if possible, as a text file on this branch or a path we can fetch.
2. N^2 per cell if you compute it. Formula used here: combined gas + radiation adiabat (optically thick), N^2 = g * (1/Gamma1 * dlnP/dr - dlnrho/dr), P = P_gas + E_rad/3, with Gamma1 for the gas-radiation mixture, from shell means; first-cell gradient one-sided. The script is `ke_growth_1106/t2_strat.py` on Caltech; say so if you need it pasted.
3. How much energy the Jacobi billing loss removes per step in the wall cell versus L, at those times, if available from existing diagnostics.

## What we'd conclude
- N^2 > 0 in the first cells: his wall is stable, and the difference to chase is where the cooling lands (billing loss).
- N^2 ~ -15e-10 over at least 2 cells for days: he has the layer too, and his 3-D stays quiet for another reason (VL2, no FOFC, DO stencil, 128^2 base).

Reply as a NOTE on this branch.
