# NOTE 2026-11-07 Orion -> Caltech: velocity fix + tight limit + normal CFL works on the 4 x 4 column, to 5.95 d (USER's proposal; your NOTE on vel_mode)

**Result: at cfl_rad 1, error_limit 1e-9 and `vel_mode 1` (arm B, dt 419 s, 40x, not the ~20x in my status note), the column is stable to 5.95 d with no wall film. It stays close to the same fix at dt 10.6 s (arm D) to 3 d. A tight limit alone fails (A). The fix alone fails (C, NaN at 2.859 d, not a wall film). There is no cost saving, and the energy budget is not closed.** Column only, unseeded. A and B were cancelled to free a job slot, not crashed.

## The change
Before `GetTgasVel` (`rad_iteration.cpp`): `if (vel_mode == 1) ph->w.SwapAthenaArray(ph->w1);`. Stage 1 uses w^n; stage 2 uses w^{n+1/2} after its radiation kick. 3 code lines, default off; deck and cfl_rad 1 regressions identical in hst and dumps. Traced side effects: `w` is otherwise read only by the rad boundary's ComToLab (consistent) and SMR flux correction (not reached in this column, needs checking in 3-D). Diff: `velfix/velfix.diff`.

| | el1e-9 (no fix, 10.6 s) | D (fix, 10.6 s) | B (fix, 419 s) | A (no fix, ~415 s) |
|---|---|---|---|---|
| KE 0.5, 1, 2, 3 d [erg] | 1.01e40, 4.14e39, 4.26e39, 1.20e39 | 1.10e40, 4.35e39, 3.23e39, 1.03e39 | 1.02e40, 3.96e39, 3.12e39, 7.0e38 | 8.5e40, 3.02e41, 3.77e41, 5.5e40 |
| v_r 22, 25, 27 at 1 d [m/s] | -23, -89, -156 | -15, -61, -109 | -20, -57, -96 | -532, -1359, -2072 |
| dE 1, 3 d [1e44, wedge] | -0.038, -0.104 | -0.007, +0.019 | +0.018, +0.042 | -1.91, -7.53 |
| T/IC cells 0-3, 3 d | 1.004 .999 1.002 .998 | 1.003 .999 1.001 .999 | 1.003 .999 1.002 .999 | 1.053 .965 .955 .938 |
| N^2 faces 0-2, 3 d [1e-10] | +3.7, +10.9, +5.0 | +5.0, +11.0, +6.1 | +4.4, +11.3, +5.5 | -101, -8.4, -0.9 |
| Jacobi mean, max | 20, 81 | 20, 81 | 458, 966 | 527, 1209 |
| core-h per star-day | 6.2 | 9.5 | 6.0 | 7.0 |

B to 5 d: T/IC .999 in cell 1, N^2 +4.5, +11.5, +5.4; A at 5 d has cell 1 T/IC 0.846 and first-face N^2 -213.

## Your three checks
1. **L_com:** band mean at 1 d, r 20-24: `el1e-9` -0.095 (G 1.09), A -42.8 (G 40.9), D 0.97, B 1.19. So 1 - G fits, and the fix restores ~1. **A small dt-dependent residual remains in B:** cell 0 alone reads 1.65 at 1 d, against 1.00 in D.
2. **Early spike:** none scaling with dt. KE at 1000 and 3000 s: B 2.2e38, 1.06e39; `el1e-9` 2.2e38, 9.9e38. Only cell 0's velocity grows with dt, ~2x from D to B. No wall-flux patch applied.
3. **The ~1 d swing:** same period and phase in B. First extrema [d]: `el1e-9` 0.47, 1.04, 1.51, 1.97; D 0.48, 1.08, 1.51, 1.98; B 0.49, 1.09, 1.51, 2.01. Later amplitudes ~10-13 % below D. A loses the swing.

## Important: 1e-9 without the fix is not the converged answer
D differs from `el1e-9` at the same dt: v_r 30-36 % smaller, and dE at 3 d goes from -0.104 to +0.019. L_com says `el1e-9` still has G ~1 at r 20-24 at 10.6 s. **So please don't treat my `el1e-9` as the reference, as my earlier note did.** That D is the correct answer rests on the mechanism and on B and D agreeing, at two dt only. Both move further from `Lr` at 1 d (-29, -109, -199). Which is right is open.

## Not established
- **Energy:** a reflecting closed-box test didn't conserve in either mode (0.5 d, dE -2.68e44 mode 1, -3.68e44 mode 0, mostly in cell 1). That box isn't a quiet test because the IC isn't balanced against that wall, so the loss could come from the BC, the solve's non-telescoping flux or gravity. It says nothing conclusive about the real deck. B's dE > 0 is unexplained. The audit ledger hasn't been ported.
- **Cost:** no saving (458 vs 20 Jacobi per solve cancels the 40x dt). 1e-7, 1e-8 not tried with the fix; C shows 1e-6 isn't enough.
- **Lateral motion:** largest column-wide lateral velocity spread 0.67 km/s at 5 d in B, not understood.
- 3-D, seeded, SMR, beyond 5.95 d. Files: `velfix/el{A,B,C,D}/`, in the `el` columns, plus `col4x4_dt.txt`.
