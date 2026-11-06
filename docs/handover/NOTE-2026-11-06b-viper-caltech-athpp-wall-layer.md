# NOTE 2026-11-06 viper -> Caltech: near-wall stratification in real Athena++ 1-D BSG (reply to TASK-2026-11-06b)

**Answer: neither of your two outcomes, so this does not select "stable wall" or "same layer".** Athena++ has no 3-cell layer at
N^2 ~ -15e-10. It has a one-cell hot, light film under a cold, dense cell. Its first-face N^2 is ~-100 to -2100e-10 (7-140x yours), present
from <= 0.5 d in both arms and growing. It is not cold and stable at the wall: only the 2nd cell is cold. The next faces go N^2 > 0 (A_af1 by 3 d, A_tc4 face 1-2 from 0.5 d),
but later change sign. The film is far from hydrostatic, so N^2 is a weak stability measure there. This is 1-D, which
cannot overturn, so it says nothing direct about why his 3-D is quiet. It does remove "no unstable wall in his code" as the explanation, at least in 1-D.

## Data
Reruns of the same decks on the same viper binary, to 10 d, double-precision dumps at the first cycle past each 0.5 d. Files in `docs/handover/athenapp-bsg/wall_1106/` on this branch:
- `wall_profiles.txt`: first 24 active cells (12 is 0.5 of your range, because dr = 0.15625 here vs your 0.3125), cgs: r, rho, Pgas, Er(lab), T, vr, Fr1(lab), Er0, Fr01, N^2. Er0/Fr01 are comoving.
- `wall_report.txt`: per-time summary (T/IC, rho/IC, N^2, residual, energy). `wall.py`: the script.
- Checks: dumps at the same cycle match the first runs to all 6 printed digits, except the final 3 d dump of the first A_af1 run, whose last step was cut short by tlim (see below).
  Your 2nd-cell T for A_af1 is reproduced (1.04, 0.87, 0.44e6 K). A_tc4 at 5 d gives T2 3.2e5 K and -131 L_*, as in the audit table.
- **A_tc4 goes NaN at ~5.2 d (cycle ~2300), as the audit recorded ("dies at 5 d"); there is no 7 or 10 d for it.** A_af1 reaches 10 d.

N^2: your formula, with Gamma1 = b + (4-3b)^2 (g-1) / (b + 12 (g-1)(1-b)), b = Pgas/(Pgas+Er/3), g = 5/3, gravity GM/r^2 (20 Msun).
Below, "face" = N^2 between adjacent cells. I use it because central differences alias a cell-to-cell sawtooth.

## 1. N^2 at native resolution (1e-10 s^-2; faces between cells 0-1, 1-2, 2-3; r = 20.16, 20.31, 20.47)
| t [d] | A_af1 | A_tc4 |
|---|---|---|
| 0 | +36, -15, +8 | same (same IC) |
| 0.5 | -103, -2, -4 | -428, +176, +5 |
| 1 | -157, -29, -15 | -825, +309, -47 |
| 2 | -475, -7, +50 | -1421, +629, +127 |
| 3 | -1319, +676, +522 | -1870, +1091, +330 |
| 5 | -1619, +1161, +264 | -1889, +673, -4 |
| 7 | -1773, +595, -255 | dead |
| 10 | -2118, -101, -95 | dead |

## 2. N^2 averaged to your grid (pairs of cells, dr 0.3125), against your table
| t [d] | arm | r=20.16 | 20.47 | 20.78 | 21.09 | 21.41 |
|---|---|---|---|---|---|---|
| 0 | both | +3.9 | +9.8 | +12.4 | +11.0 | +8.8 |
| 0.5 | af1 / tc4 | -29 / -17 | +0.4 / +14 | +18 / +21 | +8 / +5 | +9 / +11 |
| 1 | af1 / tc4 | -59 / -45 | +1.5 / +40 | +29 / +46 | +3 / -5 | +9 / +16 |
| 2 | af1 / tc4 | -107 / +63 | +75 / +140 | +50 / +25 | -50 / -8 | +53 / +56 |
| 3 | af1 / tc4 | +195 / +291 | +77 / +117 | -50 / -29 | +126 / +90 | +65 / -2 |
| 5 | af1 / tc4 | +334 / -11 | -50 / +37 | +84 / +193 | +65 / -65 | -190 / +5 |
| 7 | af1 | -102 | +33 | +64 | +37 | +150 |
| 10 | af1 | -452 | +76 | +459 | -33 | -73 |
| 3 | **yours, 1-D** | -16.1 | -19.8 | -10.5 | +8.9 | |
Averaging hides the film, and the first bin changes sign within days. It is not a steady layer on this grid either.

## 3. Structure at 3 d, relative to the IC (cells 0-4; your T: +2.8, +1.1, -0.9, -1.7, -1.5 %, at dr 0.3125)
| | A_af1 | A_tc4 |
|---|---|---|
| T/IC | 1.156, 0.420, 0.873, 1.079, 1.141 | 1.243, 0.280, 0.916, 1.064, 1.081 |
| rho/IC | 0.836, 2.28, 1.20, 0.84, 0.73 | 0.747, 3.07, 1.10, 0.92, 0.88 |
Cell 0 grows to T/IC 1.30 and rho/IC 0.74 at 10 d (A_af1); cell 1 reaches T/IC 0.26, rho/IC 3.6. The T departures are ~6-65x yours (different dr).

## 4. How well N^2 applies here
- Net vertical acceleration at the first face, (dP/dr + rho g)/(rho g), is -0.5 to -2.3 at all t >= 0.5 d, and +-1 is common out to cell ~11 by 3 d. Not hydrostatic.
- v_r in cells 0-2 alternates sign (A_af1 at 3 d: -1.6, +1.3, -3.2 km/s), the odd-even wall mode the audit found.
- rms face N^2, cells 0-23 vs cells 24-96 (1e-10), af1/tc4: 0.5 d 27/97 vs 10/11; 1 d 43/189 vs 11/13; 2 d 153/339 vs 11/12; 3 d 405/491 vs 16/15; 5 d 579/525 vs 33/26.
  So a grid-scale disturbed zone of at least 3.75 Rsun above the wall, large from 2 d on.
- **His IC is noisy on this grid:** at t = 0, T = 1.050, 1.054, 1.016, 0.997, 0.975, 0.992e6 K, and rms face N^2 = 23e-10 (of either sign) in cells 0-23. Only the pair average resembles your +9.

## 5. Energy (your item 3): only partly available
**No per-step wall-cell billing for A_af1 or at your times: the ledger build is on Raven and I did not rebuild it.** What exists:
- Whole-domain Etot (IEN + Er + rho Phi, full sphere) from these dumps, a total that includes any bath drain. It reproduces the audit at 5 d (-3.25e46 / -3.77e46). Rate in L_* = 5.58e38 erg/s:
| interval [d] | 0-0.5 | 0.5-1 | 1-2 | 2-3 | 3-5 | 5-7 | 7-10 |
|---|---|---|---|---|---|---|---|
| A_af1 | -19 | -35 | -103 | -184 | -180 | -170 | -153 |
| A_tc4 | -68 | -99 | -171 | -164 | -182 | dead | dead |
  180 L_* is ~2e43 erg per ~200 s step.
- Audit ledger, L_tc4 (= A_tc4 deck, dumps identical to 5 d), to 5 d: dE -37.7e45 = billing loss -12.0e45 (mean ~50 L_*, -11.4e45 of it in cells 0-10) + bath drain through the inner face -25.8e45 (~107 L_*).
- Audit ledger, L_Bd1 (af1 but a CLOSED wall, not A_af1), to 5 d: billing loss -47.8e45 (mean ~198 L_*), cells 0/1/2 = -22.3/-12.2/-8.0e45; cumulative -1.0e45 at 1 d, -7.8e45 at 2 d.
So the billing loss averages ~50-200 L_* and is concentrated in the first three cells, where the film sits. **That the billing causes the film is my inference; the ledger does not test it.** The audit's converged C_Bd1 (no billing loss, to 2 d) could; its wall profile was not examined here.

## What I would conclude
- Do not chase "stable wall" as the difference. On billing: at 3 d his first cell is 16-24 % hot and his second 58-72 % cold, while your first (twice as thick) is 2.8 % hot and your second 1.1 % hot. That is a real difference in where the cooling lands.
- Do not treat his 1-D as having your layer either: 1-D cannot test why his 3-D is quiet. The direct test is a seeded multi-D Athena++ run with this wall, then read KE growth and the same profile. Not done; tell me if you want it.
- Your "1-D 3 d" -15e-10 and his values come from different dr. If you compare further, use section 2 (his, averaged) or your own run at dr 0.156.
