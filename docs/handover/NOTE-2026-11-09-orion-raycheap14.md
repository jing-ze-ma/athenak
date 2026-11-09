# NOTE 2026-11-09 ORION -> viper: raycheap14 -- implicit_lag_restage: order + cost + keys-absent gate

Answer to `TASK-2026-11-09-orion-raycheap14.md`. No code changes. Binary: `fork/raycheap2-1109` head dea3e3732 (Release,
MPI, -ffp-contract=off, same recipe as addb642d0), /orion/ptmp/jinma/raycheap14_1109/athena_dea3e3732. Decks = the
raycheap13 decks (D5 = D0 + `rad_m1/implicit_lag_restage = true`; S6 = S0 + the same key). Runs on 3 p.shared nodes
(one node each: {S0, D0/D3/D5}, {S6, 2-D gate}, {S3, cs gate}); walls are on shared nodes, so +-30% noise on the
small 2-D levels.

**Summary**
- **2-D pulse: lag_restage does NOT give 2nd order.** D5 p 1.70 / 1.22 / 0.93 (D0 1.22 / 1.12 / 1.06, D3 1.96 / 1.98 /
  1.99). The error is cut at the coarse levels (d 6.8e-4 vs D0 2.5e-3, even below D3's 8.8e-4), but the gain shrinks as dt
  goes down: at the finest pair D5 d = 4.7e-5, i.e. 3.2x D3's 1.46e-5, and the slope returns to ~1. So one re-evaluation
  at the pass-0 state removes the O(dt) lag only partly; an O(dt) remainder is left.
- **Cost: D5 = 4 Picard passes every solve (D0 3, D3 7-10).** Wall at the finest level: D5 31 s, D0 26 s (1.2x), D3
  58 s (2.2x). The target (p >= 1.8 at ~1x) is met in cost but not in order.
- **cs transport-only deck (S0 recipe): S6 and S3 are bitwise identical to S0** in all 4 data variables (m1_e, m1_f1-3)
  at c 0.15 and 0.01875 (only the PAR_DUMP header differs). d 8.979e-04 8.028e-04 5.848e-04, p 0.16 / 0.46 in all three.
  On this deck the per-stage frozen pieces do not move the answer, so neither picard/pass (S3) nor restage (S6) helps.
  S3 wall is 1.08-1.12x S0 (Picard counts identical, 13.96 -> 6.04 mean). The raycheap13 Part 3 job (S0-S5, addb642d0)
  started on p.exclusive just now; it will be reported in the raycheap13 NOTE.
- **Part 3 gate: keys absent, dea3e3732 is bitwise = addb642d0** (all bin files byte-identical) on the 2-D pulse L2 (40
  cycles, m1 + hydro_u dumps) and on cs c0.15 (nlim 40, m1 / m1_lat / m1_face dumps). Also D0 and D3 on dea3e3732
  reproduce the raycheap13 d values to all digits.

## Part 1: 2-D Cartesian diffusion pulse (hesdirk2, m1, dt 0.24 ... 0.015, 1 rank)

| arm | d (4 pairs) | p (3) | Picard mean/max/non-conv per level | wall L0-L4 (s) |
|---|---|---|---|---|
| D0 defaults (step / step) | 2.541e-03 1.094e-03 5.039e-04 2.413e-04 | 1.22 1.12 1.06 | 3.00/3/0 3.00/3/0 3.00/3/0 2.99/3/0 2.99/3/0 | 7/8/19/15/26 |
| D3 recon_lag=picard + closure_lag=pass | 8.844e-04 2.279e-04 5.789e-05 1.459e-05 | 1.96 1.98 1.99 | 9.74/12/0 9.05/11/0 8.37/10/0 7.84/10/0 6.99/9/0 | 10/14/21/43/58 |
| **D5** implicit_lag_restage=true | 6.789e-04 2.085e-04 8.928e-05 4.698e-05 | **1.70 1.22 0.93** | 4.00/4/0 4.00/4/0 4.00/4/0 3.99/4/0 3.99/4/0 | 5/7/9/24/31 |

All rc 0, FATAL 0, NOT ADMISSIBLE 0, DT_RETRY 0; all end at tlim 2.4 (cycles 10/20/40/81/161); stage/BE/fallback =
9/1/0 ... 160/1/0 in every arm. Log of D5 shows closure_lag=step (the restage is under the step lags, as designed).

## Part 2: cs Ar4 transport-only (raycheap11 T1 / raycheap13 S0 recipe), 24 ranks, tlim 20

| arm | d (3 pairs) | p (2) | Picard mean/max/non-conv per level (c 0.15 ... 0.01875) | wall (s) |
|---|---|---|---|---|
| S0 | 8.979e-04 8.028e-04 5.848e-04 | 0.16 0.46 | 13.96/24/0 10.76/19/0 7.93/16/0 6.04/12/0 | 294/435/666/1037 |
| S3 recon picard + closure pass | bitwise = S0 | 0.16 0.46 | same as S0 | 330/485/722/1120 |
| S6 S0 + implicit_lag_restage | bitwise = S0 | 0.16 0.46 | same as S0 | 318/462/688/1081 |

All rc 0, FATAL 0, NOT ADMISSIBLE 0, DT_RETRY 0, cycles 68/136/271/542. Logs show closure_lag=pass in S3, step in S0/S6.

## Part 3: keys-absent gate (dea3e3732 vs addb642d0)

| case | cycles | result |
|---|---|---|
| 2-D pulse D0_L2 | 40 | bitwise (m1.00000/1, hydro_u.00000/1) |
| cs S0_c0.15, nlim 40 | 40 | bitwise (m1, m1_lat, m1_face, 00000/1) |

Dirs: /orion/ptmp/jinma/raycheap14_1109 (gen14.py, ana14.py, runs/).
