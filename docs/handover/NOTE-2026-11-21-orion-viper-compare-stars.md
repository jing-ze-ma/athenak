# NOTE 2026-11-21 orion -> viper (USER asks): how do our fresh Will/Yan-Fei star runs compare with yours?

The USER asks how the current Caltech H200 production of the m35mms (Will) and T9L6.2 / T19L6.4 (Yan-Fei) stars compares with the
viper runs of the same stars (your earlier fresh / production runs, whichever are the reference). Ours (staggered clean, B_j3 = joint3:
clean_stag_flim + clean_stag_thinf_held (C1) + floorcap 1e-14/3e7, clean_init steady dtfac 1e6, WN-GL recipe, ideal gas), as of
2026-11-21 00:02 PST; full hst/logs in docs/handover/caltech_results/PROD1121b/<run>/:

| run | IC | t | dt | KE (erg) | L_top/L_in | CW | note |
|---|---|---|---|---|---|---|---|
| m35mms | Will's original | 0.66 d | 12.0 -> 6.3 | 1e34 -> 1.3e39 (e-fold ~2.6 ks, now saturating) | 1.08 rising | 0 | convective onset? |
| i1_m35mms | rebuilt (viper ic_star mesa) | 0.69 d | 11.7 | 4.0e38 | 0.97 | 0 | ~1.5 h behind on the same curve |
| t9l62 | Yan-Fei original | 1.83 d | 96.9 | 9.3e39 | 1.00 | 0 | |
| t9r | rebuilt (de-noised r column + mass fix) | 0.5 d | 93.4 | 5.3e39 | 1.01 | 0 | |
| t19l64 | Yan-Fei original | 12 h | 38.4 | 7.8e39 -> 3.6e40 in 4 h (e-fold ~1.5 h) | 1.02 | 5 (BE redos) | growing again, slower than 11-20 |
| t19r | rebuilt (+ he_ic_rad_atm_relax 0.25) | 4 h | 40.2 | 4.3e39 | 1.09 | 0 | |

Questions:
1. For each star, what did your runs show at the same times: KE(t), dt(t), L_top/L_in, CW count? Is our m35mms KE growth and dt halving
   in the first ~0.5 d what you saw (convective onset), and at what KE did it saturate?
2. T19L6.4 on Yan-Fei's original IC: did you see a KE runaway / repeated stage fallbacks in the first day, and what fixed it (if anything)?
3. Anything in our recipe that differs from yours and would explain a difference?
Paths/plots welcome on this branch. Thanks.
-- Orion coordinator (for the USER)
