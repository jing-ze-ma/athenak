# NOTE 2026-11-24zf Caltech -> orion: RSG R_ph is continuous (rph.py on 08362-08365)

The results are in caltech_results/RSG_RC3_1124/rph_rc3.txt. Run as a CPU step on hpc-sm-02-03.
One local patch: rc3's m1 dumps name the variables vet_e / vet_f1..3, not m1_*, so the reads fall back to vet_*. Worth folding into rph.py.

Our dump times equal your continuation's (08362 is at t 7.224770e8 here too), so the comparison is at equal t:

| dump | t | R23 mean (rc3) | R23 mean (old cont.) | diff | R23 med diff | R_Teff diff | T_eff rc3 / old |
|---|---|---|---|---|---|---|---|
| 08362 | 7.22477e8 | 7.09000e13 | 7.09147e13 | -0.021 % | -0.010 % | -0.005 % | 3646.1 / 3644.9 |
| 08363 | 7.22564e8 | 7.08899e13 | 7.09062e13 | -0.023 % | -0.010 % | -0.010 % | 3645.4 / 3643.4 |
| 08364 | 7.22650e8 | 7.08762e13 | 7.08943e13 | -0.026 % | -0.010 % | -0.013 % | 3645.8 / 3642.5 |
| 08365 | 7.22736e8 | 7.08637e13 | 7.08823e13 | -0.026 % | -0.016 % | -0.014 % | 3645.5 / 3644.3 |

- R23 mean falls ~1.2-1.4e10 per dump, the same trend as the old run.
- The offset vs the continuation is a steady -0.02 to -0.03 %, inside your ~0.05 % continuity band (the 11-21 switch moved it +1.0 %).
- So R_ph is continuous. RSG keeps running.

-- Caltech
