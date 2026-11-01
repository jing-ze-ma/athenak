# NOTE 2026-11-01 Caltech -> viper (USER): max_eta / STS final, +1.5 rot: R13 FAIL; R14 within 1 % in T and B but its E_B drifts up (+0.5 -> +0.85 %); cost 0.61x; explicit loses at every cap

Final reply to TASK 6e5e226a6, part A. Interim: NOTE-...-maxeta-interim (5e1ee1c4). Same setup:
- **Binary:** fin51-r-fix 65658289, md5 c7682f07bd4c6d9563e95e3875c53a5f.
- **Start:** trc_1101/keep/off_0p50.rst (t = 0.5 rot, kept).
- **Run length:** R15, R14 and R13 ran to t = 2.0 rot (+1.5 rot). X13 stopped at +0.5 rot; X14 and X15 ran 200 cycles, for cost only.

## Cost (unchanged from the interim)
- Measured on one GPU, interleaved, 2 repeats.

| arm | s/rot | vs R15 | mean dt +0 -> +1.5 rot [s] | RKG stages/cycle |
|---|---|---|---|---|
| R15 | 5090 | 1.00 | 3.921 | 15-17 |
| R14 | 3120 | 0.61 | 3.904 | 5 |
| R13 | 2830 | 0.56 | 4.101 | 3 |
| X13 | 6240 | 1.23 | 1.382 | - |
| X14 | 60,300 | 11.8 | 0.138 | - |
| X15 | 577,000 | 113 | 0.0143 | - |

- **Explicit loses at every cap.** Its resistive dt is 1.38 s at 1e13, 3x below the hydro dt. X13 costs 2.2x R13.

## Accuracy vs R15 at +0.5 / +1.0 / +1.5 rot
Horizontal means at 1e-4..1e-1 bar.

| quantity | R14 | R13 | chaos spread (csvertex A vs C, +1 rot) |
|---|---|---|---|
| max dT/T (day / term / night) | 0.26 / 0.94 / 0.52 % | 1.08 / 1.15 / 0.88 % | 0.07 % |
| E_B(p > 1e-6) | +0.46 / +0.71 / +0.85 % | +5.2 / +7.3 / +8.3 % | 0.01 % |
| B_rms 1e-5..1e-4 bar | -0.4 / -0.2 / +0.0 % | -7.1 / -7.7 / -8.5 % | 0.0 % |
| B_rms 1e-4..1e-3 bar | +0.9 / +0.6 / +1.2 % | +6.6 / +0.9 / -0.7 % | 0.1 % |
| B_rms 1e-3..1e-2 bar | +0.2 / +0.3 / +0.4 % | +2.0 / +3.6 / +4.7 % | 0.0 % |
| B_rms 1e-6..1e-5 bar | -0.0 / +0.9 / +1.2 % | -1.2 / -0.9 / -0.7 % | 0.0 % |
| winds, max abs. diff [km/s] | 0.003 / 0.007 / 0.011 | 0.011 / 0.017 / 0.014 | 0.002 |
| winds, max rel diff (all latitudes / equatorial) | 3 / 19 / 22 % (equatorial, 1e-3 bar, where u_eq ~0.02-0.03 km/s) | 30 / 79 / 149 % (1e-3 bar, where u ~0.01-0.04 km/s) | 2.8 / 4.3 % |
| KE at +1.5 rot | -0.03 % | -0.05 % | ~0 |
| total-E drift +0 -> +1.5 rot (rel) | -2.481e-5 | -2.475e-5 | R15: -2.481e-5 |
| events dfloor / efloor / tfloor vs R15 | 0.99 / 0.92 / 0.93 | 0.93 / 0.81 / 0.85 | - |

- The relative wind differences are large only in the slow 1e-3 bar layer (mean |u| <= 0.04 km/s). At 1e-4, 1e-2 and 1e-1 bar, R14's winds are within 5 %, 1 % and 1 %.
- Joule heating and ohm_drag are not output by this deck.

## Verdict (decision rule of 6e5e226a6)
- **R13 (and X13, which equals R13): FAIL.**
  - B_rms at 1e-5..1e-4 bar falls 7-8.5 %, a systematic trend at the Fe level, growing.
  - E_B rises +5 -> +8 %.
- **R14: borderline.**
  - Passes: T within 1 % (< 2 %); B_rms per band within 1.2 %; KE and energy drift equal to R15.
  - Fails: it is not within the chaos spread. E_B drifts up steadily (+0.46, +0.71, +0.85 % at +0.5, +1, +1.5 rot), and the 1e-4..1e-3 bar B_rms sits +0.6..+1.2 % high.
  - Winds differ by <= 0.011 km/s, which is > 2 % only in the slow 1e-3 bar layer.
  - Cost is 0.61x R15.
  - If the USER accepts a ~1 % E_B / B_rms bias at the Fe level, 1e14 is the cap to take. Otherwise keep 1e15.
- **RKG vs explicit:** RKG at every cap. Explicit is never within 20 % of RKG on this state.

Run dirs: trc_1101/runs/me_R15, me_R14, me_R13, me_X13 (bins at +0.5 / +1 / +1.5 rot), and met_* (cost set).
