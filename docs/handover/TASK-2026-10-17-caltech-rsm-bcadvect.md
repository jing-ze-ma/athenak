# TASK 2026-10-17 raven -> caltech: IDORT rs_m decks, drop `implicit_bc_advect = false` (small, deck-only)

USER decision (10-17): the IDORT radiative-shock decks must use the production default `implicit_bc_advect = true`.
Raven re-baselined them; please make Caltech's copy match.

## Do
In your wavefix suite tree (`.../wavefix_1013/caltech_wavefix/runs/cart/`), for the five decks
`rs_m1.2/FIN3_n512`, `rs_m2/FIN3_n512`, `rs_m2/FIN3_n1024`, `rs_m3/FIN3_n512`, `rs_m3/FIN3_n1024` (`in.athinput`):
- delete the single line `implicit_bc_advect = false` in `<rad_m1>` (keep a `.bak_1017` copy);
- same in any generator template you use for rs_m decks.
Nothing else changes. Do not touch production runs.

## Why
With `false` (legacy), Mach 50 fails (shock drifts 29 cells, mass flux -5 %) and M2/M3 shocks creep 3-20x faster.
With `true`, all cases pass; accuracy vs the S_N transport refs is essentially unchanged.

## New baseline (Raven, b16b66db CPU, arms ref and ROL2U bitwise-equal)
| case | vs transport L1 Tg / precursor / T_rad at shock | final shock speed / v0 |
|---|---|---|
| rs_m1.2 n512 | 0.0002 / -1.4 % / -0.02 % | -1.3e-4 |
| rs_m2 n512 | 0.0006 / -0.8 % / -0.29 % | -9.3e-5 |
| rs_m2 n1024 | 0.0008 / -0.9 % / -0.37 % | -2.0e-5 |
| rs_m3 n512 | 0.0019 / -0.5 % / +0.06 % | -1.7e-4 |
| rs_m3 n1024 | 0.0018 / -0.5 % / -0.02 % | -2.2e-4 |
Full tables: branch littests-1011 e8c07abb, tests_m1/idort_rs/results_1017/rebase_rs_bcadv.txt.

## Report
A short NOTE confirming the edit (and, if you rerun them, your numbers for the five cases).
