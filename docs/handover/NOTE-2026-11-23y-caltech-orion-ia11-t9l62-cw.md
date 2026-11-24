# NOTE 2026-11-23y Caltech -> orion: first CLEAN WARNING in IA11 t9l62 (as asked)

IA11 t9l62, resumed chunk (log.1, step 7333072.26, from t9l62.00011.rst). This is the only CW so far; 0 FATAL.

log.1 line 13138:
`CLEAN WARNING cycle 6338 t2s 3: stage solve failed, the step is redone with backward Euler (stage fallback 1)`

Context, the lines before and after:
```
CLEAN FWB cycle 6338 t2s 3 faces 110938539 radiation momentum sum |dp| 2.171025e+29 x1 faces sum dp -9.220004e+27
CLEAN WBFCT cycle 6338 t2s 3 rounds 6 left 0 kept 7.783269e+29
CLEAN XT cycle 6338 t2s 3 write-back closure 1.296e+01
CLEAN WARNING cycle 6338 t2s 3: stage solve failed, the step is redone with backward Euler (stage fallback 1)
CLEAN cycle 6338 t2s 1 it 2 (sum 43) pos 0 mark 0 all 0 neg 0 r1 0 r2 5878 first 0 5878 guard 0 solves 6 cwin 0
CLEAN MCAP cycle 6338 t2s 1: 5698 face momentum terms at the cap h chi_f E in the last pass
CLEAN FCONE cycle 6338 t2s 1 faces 24043399 momentum cut sum |dF| V/c^2 8.225496e+28
```
This is a stage fallback, not a floor_pay "energy unpaid". Now (~16:27): cycle 7450, t 4.828e5, dt 7.95. IA9 degraded at about t 4.76e5.
Full log: prod_ia11_1123/t9l62/log.1 on Caltech; tell us if you want it posted.

-- Caltech
