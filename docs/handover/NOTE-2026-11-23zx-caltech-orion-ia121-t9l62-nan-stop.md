# NOTE 2026-11-23zx Caltech -> orion: t9l62 IA12.1 NaN at cycle 7926; STOPPING (USER rule: stop on NaN)

The first NaN is at cycle 7926, t2s 3, in the same sums as IA11's cycle-7822 failure (log.0 line 23536):
```
CLEAN MCAP cycle 7926 t2s 2: 4236 face momentum terms at the cap h chi_f E in the last pass
CLEAN FCONE cycle 7926 t2s 2 faces 25511442 momentum cut sum |dF| V/c^2 6.182656e+28
CLEAN FWB cycle 7926 t2s 2 faces 139692324 radiation momentum sum |dp| 3.252974e+29 x1 faces sum dp -4.216527e+27
CLEAN XT cycle 7926 t2s 2 write-back closure 8.748e-03
CLEAN OLDVEC cycle 7926 t2s 3 eint<0 87 (carried) E<0 0
CLEAN GNSW cycle 7926 t2s 3 pass 8: stalled, the gas Newton slope off (35 so far)
CLEAN cycle 7926 t2s 3 it -3 (sum 541) pos 0 mark 0 all 0 neg 0 r1 0 r2 201722 first 0 201722 guard 13 solves 200 cwin 0
CLEAN MCAP cycle 7926 t2s 3: 33273014 face momentum terms at the cap h chi_f E in the last pass
CLEAN FCONE cycle 7926 t2s 3 faces 38190253 momentum cut sum |dF| V/c^2 -nan
CLEAN FWB cycle 7926 t2s 3 faces 139692237 radiation momentum sum |dp| -nan x1 faces sum dp -nan
CLEAN WBFCT cycle 7926 t2s 3 rounds 2 left 13749516 kept -nan
CLEAN WBFDIAG cycle 7926 t2s 3 m 1 kji 112 116 91 rho 5.7642e-09 T 4.4046e+13 EO 1.3017e+07 X 9.8173e+07 E -2.9112e+07 Q 0.0000e+00 EI 2.1723e+05 d0 0.0000e+00 eg -5.3854e+10 wk 5.3896e+10 wp -1.7425e+02 dw -4.7827e+09 rd -1.2728e+08 den 4.2129e+07 take 0.0000
```
- Last cycle line before it: elapsed=6.019203e+03 cycle=7900 time=4.843363e+05 dt=1.874362e+00.
- STOP file written 19:32:48. The step ends with a final rst (post-NaN: do NOT use it).
- The IA12.1 run (prod_ia121_1123/t9l62) wrote no periodic rst before this. The last good start is still prod_ia11_1123/t9l62/rst_keep/t9l62.00012.rst.
- hpc-sm-02-09 becomes free. bsg's node 01-17 and N2 A/B are running the m35mms D1/D2/D1h arms now.
- Arms S/L (your 23zu reading) point at floor_pay_skip. Tell us the next t9l62 binary/args; a restart needs the USER's OK.

-- Caltech
