# NOTE 2026-11-23zq Caltech -> orion: t9l62 IA12.1 implicit_pos_floor **5.52e47** at cycle 7605 (still running, USER)

log.0 line 18341:
```
CLEAN WBFDIAG cycle 7605 t2s 3 m 0 kji 74 50 59 rho 5.2457e-11 T 4.8998e+13 EO 2.4354e+08 X 1.9244e+08 E -4.2648e+12 Q -5.7658e+07 EI -1.4837e+05 d0 0.0000e+00 eg -4.6114e+14 wk 4.6540e+14 wp -2.3503e+00 dw -4.6540e+14 rd -4.2650e+12 den 4.2651e+12 take 0.0000e+00 dm 2.7474e-06 -2.2099e+02 1.1922e-0
CLEAN WBFDIAG cycle 7605 t2s 3 m 0 kji 74 49 60 rho 1.0691e-11 T 2.1519e+17 EO 3.0978e+08 X 4.0752e+08 E -1.4052e+16 Q 0.0000e+00 EI 2.0840e+09 d0 0.0000e+00 eg -2.0805e+09 wk 1.4052e+16 wp -1.4590e-02 dw -1.4052e+16 rd -1.4052e+16 den 1.4052e+16 take 0.0000e+00 dm 5.4813e+02 -1.4859e-04 -2.6487e-05
CLEAN WBFDIAG cycle 7605 t2s 3 m 0 kji 74 51 59 rho 7.3218e-12 T 4.3574e+13 EO 2.5363e+08 X 1.2018e+08 E -1.1385e+09 Q 0.0000e+00 EI -2.8382e+03 d0 0.0000e+00 eg -1.5354e+10 wk 1.6746e+10 wp -2.2015e+00 dw -1.6745e+10 rd -1.2586e+09 den 1.3921e+09 take 0.0000e+00 dm 8.1353e-09 4.9347e-01 -3.8338e-05
CLEAN WBFDIAG cycle 7605 t2s 3 m 0 kji 75 49 59 rho 2.4714e-09 T 5.0622e+13 EO 1.0432e+08 X 2.1877e+08 E -1.0580e+13 Q -5.6919e+07 EI -6.3975e+06 d0 0.0000e+00 eg -2.7753e+14 wk 2.8811e+14 wp -1.1690e+03 dw -2.8811e+14 rd -1.0581e+13 den 1.0580e+13 take 0.0000e+00 dm -1.8311e-04 -1.1555e-03 -1.1938e
CLEAN WBFDIAG cycle 7605 t2s 3 m 0 kji 76 49 59 rho 3.9099e-08 T 2.5848e+13 EO 4.7613e+07 X 1.4664e+07 E -8.7541e+05 Q -4.0478e+06 EI 3.1056e+06 d0 0.0000e+00 eg -1.6821e+10 wk 1.6873e+10 wp 4.8747e+02 dw -1.6749e+10 rd -1.5539e+07 den 5.2535e+07 take 0.0000e+00 dm 3.0542e-03 -6.0391e-04 -3.6691e+01
CLEAN XT cycle 7605 t2s 3 write-back closure 4.415e+07
CLEAN WARNING cycle 7605 t2s 3: implicit_pos_floor created 5.524926e+47 of E (neither the gas nor the column covered the raise)
CLEAN WARNING cycle 7605 t2s 3: stage solve failed, the step is redone with backward Euler (stage fallback 14)
CLEAN cycle 7605 t2s 1 it 5 (sum 118) pos 0 mark 0 all 0 neg 0 r1 0 r2 3668 first 0 3668 guard 0 solves 9 cwin 0
CLEAN MCAP cycle 7605 t2s 1: 5860 face momentum terms at the cap h chi_f E in the last pass
CLEAN FCONE cycle 7605 t2s 1 faces 26544369 momentum cut sum |dF| V/c^2 1.488033e+30
```
- That is ~4x the hydro.hst column-7 total (~1.335e47). The same stage is then redone with BE ("stage fallback 14", same cycle 7605, t2s 3).
- hydro.hst column 7 stays at 1.33501545e47-1.33501554e47 across these rows, so the created energy looks discarded with the redone stage, not kept. Please confirm.
- Warnings since 7457: stage fallbacks 13-16 (7602, 7605, 7607, 7616). 0 NaN, 0 FATAL. At cycle 7600: t 4.8371e5, dt 5.40.
- Still running per the USER. We stop on dt < 0.01, NaN/inf or FATAL. If you think it should stop, say so and we'll ask the USER.

Last 3 hydro.hst rows:
```
   4.8380593821398274e+05   4.6483349543052279e+00   2.3310565725611246e+32   1.8764651824339854e+36   2.0201122709619678e+34   4.6087101579947073e+34   1.3350157015158322e+47   1.9889921032895589e+44
   4.8380708316368231e+05   1.1449496995572788e+00   2.3310565725611246e+32   1.8767811124140103e+36   2.0190035052492155e+34   4.6087864195999183e+34   1.3350152479584901e+47   1.9889996877082129e+44
   4.8380937306308141e+05   2.2898993991145575e+00   2.3310565725611282e+32   1.8774487344495968e+36   2.0168205010746609e+34   4.6089187212286442e+34   1.3350152358404610e+47   1.9890142785024016e+44
```
-- Caltech
