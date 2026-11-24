# NOTE 2026-11-23zp Caltech -> orion: t9l62 IA12.1 implicit_pos_floor 9.96e42 at cycle 7457; USER says KEEP RUNNING

log.0 line 16040:
```
CLEAN WBFCT cycle 7457 t2s 3 rounds 6 left 5 kept 1.914113e+40
CLEAN WBFDIAG cycle 7457 t2s 3 m 0 kji 72 45 63 rho 1.1621e-10 T 3.1617e+16 EO 6.2421e+07 X 1.9416e+07 E -4.6107e+08 Q -1.7265e+08 EI 3.0849e+07 d0 0.0000e+00 eg -2.5337e+07 wk 7.2149e+08 wp 1.6461e+01 dw -7.2126e+08 rd -4.8049e+08 den 6.9615e+08 take 0.0000e+00 dm 3.2849e-05 1.1397e-04 -3.7666e-01 
CLEAN WBFDIAG cycle 7457 t2s 3 m 0 kji 73 45 62 rho 5.5413e-12 T 9.4967e+16 EO 4.5660e+07 X 3.8267e+07 E -2.9043e+12 Q -1.9088e+08 EI 9.9137e+07 d0 0.0000e+00 eg -9.8348e+07 wk 2.9046e+12 wp 4.9417e-01 dw -2.9046e+12 rd -2.9043e+12 den 2.9045e+12 take 0.0000e+00 dm -5.6721e+00 -3.4588e-06 -3.7497e-0
CLEAN WBFDIAG cycle 7457 t2s 3 m 0 kji 73 44 63 rho 3.7042e-11 T 5.8229e+16 EO 9.7690e+07 X 6.8296e+07 E -2.5326e+10 Q -3.8042e+07 EI 9.6154e+07 d0 0.0000e+00 eg -9.2919e+07 wk 2.5555e+10 wp 1.3338e+00 dw -2.5555e+10 rd -2.5394e+10 den 2.5462e+10 take 0.0000e+00 dm 1.3244e-05 1.3867e+00 -2.4586e-05 
CLEAN WBFDIAG cycle 7457 t2s 3 m 0 kji 74 45 63 rho 1.3330e-10 T 3.1271e+16 EO 6.0954e+07 X 1.6420e+08 E -4.2320e+08 Q 0.0000e+00 EI 3.4065e+07 d0 0.0000e+00 eg -2.7813e+07 wk 5.1197e+08 wp 2.8143e+01 dw -5.1377e+08 rd -5.8740e+08 den 4.8415e+08 take 0.0000e+00 dm 1.6387e-04 3.1397e-04 -4.0669e-01 v
CLEAN XT cycle 7457 t2s 3 write-back closure 1.587e+04
CLEAN WARNING cycle 7457 t2s 3: implicit_pos_floor created 9.955100e+42 of E (neither the gas nor the column covered the raise)
CLEAN WARNING cycle 7457 t2s 3: stage solve failed, the step is redone with backward Euler (stage fallback 9)
CLEAN cycle 7457 t2s 1 it 3 (sum 60) pos 0 mark 0 all 0 neg 0 r1 0 r2 4341 first 0 4341 guard 0 solves 7 cwin 0
CLEAN MCAP cycle 7457 t2s 1: 5186 face momentum terms at the cap h chi_f E in the last pass
CLEAN FCONE cycle 7457 t2s 1 faces 27732626 momentum cut sum |dF| V/c^2 3.607324e+29
```
- This is the 2nd "implicit_pos_floor created" in the run (the 1st is 6577, 7.8e41).
- Last cycle line: elapsed=3.753120e+03 cycle=7450 time=4.827008e+05 dt=9.747017e+00. There are 37 CLEAN WARNINGs in total, 0 NaN, 0 FATAL.
- By 23zk this is a stop condition. The USER decided on Caltech to **keep it running**. We still stop on dt < 0.01, NaN/inf or FATAL, and we post (without stopping) any further pos_floor or a skip_de jump > 1e41 followed by a fallback.

-- Caltech
