# NOTE 2026-11-23zi caltech -> orion: IA12.1 build fa59bffd6 + gate PASS; t9l62 restarted from 00012; first ~350 cycles: dt collapse to 0.013 + 15 DT_RETRY. Run NOT stopped (needs USER)

## Build + gate
- caltech-ia121-fa59bffd6 incremental (2 TUs), CUDA 12.9.41 + g++ 13.2, HOPPER90, -ffp-contract=off: **wt_orion1121/bin/athena_ia121_fa59bffd6 md5 9b4e69b95745378ec87bbbfe2287667c**
- Tiny BSG gate: **PASS**, pass 1 3.686e-11, converged after 3 passes, F/c 2.9006e+42 (caltech_results/H200CHK/gate_ia121_fa59bffd6.txt).

## t9l62 restart (17:47, step 7333072.39, hpc-sm-02-09, 4 ranks)
- `-r prod_ia11_1123/t9l62/rst/t9l62.00012.rst` (read in place; rst/ and rst_keep untouched) `-d prod_ia121_1123/t9l62 -t <END-10min>` + the keys of
  IA121_t9l62_restart.args (md5 eed6ae54 OK). Full command line in caltech_results/PROD_IA121_1123/t9l62/chunks.txt. Launcher prod_ia121_1123/launch_t9l62.sh
  (KEEPRST 3, only in the new dir). All 13 absent keys reported "added from the command line"; 0 FATAL.

## First cycles (to cycle 6800, t 4.7632e5, 18:17)
| | value |
|---|---|
| CLEAN WARNING | **35** |
| stage fallbacks | **21** |
| DT_RETRY | **15** (cycle 6469 x3, 6498 x2, **6577 x8**, 2 recovered) |
| FATAL / NaN | 0 / 0 |
| floor_pay "energy unpaid" | 0 |
| MCAP lines | 723 |
| dt | 3.52 (start) -> 0.161 (c6500) -> **min 0.0130 at t 4.75395e5** -> recovered to 7.61 (c6800) |
| s/cycle | 4.57 (c6500-6800, incl. retries; IA11 2.3) |
| Ltop/Lin | 2.556 -> 2.544 |
| Lbot/Lin | 0.533 -> **0.373** (falling, same trend as IA9 arm A) |
| GPU memory | 104.3-104.6 GB per GPU |

- Also at cycle 6577: `implicit_pos_floor created 7.842411e+41 of E (neither the gas nor the column covered the raise)`. Context:
```
CLEAN WBFDIAG cycle 6577 t2s 2 m 0 kji 77 48 53 rho 6.4526e-11 T 1.2530e+15 EO 2.6032e+08 X 3.6881e+08 E -2.7055e+11 Q -2.6568e+08 EI 1.3114e+05 d0 0.0000e+00 eg -9.8914e+03 wk 2.7107e+11 wp 5.3914e+01 dw -2.7107e+11 rd -2.7091e+11 den 2.7107e+11 tak
CLEAN WBFDIAG cycle 6577 t2s 2 m 0 kji 77 48 54 rho 5.2284e-11 T 5.2977e+15 EO 1.0008e-04 X 7.8800e+07 E -1.9028e+11 Q -2.1487e+08 EI 4.3285e+05 d0 0.0000e+00 eg -1.7375e+04 wk 1.9050e+11 wp 2.1347e+01 dw -1.9049e+11 rd -1.9036e+11 den 1.9050e+11 tak
CLEAN WBFDIAG cycle 6577 t2s 2 m 0 kji 77 48 56 rho 3.0206e-08 T 3.3412e+13 EO 1.0000e-04 X 4.0699e+07 E -1.6674e+04 Q -8.4967e+05 EI 2.8308e+07 d0 0.0000e+00 eg -8.6559e+07 wk 8.7428e+07 wp -1.7854e+03 dw -9.0280e+07 rd -4.0716e+07 den 8.6813e+05 ta
CLEAN XT cycle 6577 t2s 2 write-back closure 2.403e+03
CLEAN WARNING cycle 6577 t2s 2: implicit_pos_floor created 7.842411e+41 of E (neither the gas nor the column covered the raise)
CLEAN WARNING cycle 6577 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 7)
CLEAN GNSW cycle 6577 t2s 1 pass 10: stalled, the gas Newton slope off (9 so far)
CLEAN WARNING cycle 6577: backward-Euler solve not converged (dt_retry, 8 so far; last change 6.324e-05)
```
- Largest floor_audit floor_de line so far: `<hydro> floor_audit cycle 6823 [update] cell-c2p 4.979477709e+10 e_c<=0 2554 kedom 4376 sub_tfloor 13994543 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 13994543 (kedom 2568) floor_de 1.785496353e+41 <r> 5.514832378e+13 <r>_de 2.340640819e+12 skip(dfl`
- Files (caltech_results/PROD_IA121_1123/t9l62/): wbfdiag_to_c6800.txt (16 CLEAN WBFDIAG lines), floor_audit_to_c6800.txt (1104 lines),
  cw_dtretry_to_c6800.txt, hst, chunks.txt. Later WBFDIAG / floor_audit through 4.9e5 come with the watch.

## Alarm status
By your 23zf criteria this start is already an ALARM (dt < 1, DT_RETRY > 3/h). It has recovered since (dt 7.6, no retry since cycle 6577).
**I did NOT stop the run**: stopping a production needs the USER's direct OK, which COORD is asking for. The watch (prod_ia121_1123/t9mon_ia121.sh -> NOTE-2026-11-23zg-caltech-orion-ia121-t9l62-watch.md) is running from now on.
It posts a 2-h line, and posts + exits on dt < 1, DT_RETRY > 3/h, NaN/inf or unpaid floor energy.

## Other
- bsg (IA11) FINISHED at 17:29 PST at tlim (cycle 32156, t 4.92e6, rc 0). Its node hpc-sm-01-17 is free for scaling t9l62 to 8 H200 if the USER wants it.
  For now that node is running the m35mms dt diag (done, NOTE 23zh) and then the M35ZAMS 128^2 smoke (started 18:05).

-- caltech
