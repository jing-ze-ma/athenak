# NOTE 2026-11-21iu orion -> caltech: USER asks to STOP t19l64 (original IC); then C3h test from its rst on that node

1. USER (Orion session, ~01:40 PST): stop production t19l64 (J4, original IC, 02-04) — it is repeating the 11-20 runaway (KE 9e41,
   L_top 1.55, CW bursts). Per your rule the USER confirms in your session; please ask there and stop it with a final rst (keep its rsts).
   t19r (rebuilt) and all other production keep running.
2. Then on the freed node, a TEST (not production): build B_c3h = caltech-runs/c3h 7c9ebb817 (= joint4d + clean_stag_thinf_held_gate,
   front side only; default off) with CUDA 12.9 + g++ 13, incremental, smoke gate. Restart t19l64 from its latest rst with the J4 keys +
   rad_m1/clean_stag_thinf_held_gate = true + rad_m1/clean_nl_print = true, 300 cycles, rst output off (output*/last_time = 1e30 is not
   needed on Caltech). Report CW count, solves/stage, dt, KE, L_top vs the stopped run over the same cycles. Results to caltech_results/C3H/.
   (On viper, C3h on hegiant.00184 runs 36+ cycles with 0 CW and full dt where C1 has 7 CW and dt 6.)

-- Orion coordinator (for the USER)
