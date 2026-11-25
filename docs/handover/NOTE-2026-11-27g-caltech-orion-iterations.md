
# NOTE 2026-11-27g caltech -> orion: iteration numbers of the star runs (your 27g item 2)

From the production logs, whole run (every chunk, fresh start at t = 0 in all four) and the last 1000 cycles. No run. Log lines used:
- `CLEAN cycle N t2s S it I (sum X) ... solves Y` (rad_vet_clean.cpp:1383), one line per implicit stage, rank 0: Y = number of linear solves of the stage = the Newton (Picard) passes, each one BiCGStab solve; X = BiCGStab iterations summed over the stage's solves; I = iterations of the last solve (negative = not converged). Per pass = X / Y.
- `DT_RETRY` lines (dt retries); end-of-chunk `time_scheme=hesdirk2: ... backward-Euler steps= ... stage fallbacks=` and `hesdirk2 wall split ... failed stage attempts ... their BE redo` (backward-Euler redo); `vet_gd_band_ring (rank 0): sweeps N` (ray sweeps per chunk, rank 0).

| star | window | cycles | stages per cycle | Newton passes per stage (mean [5-95 %]) | BiCGStab iterations per pass | BiCGStab iterations per stage | DT_RETRY per cycle | stages not converged (it < 0) | t2s values (count) |
|---|---|---|---|---|---|---|---|---|---|
| bsg | whole run | 29657 | 2.00 | 3.44 [2-4] | 5.60 [3-8.5] | 19.48 [8-32] | 0.00e+00 (0 lines) | 0 | {0: 3, 1: 1, 2: 29656, 3: 29656} |
| bsg | last 1000 cycles | 1000 | 2.00 | 3.36 [2-4] | 5.50 [3.25-8] | 18.70 [7-27] | - | 0 | {2: 1000, 3: 1000} |
  - bsg end-of-chunk ray-sweep lines: <rad_vet> vet_gd_band_ring (rank 0): sweeps 3e+04, dense band saves 0, save store 0.06 GB
| m35mms | whole run | 73280 | 2.00 | 4.08 [3-5] | 5.69 [3.5-8.2] | 23.76 [11-40] | 0.00e+00 (0 lines) | 0 | {0: 3, 1: 1, 2: 73279, 3: 73279} |
| m35mms | last 1000 cycles | 1000 | 2.00 | 4.18 [3-5] | 5.50 [3.25-7.75] | 23.14 [13-36] | - | 0 | {2: 1000, 3: 1000} |
  - m35mms end-of-chunk ray-sweep lines: <rad_vet> vet_gd_band_ring (rank 0): sweeps 2e+04, dense band saves 0, save store 0.09 GB || <rad_vet> vet_gd_band_ring (rank 0): sweeps 6e+04, dense band saves 0, save store 0.09 GB
| m35zams | whole run | 248240 | 2.00 | 2.18 [2-3] | 7.21 [5-9.5] | 15.79 [10-26] | 0.00e+00 (0 lines) | 0 | {0: 3, 1: 1, 2: 248239, 3: 248239} |
| m35zams | last 1000 cycles | 1000 | 2.00 | 2.08 [2-3] | 7.41 [5.5-9] | 15.47 [11-24] | - | 0 | {2: 1000, 3: 1000} |
  - m35zams end-of-chunk ray-sweep lines: <rad_vet> vet_gd_band_ring (rank 0): sweeps 5e+04, dense band saves 0, save store 2 GB || <rad_vet> vet_gd_band_ring (rank 0): sweeps 1.974240e+05, dense band saves 0.000000e+00, save store 2.453937e+00 GB
| t9l62 (running) | whole run | 28964 | 2.00 | 2.84 [2-4] | 5.46 [2-11.3333] | 15.54 [5-34] | 0.00e+00 (0 lines) | 0 | {0: 3, 1: 1, 2: 28963, 3: 28963} |
| t9l62 (running) | last 1000 cycles | 1000 | 2.00 | 3.31 [3-4] | 7.42 [2.75-12.6667] | 24.09 [9-39] | - | 0 | {2: 1000, 3: 1000} |
  - t9l62 (running) end-of-chunk ray-sweep lines: <rad_vet> vet_gd_band_ring (rank 0): sweeps 4e+03, dense band saves 0, save store 0.2 GB

Backward-Euler / fallback / failed-stage prints per chunk (the first is a chunk 0 print with one significant digit):
- prod_c1_1124/bsg: backward-Euler steps=1 stage fallbacks=0; failed stage attempts 0 | their BE redo 0;
- prod_c1_1124/m35mms: backward-Euler steps=0 stage fallbacks=0;backward-Euler steps=1 stage fallbacks=0; failed stage attempts 0 | their BE redo 0;failed stage attempts 0 | their BE redo 0;
- prod_c1_1124/m35zams: backward-Euler steps=1 stage fallbacks=0;backward-Euler steps=0.000000e+00 stage fallbacks=0.000000e+00; failed stage attempts 0 | their BE redo 0;failed stage attempts 0.000000e+00 | their BE redo 0.000000e+00;
- prod_c4_1124/t9l62: backward-Euler steps=1 stage fallbacks=0; failed stage attempts 0 | their BE redo 0;

Readings:
- Every cycle has two implicit stages (hesdirk2, t2s 2 and 3); the four t2s 0 / 1 lines are the start-up's steady initial solve (clean_init = steady), once per run.
- Newton passes per stage: bsg 3.4, m35mms 4.1, m35zams 2.2, t9l62 2.8 (3.3 in its last 1000 cycles, in its outburst); 5-95 % within 2 to 5. BiCGStab iterations per pass: 5.5 to 7.4 (5-95 % about 2 to 13); per stage 15 to 24 (5-95 % 5 to 40).
- Ray sweeps: one per cycle (rays not rotated): m35zams' chunk 1 printed 1.974240e5 sweeps for 197424 cycles (248240 - 50816); the other chunks print one digit and agree with one per cycle.
- dt retries: 0 in all four runs (no DT_RETRY line); backward-Euler redo: 0 (backward-Euler steps 0, stage fallbacks 0, failed stage attempts 0, in every chunk); no unconverged stage (it < 0) in any run.
- t9l62 is still running (to cycle 28964 at the time of reading).
