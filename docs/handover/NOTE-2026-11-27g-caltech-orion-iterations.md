
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

## Update (15:55): t9l62 finished (15:51:34, cycle 29228), t9mlt so far; and a correction
| star | window | cycles | stages per cycle | Newton passes per stage (mean [5-95 %]) | BiCGStab iterations per pass | BiCGStab iterations per stage | DT_RETRY per cycle | stages not converged (it < 0) | t2s values (count) |
| t9l62 | whole run | 29228 | 2.00 | 2.84 [2-4] | 5.48 [2-11.3333] | 15.61 [5-34] | 0.00e+00 (0 lines) | 0 | {0: 3, 1: 1, 2: 29227, 3: 29227} |
| t9l62 | last 1000 cycles | 1000 | 2.00 | 3.29 [3-4] | 7.51 [3-12.6667] | 24.21 [10-39] | - | 0 | {2: 1000, 3: 1000} |
| t9mlt (running) | whole run | 22648 | 2.00 | 3.54 [2-5] | 4.86 [2.75-7.33333] | 17.33 [7-28] | 0.00e+00 (0 lines) | 0 | {0: 3, 1: 1, 2: 22647, 3: 22647} |
| t9mlt (running) | last 1000 cycles | 1000 | 2.00 | 3.70 [3-4] | 5.44 [3.6-7.5] | 20.24 [12-30] | - | 0 | {2: 1000, 3: 1000} |
  - t9l62 end-of-chunk ray-sweep lines: <rad_vet> vet_gd_band_ring (rank 0): sweeps 4e+03, dense band saves 0, save store 0.2 GB || <rad_vet> vet_gd_band_ring (rank 0): sweeps 3e+04, dense band saves 0, save store 0.2 GB
  - t9mlt (running) end-of-chunk ray-sweep lines: none
- Correction: the end-of-chunk prints say backward-Euler steps = 1 in each chunk (not 0 as I wrote): one backward-Euler step at the start of each chunk (the first step of a fresh start or of a restart; the four t2s 0 / 1 lines above are the steady initial solve of the fresh start). Stage fallbacks 0, failed stage attempts 0, their BE redo 0 in every chunk. So the backward-Euler redo fraction is 0, and the backward-Euler steps are 1 per chunk start (bsg 1 in 29657 cycles, m35mms 2 in 73280, m35zams 2 in 248240, t9l62 2 in 29228).
- t9mlt has no end-of-chunk print yet (running); its CLEAN cycle lines give 3.5 Newton passes per stage, 4.9 BiCGStab iterations per pass, 17 per stage, 0 dt retries, 0 unconverged stages so far.

- Correction of the correction: per chunk the prints are bsg log.0 1; m35mms log.0 1, log.1 0; m35zams log.0 1, log.1 0; t9l62 log.0 1, log.1 0. So exactly ONE backward-Euler step per run, in the first chunk (the fresh start's first step), and none after the restarts; the per-run numbers I gave just above (2 for m35mms, m35zams, t9l62) are wrong: 1 each.
