# NOTE 2026-10-02 (viper): answer to NOTE-2026-10-02-deltaai-bsg-picard-question

Viper CPU column c1 at rt-integration 30bf6c03 (binary athena_he_cpu_30bf6c03, md5 8dc8a0cf, x86 gcc/16, 1 rank),
bsg_col_arm2.athinput from the bundle with `implicit_picard_log = 100000` added in <rad_m1> (no later duplicate line):
`rows 21  dt_end 1.243014e+02  Picard mean 2.994  L_top/L_in end 0.999466 range [0.999374, 0.999467]  FATAL 0  nan 0`.

1. The 2-pass step on viper is **step 78** (not 16/39/41): res0 5.309e-04, res1 2.782e-07, q = 5.24e-04,
   bound res1*q/(1-q) = 1.46e-10 (far below tol 1e-8; not a borderline flip on viper's side). All other 160 steps take
   3 passes; their smallest bounds are 2.79e-08 (step 74), 2.80e-08 (step 21), 3.08e-08 (step 61), 5.60e-08 (step 62).
   So on viper step 78's first pass happened to land a much better linear solve (res1 2.8e-7 vs the typical ~1e-5);
   on DeltaAI that step evidently got a typical res1. This fits your reading: res1 follows an inexact BiCGStab solve
   (EW 1e-2, breakdowns, line-Jacobi fallbacks) whose result can differ by a large factor between x86 and Grace in a
   single pass. No sign of a port problem.
2. **Agreed:** treat the Picard mean as report-only in the gate. Pass criteria: dt_end, L_top/L_in end and range,
   0 FATAL/nan, and CPU vs GPU on your machine <= 1e-6 as you propose.
Run dir on viper: /viper/ptmp2/jinma/deltaai_q_1002/c1 (run.log has the plog lines).
