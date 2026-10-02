# NOTE 2026-10-02 (DeltaAI): BSG port gate -- one question on the column Picard mean (2.994 vs 3.000)

Reply to NOTE-2026-10-02-bsg-code-ready.md. Status report follows separately (NOTE-2026-10-0x-deltaai-bsg.md).

## Status
- **Builds at rt-integration 30bf6c03:** GPU athena_hes_gpu_30bf6c035ce5 (CUDA sm_90 + Grace, MPI)
  md5 2e9c93d04e08ba3189d6c23ad25cb446. CPU athena_hes_cpu_30bf6c035ce5 (host-only Grace, gcc 14, no MPI)
  md5 bbaf2931afdb911973d1f60d47e4c385.
- **Inputs:** the bundle inputs are byte-identical to `docs/handover/bsg_1001_bundle`. The 3-D input differs only in
  `<output1>` dt = 1000.
- **CPU column c1:** bsg_col_arm2, 1 rank, gate_compare.py:
  `rows 21  dt_end 1.243014e+02  Picard mean 3.000  L_top/L_in end 0.999466 range [0.999374, 0.999467]  FATAL 0  nan 0`
  Everything matches your numbers to the printed digit **except Picard mean: 3.000 here, 2.994 on viper.** On DeltaAI
  all 161 steps take exactly 3 passes; on viper one step takes 2.
- **Next steps:** the GPU gate (g1, g2, 20-cycle 3-D smoke) is queued. Production waits for it: a 24 h ghx4 job plus at
  most 3 interactive 2 h links, using the READY marker.

## What we checked
- **Stopping rule:** vet_col has fixcl, so implicit_conv_est is on. A step stops at pass 1 once
  q = res1/res0 < 0.5 and res1*q/(1-q) < implicit_tol (1e-8).
- **DeltaAI margin (`implicit_picard_log`):** the smallest value of that bound over the 161 steps is
  **5.6e-8, at step 16** (res0 3.15e-3, res1 1.32e-5, q 4.2e-3). Next are steps 39 and 41 (res1 about 3.1e-5).
- **Looser tolerance:** with `implicit_tol = 6e-8`, 2 steps stop at 2 passes (mean 2.987). dt_end, L end and the
  range are unchanged in every printed digit, and c1 vs this run differs by at most 5e-11 in the hst.
  **Our reading:** a borderline flip of the pass-1 estimate. res1 follows an inexact BiCGStab solve (EW 1e-2, about
  200 breakdowns and 14 line-Jacobi fallbacks per run), so it can differ by a factor of a few between x86 and Grace.
  We see no sign of a port problem.

## Question for viper
Please confirm with one viper CPU run of bsg_col_arm2 at 30bf6c03 with `implicit_picard_log = 100000`.
Set it in the `<rad_m1>` block by editing in place: a second key added at the top of the block is overridden by the
later line. We need:
1. **Which step takes 2 passes** on viper, and its res0, res1 and the bound res1*q/(1-q). We expect step 16, 39 or 41.
2. **Whether you agree** that we treat the Picard mean as report-only in the gate. The pass criteria stay dt_end,
   L_top/L_in and c1 vs GPU <= 1e-6.
