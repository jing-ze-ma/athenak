# NOTE 2026-09-29 (viper): m1-perf-0928 HIP build + He-box check (TASK-2026-09-29-viper-m1-perf-build.md)

**Status: build OK, check PASSED.** `fork/m1-perf-0928` 3058c6bd built incrementally with
`builds/build_inc_viper.sh box_gpu72` (gcc/16 rocm/7.2 hipcc 7.2.4, openmpi_gpu/5.0, GFX942_APU,
HIP_MALLOC_ASYNC=OFF; 31 objects, 170 s, **no HIP compile error or warning** in `rad_m1_{implicit,krylov}`),
binary `/viper/ptmp2/jinma/builds/bin/athena_box_gpu72_3058c6bd`, md5 523acc75a8f5bd1a73b27256c595e482.
He box = `bench-2026-09-29-hebox` input (fresh start, hesdirk2, cfl 0.3), `time/nlim=200 output1/dt=1.0`,
apudev job 12020920 (1 node, HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1), arms interleaved a b a b at 1 GPU
then 2 GPUs: (a) input as shipped (new default `implicit_opac_newton` on), (b) `implicit_opac_newton = false`
added under `<rad_m1>`. All 8 runs rc 0, FATAL 0, NON-CONVERGED 0, fallbacks 0, nan 0. **Picard mean/solve
2.030 (a) vs 6.759 (b)** (399 solves; max 4 vs 10), BiCGStab it/solve 11.26 vs 45.65. **ms/cycle** (median of
10-cycle windows, cycles 50-200, r1/r2): 1 GPU **48.97/48.87 (a) vs 117.31/117.73 (b) = -58 %**; 2 GPUs
**39.75/39.41 (a) vs 98.66/99.46 (b) = -60 %** (dt 0.16125 identical). **hst (a) vs (b)**, max over 34 rows
of |a-b| / max|a|: mass 2.0e-16, tot-E 4.8e-11, 1-mom 2.9e-8, 1-KE 2.7e-7, 2-/3-KE 2.5e-8; the near-zero
2-/3-mom differ by 1.8e-4 / 1.3e-4 of their own tiny scale (= 1e-14 of the 1-mom scale, noise); the same at
2 GPUs. Anomaly (minor): 1-KE 2.7e-7 is just above the DeltaAI 2e-7 figure, consistent with solver-tolerance
differences (a vs a at 1 vs 2 GPUs: 1-KE 2.7e-10; a repeat r1 vs r2: bitwise). Run tree, scripts
(`check.sub`, `cmp_hst.py`) and logs: `/viper/ptmp2/jinma/m1perf_check_0929/`.
