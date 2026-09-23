# prof_0923: GPU cost profile of the correlated-k RT (semi vs implicit T4) against hydro

Goal (user): bring the implicit ck RT (T4) down to the cost of the hydro. This note is a
profile plus a ranked plan. No code was changed.

## Setup

- Binary: `git archive be02c647` (rt-integration HEAD), production recipe (HIP gfx942 APU,
  MPI, Release, PROBLEM=deep_hot_jupiter_rt). The kokkos symlink points to
  bench/wt_rgbox/kokkos. Binary: /viper/ptmp2/jinma/ckprof_0923/athena.gpu,
  md5 6f6f9908bc7174d8f7cb712340750262.
- State: bench/cs_hyd4_prod/rst/dhj.00124.rst, read in place. t = 1.891e7 s,
  **rotation 62.0**; md5 d1a55983b0cd2c105b2cb05dd9be02d6.
- Input: bench/cs_hyd4_prod/deep_hot_jupiter.athinput with the ck_impl keys added at
  their default values (/viper/ptmp2/jinma/ckprof_0923/prof.athinput). This is needed
  because the production restart was written by 1d22a118, which predates ck_impl_lin /
  fuse / jac_lin / cvsec, and AthenaK refuses command-line keys that are not already in
  the input. With that input and no overrides, the run is the production run.
- apudev, 2 ranks / 2 MI300A (228 CU, 912 SIMD, 8 waves/SIMD each), HSA_XNACK=1,
  HSA_NO_SCRATCH_RECLAIM=1. Mesh: 24 MeshBlocks of 128x16x16, 12 per rank, 3072 columns
  per rank.
- Arms (/viper/ptmp2/jinma/ckprof_0923/common.sh):
  - **h** = hydro only: `problem/user_srcs=false`. No input switch disables only the
    two-stream RT: SourceFunc always calls picket_fence_two_stream_RT. So this arm also
    drops the gravity/rotation source kernel (`usrsource`, 0.93 ms/cycle in the traces).
    It runs sanely for 220 cycles, but dt is 11.65 s against 18.1 s in the RT arms, so
    compare per cycle only.
  - **s** = production semi-implicit.
  - **t8** = T4 (ck_implicit, arat 1e30, frozen_op, lin, lin_thr 1, fuse, jac_lin, cvsec),
    tol 1e-8.
  - **t10** = t8 + tol = dtol = 1e-10, maxit 20.
  - **t2c / t2n** = frozen_op without lin, with frozen_cof on / off. This is the only way to
    isolate the expm1 cost: ck_impl_frozen_cof already defaults to true and ck_impl_lin
    refuses to run without it, so "T4 + frozen_cof" is T4.
- Timing: job 11953542 (time2/). Arms interleaved h s t8 t10 | t10 t8 s h | h s t8 t10, 120
  to 220 cycles each. ms/cycle comes from the `elapsed=` lines (ndiag = 20, the first 20
  cycles are skipped). t2c/t2n come from job 11953274 (time/), same binary and restart;
  the h/s arms of that job agree with job 11953542 to 1.5 %.
- Traces: rocprofv3 kernel + kokkos + HIP API + memcpy traces, about 21 cycles each.
  Kernels are attributed to Kokkos labels through the HIP launch correlation id. The
  window starts after the first cycle, and figures are rank 0. s/h/t2c/t2n: job 11953275
  (prof/); t8/t10: job 11953543 (prof2/). Scripts: ana_time.py, ana_prof.py, seq.py in
  /viper/ptmp2/jinma/ckprof_0923.

## 1. Cost table (untraced wall time, median of 3)

| arm | ms/cycle | x hydro | RT = arm - h (ms) | RT x hydro | passes/call (mean, max) | RT calls/cycle |
|---|---|---|---|---|---|---|
| h (hydro only) | 15.39 | 1.00 | - | - | - | 0 |
| s (semi, production) | 62.79 | 4.08 | 47.4 | 3.1 | 1 | 2 |
| t8 (T4, tol 1e-8) | 125.61 | 8.16 | 110.2 | 7.2 | 4.00, 5 | 2 |
| t10 (T4, tol 1e-10) | 142.77 | 9.28 | 127.4 | 8.3 | 5.96, 7 | 2 |
| t2c (frozen_op, cof on, no lin) | 382.9 | 24.9 | 367 | 23.9 | 5.49, 6 | 2 |
| t2n (frozen_op, cof off, no lin) | 413.6 | 26.9 | 398 | 25.9 | 5.49, 6 | 2 |

Repeat spreads: h 15.34 to 15.46, s 62.44 to 63.38, t8 125.5 to 127.5, t10 142.5 to 143.4.

GPU kernel busy time per cycle in the traces, split at the `rt_*` / `ck_*` labels:

| arm | busy | RT kernels | non-RT kernels | host gap (untraced wall - busy) | launches/cycle | hipDeviceSynchronize/cycle |
|---|---|---|---|---|---|---|
| h | 11.06 | 0 | 11.06 | 4.3 | 94 | 88 |
| s | 51.42 | 39.36 | 12.06 | 11.4 | 169 | 214 |
| t8 | 109.47 | 97.23 | 12.24 | 16.1 | 305 | 368 |
| t10 | 127.93 | 115.62 | 12.31 | 14.8 | 399 | 499 |

The non-RT GPU work is the same in every arm (12.1 to 12.3 ms; h is lower by the missing
usrsource). The RT excess over hydro therefore consists of the RT kernels plus about 7 to
12 ms/cycle of extra host gaps (launches and syncs).

## 2. Per-pass breakdown (ms per call, i.e. per RK stage, per rank; from seq.py on the traces)

**Semi** (2 sweeps per cycle): rt_pre_tp 0.18 + rt_pre_opac 0.27 + **rt_chain_ck 18.5** +
rt_apply 0.70 + cut/geom 0.05 = **19.7 ms per sweep**.

**T4, tol 1e-8** (4.00 passes per call, 2 calls per cycle; kernel sequence read from the
trace):

| component | pass 0 (storing) | each later pass (everything opacity-derived frozen) |
|---|---|---|
| (a) opacity lookups + expm1 + beam + store, inside rt_chain_ck | 23.5 (whole kernel) | 0 (removed by frozen_op + frozen_cof + lin) |
| ck_lin_build (factorisation store) | 3.0 | 0 |
| (b) Planck/T update: rt_pre_tp (+ rt_pre_opac, still launched) | 0.2 + 0.3 | 0.2 + 0.1 |
| (c) sweep arithmetic: rt_chain_ck_lin1 + lin_sum | - | 1.9 + 0.6 = 2.5 |
| (c') Jacobian: rt_chain_ck_jlin + ck_jlin_sum, rebuilt on EVERY pass | 2.6 + 0.2 | 2.6 + 0.2 = 2.8 |
| (d) tridiagonal + Newton step + reductions: ck_impl_fused | 0.3 | 0.3 (0.1 once the columns have converged) |
| apply de: rt_apply | 0.7 | 0.7 |
| **total** | **30.8** | **6.6** |

This gives 30.8 + 3 x 6.6 = 50.6 ms per call and 101 ms per cycle, against 97.2 ms of RT
kernels measured. **Pass 0 is 61 % of the T4 RT cost.** A pass with everything
opacity-derived frozen costs 6.6 ms, i.e. 1/3 of a semi sweep. 42 % of that pass is the
Jacobian rebuild and 38 % is the linear re-apply.

Inside pass 0 (from the t2 arms, whose chain kernel instantiation also assembles the
Jacobian):
- Storing pass: 39.4 ms. Frozen pass with stored coefficients (t2c): 25.0 ms. Frozen pass
  recomputing expm1 (t2n): 29.3 ms.
- So **expm1 = 4.3 ms per sweep** (what frozen_cof removes, 15 % of that sweep).
- Table lookups + beam + store = about 10 ms.
- T4's pass-0 kernel (23.5 ms) is the semi sweep (18.5 ms) plus about 5 ms of storing.
- The t2 frozen sweep (25 ms, full recurrences + atomics) is 10x the T4 linear kernel
  (2.5 ms). That ratio is what the lin kernel bought.

Column skipping (t8v hist, mean active columns out of 3072 after each pass): 3072, 3072,
2053, 4, 0. At tol 1e-10 (t10v): 3072, 3072, 3072, 2117, 735, 112, 0. The final pass of each
call only confirms convergence (0 active), yet jlin + rt_apply still cost their full 3.3 ms
on it. The skip removes work inside the chain kernels only, not the Jacobian or apply
launches.

Tolerance cost: going from 1e-8 to 1e-10 adds 1.96 passes per call, i.e. +17.2 ms/cycle
measured (4.4 ms per extra pass per stage, below 6.6 because of the skip in lin1).

## 3. ck table size and occupancy

- Table data/exo_fms_ck/ck/Premixed_1x_g8_11.txt: **11 bands x 8 g = 88 chains** (ck_nquad
  = 1), on a 38 T x 34 p grid. CK_NB and CK_NG are compile-time constants
  (correlated_k.hpp:74-75).
- rt_chain_ck: one thread per (MeshBlock, block of RT_NB = 4 chains, k, j). That is
  12 x 22 x 256 = 67 584 threads = 264 workgroups = **1056 waves on 912 SIMDs (1.16 waves
  per SIMD, 14 % of the 7296 wave slots)**. Private segment: 31 KB per lane (T4 / semi
  instantiation), 45 KB (t2 instantiation), i.e. a scratch-resident per-thread column.
- **The chain kernel sits at the occupancy floor.** Its time is the serial per-thread
  latency (4 chains x 128 cells of recurrences with scratch traffic), not throughput. The
  total chain count only sets how many threads run side by side.
- So cutting g-points (88 -> 44 chains, i.e. 528 waves) removes threads but not per-thread
  latency. The expected gain on rt_chain_ck is only the 1056-vs-912 wave tail, at most
  about 15 %. This is inferred from the launch geometry and was **not measured**: no
  reduced table exists. By the same argument, the gain would come from more threads per
  chain block, not fewer chains.
- rt_chain_ck_lin1 (lin_thr = 1, one thread per chain): 1056 workgroups = 4224 waves (4.6
  per SIMD), so it is closer to throughput-bound and would scale roughly with the chain
  count.
- Memory, rocm-smi "VRAM used" per GPU on the apudev node: 26.4 GB for t8/t10/t2n and
  27.8 GB for t2c. For s and h it reads 16 MB, so the APU counter does not see the semi
  allocations (the T4 operator stores do show up). On the running production cs_mhd_prod4
  (job 11941995): 6.76 GB per GPU and 14.0 to 14.6 GB host RSS per rank. Capacity is
  118 GB per GPU, so memory is not a constraint for any option below.

## 4. Ranked plan toward RT cost ~= hydro cost (15.4 ms/cycle)

Baseline: T4 1e-8 RT = 110 ms/cycle wall = 97 ms of RT kernels + about 13 ms of extra host
gap. The savings below are **estimates** built from the measured components above. None
was measured as a combination.

| # | item | mechanism | est. saving (ms/cycle, T4 1e-8) | notes |
|---|---|---|---|---|
| 1 | **RT once per cycle** (`ck_impl_once`, exists) | 2 RT calls/cycle -> 1 (rk2) | about 50 (-> ~55 RT) | already implemented; refused with rt_rad_force; needs an accuracy gate. Not timed here |
| 2 | **Opacity/operator frozen across steps** (the user's position), refreshed every k steps | skip pass 0 (30.8 ms per call) on non-refresh steps: a call becomes n x 6.6 ms | with #1 and k = 4: 55 -> 4 x 6.6 + 24.2/4 = 32 | needs the stored operator kept across calls (26 GB is fine) and a refresh rule; the beam deposit Qb is kept too |
| 3 | **Jacobian reuse within a call** (build once at pass 0, quasi-Newton; `ck_impl_reuse_jac` exists for the old path) | drop jlin + jlin_sum on passes >= 1: 2.8 ms each | 3 x 2.8 = 8.4 per call; with #1: ~8 | may add passes; the jac_lin agent (ck-jac) owns this |
| 4 | **Fewer passes**: better initial guess (`ck_impl_warm`, exists), or stop without the confirmation pass | the last pass has 0 active columns but costs 3.3 ms; 4 -> 3 or 2 passes | 3.3 to 10 per call | warm start is not in the restart; the confirmation pass could be replaced by a converged-count reduction inside the fused kernel |
| 5 | **Pass-0 kernel occupancy** (split the 4-chain thread block, RT_NB 4 -> 1, or a band-parallel pass 0) | 1.16 waves/SIMD, 31 KB scratch per lane; 4x the threads at 1/4 of the work each | up to ~15 of 23.5 per call if latency-bound (not measured) | also the semi production sweep (18.5 ms): benefits semi directly |
| 6 | **Launch/sync fusion** | 305 launches and 368 hipDeviceSynchronize per cycle vs 94/88 for hydro; host gap 16 vs 4 ms; merge pre_tp + pre_opac + apply + fused + sums per pass | about 5 to 10 | pre_opac is still launched on frozen passes (0.1 ms, avoidable) |
| 7 | **Column skipping at kernel granularity** | skip jlin/apply for converged columns; pass 3 has 2053/3072 active, the last pass 0 | about 2 to 4 per call | colskip already exists inside the chain kernels |
| 8 | **Reduced g-points** (8 -> 4 per band) | chain kernel at the occupancy floor | <= 15 % of pass 0; ~50 % of lin1 (2.5 ms/pass) | accuracy cost; low priority unless combined with #5 |
| 9 | **Mixed precision** (RT_FP32, today refused by ck_impl_lin) | on MI300A FP32 is ~2x FP64 vector rate, but the kernels are latency/scratch bound; halving the 31 KB scratch is the real gain | <= 20 % of the chain kernels | Newton at 1e-10 in FP32 is doubtful; last |
| 10 | **RT subcycling** (RT every k hydro steps, heating rate held) | cost / k on top of #1 | with #1-#4 and k = 2: halves the rest | the implicit solve is stable, so this is accuracy-limited only; needs a gate against v_MLT/energy budget |

**Combination that reaches ~hydro.** From 110 ms:
1. #1 once per cycle: ~55.
2. #2 operator frozen for k = 4 steps: ~32.
3. #3 Jacobian reuse: ~26.
4. #4 one pass fewer (drop the confirmation pass): ~19.
5. #6 fusion/sync: **~12 to 15 ms/cycle, i.e. RT ~= 0.8 to 1.0 x hydro.**

Items #5 and #10 are the reserve (#10 with k = 2 alone would halve it again). Without
freezing the operator across steps (#2), pass 0 alone (2 x 30.8 per cycle, or 30.8 once per
cycle) stays at 2 to 4x hydro. So **#1 + #2 are the essential pair**, and #5 is the only
route that also cuts the semi production cost (47 ms/cycle of RT today, 3.1x hydro).
