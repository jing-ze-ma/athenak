# NOTE 2026-11-01 Caltech -> viper (USER): face-field halo fix -18.6 % on r1m (bitwise); fused Strang RKG + ck Strang in progress

Seen: fin51-r-fix = aea58d13 (followup-1102 merged) and tusplit2-test-1102 (code motion; not tested here). Strang results: e975702d / 72d48003.

## 1. halofc-1102 @ d999f5ee (pushed, NOT merged; base aeafae97)
- **Cause of the r1m halo cost.** The cubed-sphere face-field pack (PackAndSendFC) takes 1.2-1.5 ms/call x 23 calls/cycle = 27.6 ms of 81 (34 %).
  - ncu: 2720 B stack + 128 registers per thread, because the nested by-reference seam lambdas spill to local memory.
  - 13 % occupancy. A few teams carry ~3600 seam samples each while most are empty.
  - It is not bandwidth and not launch latency. mhdspeed A made it worse by repeating that setup per chunk.
- **Fix.** The segment table is built once per buffer (the old per-team setup, verbatim). Then two flat kernels, one thread per sample:
  - FCPlainPack: 9.7 us/call, 81 % occupancy;
  - FCSeamPack: ~335 us/call, 24 % occupancy.
  - Cubed sphere only; every other mesh keeps the old kernel text unchanged.
  - Expressions are copied in form and order, under --fmad=false.
- **Gates** (H200, all against the 0abb64b6 references):
  - r1m 150 + 300 cycles: rst + mhd.hst IDENTICAL;
  - 150+150 continuation: rst IDENTICAL;
  - hydro 150: rst + hydro.hst IDENTICAL;
  - 2-GPU MPI, 100 cycles: IDENTICAL;
  - CPU 3 ranks: IDENTICAL. RKG stage count was 0 there, so on CPU the pack is exercised only in the RK halos.
  - user.hst differs, as it does between two base runs (floor-ledger order).
- **Timing** (1 H200, off_0p50.rst, ~300 cycles/arm, interleaved x2, empty GPU): 82.71 / 82.98 -> 68.21 / 66.72 ms/cycle = **-18.6 %**.
- **Please:** run the HIP/MI300A bitwise gate of d999f5ee (key-free path, cs MHD) and merge it into fin51-r-fix if it passes.
  It is not yet rebased onto aea58d13. It touches only src/bvals/bvals.hpp and bvals_fc.cpp, so we expect a clean merge.

## 2. Strang follow-up (USER)
- hfcstrang-1102 @ fcc40f98 = halofc d999f5ee + rkgstrang e204131b (clean merge, pushed).
- **In progress here: `<mhd>/rkg_strang_fuse`.**
  - The second RKG half of cycle n and the first half of cycle n+1 are fused into one super-step over (dt_n+dt_{n+1})/2.
  - Unfused halves are used only at the first and last cycle, before outputs and restarts, and around dt_retry.
  - Result: Strang at about Lie cost, 1 RKG per cycle instead of 1.41x.
- **Also in progress: `<problem>/ck_strang`.**
  - Today the ck RT is a Lie split: user_split_once, full dt after the RK stages, plus ck_impl_every 4.
  - Design: S = ck(dt/2) . RK . ck(dt/2) inside, with the fused RKG outside. That is 2nd order in both splits for 1 RKG + 2 ck calls per cycle.
  - Accuracy arms full / full4 against lie / strang / strang4 / seed; timing lie / strang / fused / full.
- **Why ck too.** Our test showed that Strang on the RKG alone cuts the r1m error only 10-60 % at production cfl. The ck Lie split, the ck cadence and the tracer source are other first-order pieces; the full arm shows which of them remain.
- **r1m: NOT launched.** The user is waiting for these results.
