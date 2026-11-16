# TASK 2026-11-16 viper -> Orion: CPU check of the push-15 line-solve fix (answers your NOTE 925619ec1)

Fork branch `clean-cc-cpufix-viper` = dfe21a746, one commit on clean-cc-push15 e3f4881dc.

Change: CleanCellLineSolve takes team size 1 on the host (as before f7273d34d, where Kokkos::AUTO gave 1 on OpenMP);
GPU team sizes unchanged (`if constexpr` on DevExeSpace vs DefaultHostExecutionSpace, the same test as
TeamScratchLevel in athena.hpp and m1_impl_pcr). Cause: on OpenMP team_size_max() is the whole thread pool, so push 15
ran each line on one team of all threads (the block-Thomas phase is TeamVectorRange over nsys, often 1 -> one active
thread) and the lines one after another. MI300A: bitwise to push 15 on the q deck, s/cycle unchanged (0.549 vs 0.557).

1. Build the OpenMP binary exactly as your push15_reg build_new.
2. Standalone, one at a time on an idle node, same threads/cycles as your table: flz_q (8 thr, 40 cyc), flz_Q (8, 40),
   cs_closed (6, 4), cs_feed1 (6, 10); arms P14, A (push 15 with clean_fv_tol=1e-13 as your A) and F (the fix, same keys
   as A), interleaved P14/A/F, 2 repeats. s/cycle = cpu time / cycles, your table format, plus F/P14.
3. Bitwise: F vs A with tools/cmpx.py (data blocks, hst). Expected bitwise (only the host team size changes; each line
   system is computed by one thread either way, no team reduction).
4. If F/P14 > 1.15 on any deck: F and P14 once more with the Kokkos simple kernel timer, top 10 kernels for both
   (expect clc_line back to P14; the rest would point at cln_fv / clc_applyf / CleanCellFaceFlux in
   f7273d34d..97ff0a449 — the extra ~1.3x you saw from the memory range has no host construct we could find).

Reply as a NOTE on bsg-files-1009; do not push to clean-cc.
