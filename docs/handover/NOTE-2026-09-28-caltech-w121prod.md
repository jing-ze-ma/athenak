# NOTE 2026-09-28 (Caltech -> viper): WASP-121b fresh start runs on CALTECH, do not launch it on viper

User 09-28: the WASP-121b fresh start of SESSION-2026-09-28-viper.md section 4 runs on Caltech, **both arms (1x and
10x), 2 H200 each**. Viper: do NOT launch w121prod 1x/10x fresh starts.

- Commit: rt-integration 6faba555 (src = 6d690e09), incremental build (build_inc.sh), Kokkos 4.6.02.
- Inputs: `/resnick/groups/carnegie_poc/jingze/w121prod_0928/in/w121prod_{1x,10x}.athinput` = the bench-2026-09-28
  inputs (= viper w121prod_0927) with nlim/ndiag reverted, Newton fixes kept (1x maxit 16; 10x maxit 24,
  dtmax 0.25), ck_impl_xstep 8 (spin-up), current IC, closed wall, top sponge on / bottom off,
  `flux_hst_floor = true`, outputs hst+log 10/rot, bin every 2 rot, rst every 0.5 rot.
- Goal: tlim 300 rot (3.304605e7 s). Expected ~40 rot/h (1x) and ~27 rot/h (10x) on 2 H200 (bench).
- Chain: build 3607758 -> smoke 3607759/3607760 (50 cycles, fail on FATAL/NOT-CONVERGED) -> prod 3607761 (1x),
  3607762 (10x). Run dirs `/resnick/groups/carnegie_poc/jingze/w121prod_0928/{1x,10x}`.
- For the nx1 256 remap at rot 300: the rst files will be copied to viper (or the remap run here).

**Update 09-29 (Caltech, NOTE-2026-09-29-dhj-rsolver.md):** both inputs switched to `rsolver = lhllc` BEFORE any run
started (build 3607758 done, 6faba555 md5 in the job logs; smokes 3607759/60 were still queued; no hllc run exists, so
nothing to stop or rename). Same job ids: smoke 3607759 (1x) / 3607760 (10x) -> prod 3607761 (1x) / 3607762 (10x).
**09-28 16:10:** smoke folded into the production jobs (50 cycles first, same allocation; exit on FATAL / NOT-CONVERGED), so the old ids are cancelled. New: **3608646 (1x), 3608647 (10x)**. **Smoke 09-28 16:45 (H100, jobs 3609784/3609785, same binary and inputs, lhllc): both arms rc 0, 0 FATAL, 0 NOT-CONVERGED, 0 NaN.** The production jobs repeat the smoke on H200 before starting.
