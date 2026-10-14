# TASK Caltech (DR, e1div-1014 part 3): dt cut-and-retry (fin4-dtretry) -- build, exactness gate, and the E1 infall target run

USER GO (via coordinator): a general dt cut-and-retry for the implicit scheme, added to FIN4. Branch **origin/fin4-dtretry @ 5b7e1bc5**
(= FIN4 d505d21f + dt_retry). Analysis/validation only: do NOT touch e1_1013/prod, its rst, or any production chain.
Raven results so far: never-triggered runs BITWISE vs FIN4 (full E1, BSG P, 4 IDORT decks); forced retries with factor 1.0 BITWISE
(q1D E1, reduced 3-D E1 on 4 A100 with CUDA-IPC halos, BSG, IDORT); restart with a live dt cap BITWISE; a Picard-limited E1 that FIN4
kills at cycle 0 survives (3 cuts, recovered). This TASK adds the full-size H200 gate and the real target: the E1 infall from rst 00022.

## What dt_retry does (keys in <time>; nothing to set for a fresh run, it is on by default; a restart of an OLD rst needs time/dt_retry=true)
- Every cycle is checkpointed in memory; a cycle whose Picard solve ends NON-CONVERGED (stage or backward Euler), DIVERGED, or a BE solve
  writing eint <= 0, is put back exactly and taken again with dt x 0.5 (up to 8 cuts, then FATAL as before); the cut dt is a cap that grows
  x1.25 per cycle after 10 good cycles and is dropped once the CFL dt is below it. Every event is logged: `DT_RETRY cycle=... t=... dt A -> B
  reason=... resid=... cut n`, `DT_RETRY recovered ...`. End of run: `<time> dt_retry: N checkpoints (s, MB), M retries, cap`.

## 1. Build (same toolchain/flags as your d505d21f build; incremental from that build dir is fine)
`git fetch; git checkout 5b7e1bc5` -> `athena_gpu_he_star_m1_5b7e1bc5_nofma` (PROBLEM=he_star_m1, MPI, CUDA HOPPER90, -ffp-contract=off). Report md5.

## 2. Gate on the production E1 deck (e1_1013/agcA_E1_final_wbrmax.athinput d29e4615), 1 node x 2 H200, your E1 launch, FRESH, nlim 100
Keys for all arms: `time/nlim=100 time/ndiag=1 output2/dt=1.5e4 output3/dt=1.5e4 output4/dt=1e30 output5/dt=1e30`.
New keys cannot be added on the command line of a fresh run: make deck copies with lines inserted after `<time>`:
`awk '{print} /^<time>/{print "dt_retry_force_cycle = 22"; print "dt_retry_force_fac = 1.0"}' deck > deck_F22` (and F88 with 88).
- R: d505d21f binary, deck as is.   N: 5b7e1bc5 binary, deck as is.   F22, F88: 5b7e1bc5 binary, deck_F22 / deck_F88.
- Compare (python3 docs/handover/e1div-1014/cmpruns.py A B: hst rows + every bin dump's data): R vs N, R vs F22, R vs F88 -> expect
  `CMPRUNS BITWISE` each. F22/F88 logs must show one `DT_RETRY ... reason=forced` line each.
- Peak GPU memory per arm (nvidia-smi --query-gpu=index,memory.used -lms 500 during each run): R vs N (expected + ~2.3 GB/rank).
- s/cycle R vs N (cycles 10-100).

## 3. Target: the E1 infall from rst 00022 with dt_retry (dir e1div_rst/D1, rst 00022 hard-linked as for R0/R1)
`e1launch.sh <5b7e1bc5 bin> -r rst/agcar3d.00022.rst <exactly R0's keys> time/dt_retry=true time/tlim=1.30e7 -d e1div_rst/D1 -t 05:40:00`
(the restart of an old rst takes FIN3/dad00933 methods (FIN4 keys off without them in the file), so D1 differs from R0 only by dt_retry).
Report: rc; end t/cycle; all `^DT_RETRY` lines (count, first 20, last 20); FATAL/DIVERGED; NON-CONVERGED count; Picard mean/max;
dt every 100 cycles (`grep "cycle=" run.log`); s/cycle; the end-of-run dt_retry line; e1div_prof.py on D1's bins (push prof_* only).
Push to docs/handover/e1div-1014/caltech_out_dtretry/ + NOTE-2026-10-14-caltech-e1div-dtretry.md. If R0/R1 (part 2) finished, their reports too.
