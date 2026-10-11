# TASK raven -> DeltaAI: BSG true repro R / RC PRODUCTION to t = 6.48 d (race with Caltech + Raven; user GO 10-11)

From: Raven coordinator. User: "you can try delta ai". Caltech's H200 start is ~2.4 d out, Raven's gpu partition ~1 d out
(Raven jobs 31078702 bsgrec_R, 31078703 bsgrec_RC). **First arm to START wins for that arm**: push a STARTED line in your NOTE
as soon as link 1 of an arm runs; Raven/Caltech then cancel their copy of that arm. Do the smoke TASK
(TASK-2026-10-11-deltaai-bsgrec-smoke.md) first, in the same session: its binary, files and run decks are the ones used here.

## Runs
- Binary: the prodfinal-1011 @ 530c38e2 he_star_m1 GH200 build of the smoke TASK (report md5).
- Decks: docs/handover/bsgrec-1011/bsgR.athinput.in and bsgRC.athinput.in (@FILES@ -> your files dir), unchanged.
  Command line only: `time/tlim=5.6e5` (6.48 d) and `rad_m1/vet_gd_twin_lowmem=false` (bitwise vs lowmem per your
  NOTE-2026-10-09-deltaai-bsg-mem2; faster when memory allows). Both keys are in the deck.
- Layout per arm: **2 nodes x 4 GH200 = 8 ranks, 2 MeshBlocks per GPU**, 1 rank per GPU, MPICH_GPU_SUPPORT_ENABLED=1.
  Your 10-09 mem2 numbers (be/hr): 1 block/GPU 0.25 s/cycle, 4 blocks/GPU 0.87 s/cycle; hesdirk2 costs ~1.35x be.
  Expect ~0.5-0.7 s/cycle, ~6500 cycles -> ~1-1.3 h per arm. If 2 blocks/GPU does not fit, use lowmem=true.
- Links: <= 1 h each (ghx4-interactive), or one longer job on ghx4 if your allocation allows; chain restarts from the
  newest rst with rst_info.py last_time keys (docs/handover/bsg-hrdet-1009). rst every 1 d as in the deck.
  The two arms may run one after the other or at the same time (4 nodes total) as the queue allows.
- Stop rule: rc != 0, FATAL, NaN or NON-CONVERGED -> stop that arm, record cycle, t, last 40 lines of run.log; no reruns with
  changed keys.

## Report (NOTE-2026-10-11-deltaai-bsgrec-prod.md on bsg-files-1009, append)
- STARTED per arm (job id), then after each link: t, cycle, dt, s/cycle, Picard mean/max, FATAL/NaN/NON-CONV counts,
  `time2_pred_cap ... limited` and NOT ADMISSIBLE counts, newest rst, GPU peak.
- At the end keep the run dirs (hst, bin dumps 0.5 d, rst) and post their paths; Raven does the hr | R | RC analysis
  (L_top/L, L(r) at 6.0 d, budget residual). If you can, also run budget_bsg.py-style items 1-4 of the Caltech 10-09 NOTE.
- GPU-h used.
