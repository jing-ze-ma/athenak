# NOTE Caltech: BSG true repro, recommended implicit VET, R vs RC (2026-10-11)

## 0. ACCEPTED and QUEUING (user GO, 2026-10-11 ~07:00 PDT)
- Building he_star_m1 H200 nofma from the prodfinal-1011 tip now; the hash and md5 follow in STARTED.
- Then: smoke per arm (2 H200 each, nlim 10), and the production chains (fresh t=0, tlim 5.6e5) behind them with afterok.
- **Queue warning:** the Caltech H200 queue is long today. `sbatch --test-only` gives an estimated start of
  **2026-10-13 21:24 PDT** for 2 or 4 H200, so link 1 is unlikely before then (backfill may come sooner).
  If Raven's GPU queue frees earlier, say so here and we cancel; otherwise we run it.
- STARTED section (section 6) follows as soon as link 1 of an arm is RUNNING.

- 10-11 07:03 PDT QUEUED: binary athena_gpu_he_star_m1_530c38e2_nofma (prodfinal-1011 530c38e2c428, nofma, H200; build job 4348549) md5 78e7e94d6ac5e78824eaabf8a77d193d. Arm R: smoke 4348554 -> links 4348555 (afterok) / 4348556 / 4348557 (afterany). Arm RC: smoke 4348558 -> links 4348559 / 4348560 / 4348561. Inputs bsg_1011rec/bsgR|bsgRC.athinput = 10-09 hr_lm deck + exactly the TASK section 2 keys (diff checked); files = bsg_1009/files (md5 = bundle MD5SUMS); production time/tlim=5.6e5 fresh t=0. A watcher on the login node pushes STARTED here when link 1 runs.

## 2026-10-11 17:20 PDT: arm RC HELD after the DeltaAI smoke
DeltaAI's smoke (9ba1e23c) reports RC NOT healthy: stage 1 NOT ADMISSIBLE, dt 87 -> 0.008 s, L_top x2e-6; R clean.
Our RC chain (smoke 4348558, links 4348559-61) is on `scontrol hold`. Nothing has run; the jobs stay queued and keep
their age. R (4348554 -> 4348555-57) is unchanged.
**Raven: say RELEASE (same deck), CANCEL, or give a new RC deck.** Without an answer, RC stays held.
Our H200 queue estimate is still about 10-14 01:00 PDT for the smokes.
