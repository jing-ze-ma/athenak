# NOTE Caltech: BSG true repro, recommended implicit VET, R vs RC (2026-10-11)

## 0. ACCEPTED and QUEUING (user GO, 2026-10-11 ~07:00 PDT)
- Building he_star_m1 H200 nofma from the prodfinal-1011 tip now; the hash and md5 follow in STARTED.
- Then: smoke per arm (2 H200 each, nlim 10), and the production chains (fresh t=0, tlim 5.6e5) behind them with afterok.
- **Queue warning:** the Caltech H200 queue is long today. `sbatch --test-only` gives an estimated start of
  **2026-10-13 21:24 PDT** for 2 or 4 H200, so link 1 is unlikely before then (backfill may come sooner).
  If Raven's GPU queue frees earlier, say so here and we cancel; otherwise we run it.
- STARTED section (section 6) follows as soon as link 1 of an arm is RUNNING.
