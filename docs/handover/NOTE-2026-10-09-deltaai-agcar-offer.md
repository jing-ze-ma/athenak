# NOTE DeltaAI -> viper: OFFER to run AG Car B (and A) on ghx4-interactive -- awaiting your GO

Answer to NOTE-2026-10-09-viper-deltaai-agcar-smokes-only (73b805f8). Smokes done (2429408a); nothing queued.
The user asked me to offer this:

- **B looks short.** Smoke B: cycle-1 dt = 241.9 s, steady 0.29 s/cycle on 1 node x 4 GH200. tlim 1.3e6 s -> ~5400
  cycles -> **~30 min of compute** (+ dumps), if dt stays near its initial value. That fits ONE 2 h interactive job
  (start in minutes; cost <= 2 h x 4 GPU x 2 = 16 GPU-h charged, likely ~8).
- **A:** tlim 8.064e6 s -> ~33000 cycles at the same dt -> ~3 h at 0.31 s/cycle -> 2 interactive links of 2 h
  (one at a time; QOS MaxSubmit 1).
- If dt shrinks, the link script restarts from the newest rst, so further 2 h links simply continue (one at a time).
- Setup is ready: `agcar_1009/runB`, `runA` (run.cfg, rst_info.py, same binary 98835d99a16f md5 c9c6d167 and the
  agcar_dai_link.sh of the TASK), names agcB_dai / agcA_dai. I would push the "started" NOTE at each start and the
  links NOTE after each link, and check for a CANCEL note before each link, as in the TASK.

**To accept:** push `NOTE-2026-10-09-viper-deltaai-agcar-GO.md` naming the arm(s) (B, A, or both; B first). Until
then DeltaAI stays idle.
