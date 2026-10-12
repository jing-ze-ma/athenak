# NOTE DeltaAI: full-size BSG timing of bsgaccel + conservative2 health check (TASK-2026-10-12-deltaai-bsgaccel-fullbsg)

**Picked up** 2026-10-11 19:50 CDT. Building NEW 0197f020 then BASE 422c2af7. Deck = bsgR (bsgrec smoke) + the 6 named-default
`<rad_m1>` keys. Order in one job: **C4, C4g first**, then N4, N4n, N2 (NEW), then B4, B2 (BASE). Job id + STARTED line follow.

- 19:57 CDT: NEW built (athena_hes_gpu_0197f0202e10, md5 d36ed4ff9aed3bc092ab03e42d425c93); BASE 422c2af7 building (job waits up to 10 min for it after the NEW runs).
  Job **3363179** PENDING, priority 2173, Slurm start estimate **21:53 CDT (~2 h)**: at the TASK's SKIP limit; I keep it queued unless you post WITHDRAWN.

**STARTED**: job 3363179 RUNNING on gh152 since 20:09:04 CDT (C4 first). Raven: no need to run this.

## URGENT, C4 (conservative2, default precond, lowmem): NOT healthy on the full BSG (posted 20:13 CDT, job still running)
rc 0, FATAL 0, but **NON-CONVERGED = 6**, Picard mean 34.9 / max 200 (48 solves), and hesdirk2 stage 1 is **NOT ADMISSIBLE on every
cycle 2-23** (`stage steps=1 backward-Euler steps=23 stage fallbacks=22`); vimp positivity fallbacks 41, **min E from the solve
-4.05e6**; gas Newton fallbacks 4.5e8 (1.6e-2 per cell-pass; R4 matched 1.2e-3).
**dt collapses** at cycle 3 from 87.37 s to 0.00836 s, doubles back to 2.14 s by cycle 11, collapses again to 0.0015 s at cycle 12
and to 0.019 s at cycle 19 (t only reaches 266.8 s in 24 cycles vs 2097 s at full dt).
L_top: 1.6716e38 (c1) 4.457e37 (c2) 1.854e37 (c3), then 1.77-2.69e37 (R4 matched: 4.32-4.75e37 at cycles 3-10). Unlike the old
`conservative` (RC4: L_top -> 8.7e31) it does not lose the luminosity, but it does not take a single hesdirk2 step after cycle 1.
Identical to R4 / RC4 through cycle 1 (L_top 1.67165e38). Peak 49 789 MiB per GH200 (4 blocks/GPU). C4g and the timing runs follow.
