# NOTE DeltaAI: conservative2 + implicit_cs2_hll = drag, full BSG + AG Car A (TASK-2026-10-12-deltaai-consfix-c2drag) -- DONE: OK on both

**Picked up** 20:50 CDT. consfix-1011 **8cba89a6** built: `athena_hes_gpu_8cba89a695bf`, md5 ec449854a7c97c084601499aee211e37.
Decks resolved (BSG files = bsgrec smoke files; AG Car A files md5 OK), job.sub as provided (dry run: 5 runs).
**Queued behind** the agcarA-ramp link 1 (job 3363376, running since 20:52; ghx4-interactive allows 1 submitted job per user):
submitted the moment that job ends (expected by ~21:50 CDT). STARTED line follows.

- 21:50 CDT: job **3363531** PENDING, Slurm start estimate **23:41 CDT** (~2 h, at the TASK's SKIP limit). Kept queued unless you post WITHDRAWN. The agcar-reatm job follows it (1 job per user on ghx4-interactive).

- 22:00 CDT: to start sooner, consfix (3363531, cancelled while pending) and agcar-reatm now run in ONE 30-min job **3363656** (reatm first, then consfix; scripts unchanged except the stop-launching time). Slurm estimate 23:56 CDT; a 30-min job may backfill earlier. STARTED follows.

- 22:06 CDT: reatm withdrawn, so job 3363656 runs only consfix; its wall limit lowered to 20 min (estimate still 23:56 CDT).

**STARTED**: job 3363656 RUNNING on gh056 since 22:31:54 CDT (consfix runs only).

## RESULT (posted ~22:40 CDT): job 3363656, gh056, consfix part 22:31:59-22:34:35, 5/5 runs rc 0, FATAL 0, NaN 0
**Verdict: conservative2 (implicit_cs2_hll = drag, 8cba89a6) OK on the full-size BSG and OK on full-size AG Car A.**
- **BSG:** no NOT ADMISSIBLE (0 in bC2 and bC2g; the old conservative and conservative2@0197f020 had it on every cycle from 2);
  stage fallbacks 0; dt path identical to bM (87.37 -> 87.33 s over 30 cycles); NON-CONV 0; min E 7.76 (bC2) / 7.87 (bC2g) as matched.
  vs bM at cycles 10 / 20 / 30: L_top 1.3e-9 / 4.4e-10 / 7.5e-11, E_rad 1.2e-11 / 1.9e-11 / 1.8e-11, Etot 1.1e-14 / 1.3e-14 / 1.4e-14
  (bC2g: 9.5e-10 / 4.3e-10 / 7.8e-11, 1.2e-11 / 1.9e-11 / 1.8e-11, 1.1e-14 / 2.1e-14 / 1.2e-13).
- **AG Car A:** aC2 vs aM, max relative difference per user.hst column over the 61 rows: <= 2.5e-10 for every luminosity and
  energy column (L_bot 2.5e-10, L_top 1.4e-11, E_rad 4.2e-12, Etot 2.9e-16); KE_int 2.4e-7, KEr_int 2.5e-7, Mr_int 3.6e-7
  (columns that start at 0: normalised by max|aM|); **v1sq_wall 3.7e-6** (also normalised; a near-zero wall diagnostic) is the
  only column above the expected 1e-6. NOT ADMISSIBLE 0, NON-CONV 0, Picard identical to aM.
- No DBG_T2_ADMISS blocks were printed (no stage was ever not admissible).
- No `implicit_gas_src` banner is printed by 8cba89a6; conservative2 was active (user.hst differs from matched in every run pair).

| run | Picard mean / max (solves) | stage steps / BE / fallbacks | NOT ADM | gas Newton fallbacks per cell-pass | s/cycle steady mean (median) | setup |
|---|---|---|---|---|---|---|
| bM | 2.81 / 7 (59) | 29 / 1 / 0 | 0 | 9.12e-3 | 0.613 (0.562) | 2.4 s |
| bC2 | 2.81 / 7 (59) | 29 / 1 / 0 | 0 | 9.12e-3 | 0.630 (0.579) | 2.3 s |
| bC2g | 2.58 / 7 (59) | 29 / 1 / 0 | 0 | 1.00e-2 | **0.421** (0.383) | 2.3 s |
| aM | 2.63 / 6 (127) | 59 / 1 / 0 | 0 | 3.43e-3 | 0.186 (0.177) | 19.8 s |
| aC2 | 2.63 / 6 (127) | 59 / 1 / 0 | 0 | 3.43e-3 | 0.193 (0.181) | 20.3 s |
(steady = cycles > 2, not multiples of 10, not the last.) conservative2 costs +3-4 % per cycle vs matched (BSG 0.630 vs 0.613, A 0.193 vs 0.186).

L_top per cycle (all runs) and dt per cycle: deltaai-consfix-1012/analysis.txt. bC2 L_top rows = bM to 5 digits at every row.

GPU-h: 2.7 min x 4 GH200 = **0.18 GPU-h** (the withdrawn reatm part used nothing).
Files: In git docs/handover/deltaai-consfix-1012/ (analysis.txt, analyse.py, combo.sub). On DeltaAI: /work/nvme/bivj/jma20/consfix_1011/runs/3363656/.
