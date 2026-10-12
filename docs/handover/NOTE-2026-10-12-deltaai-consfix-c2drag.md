# NOTE DeltaAI: conservative2 + implicit_cs2_hll = drag, full BSG + AG Car A (TASK-2026-10-12-deltaai-consfix-c2drag)

**Picked up** 20:50 CDT. consfix-1011 **8cba89a6** built: `athena_hes_gpu_8cba89a695bf`, md5 ec449854a7c97c084601499aee211e37.
Decks resolved (BSG files = bsgrec smoke files; AG Car A files md5 OK), job.sub as provided (dry run: 5 runs).
**Queued behind** the agcarA-ramp link 1 (job 3363376, running since 20:52; ghx4-interactive allows 1 submitted job per user):
submitted the moment that job ends (expected by ~21:50 CDT). STARTED line follows.

- 21:50 CDT: job **3363531** PENDING, Slurm start estimate **23:41 CDT** (~2 h, at the TASK's SKIP limit). Kept queued unless you post WITHDRAWN. The agcar-reatm job follows it (1 job per user on ghx4-interactive).

- 22:00 CDT: to start sooner, consfix (3363531, cancelled while pending) and agcar-reatm now run in ONE 30-min job **3363656** (reatm first, then consfix; scripts unchanged except the stop-launching time). Slurm estimate 23:56 CDT; a 30-min job may backfill earlier. STARTED follows.
