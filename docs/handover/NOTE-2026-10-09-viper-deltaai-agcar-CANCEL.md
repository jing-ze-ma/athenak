# NOTE viper -> DeltaAI: CANCEL AG Car B and A production (user 10-09)

**Arms: B and A.** Reason (user 10-09): DeltaAI is the short-test machine, not a production machine
(NOTE-2026-10-09-deltaai-role-short-tests, rt-integration 6e2014eb). AG Car production stays on viper / Raven /
Caltech.

Do (step 4 of TASK-2026-10-09-deltaai-agcar):
1. `touch /work/nvme/bivj/jma20/agcar_1009/runB/CANCEL /work/nvme/bivj/jma20/agcar_1009/runA/CANCEL`
   (mkdir the run dir first if it does not exist). Do not submit any more links for either arm.
2. Job **3349057** (agcB_dai link 1) and any other agcB_dai / agcA_dai job:
   - PENDING: `scancel --state=PENDING <id>` (fresh `squeue -u $USER` in the same command).
   - RUNNING: do NOT kill it. Let it run to its code wall limit so it writes its rst, then it stops (CANCEL makes
     any later link exit 0).
3. Push `NOTE-2026-10-09-deltaai-agcar-cancelled.md` on this branch: per arm, job id(s), state at cancel
   (pending / ran), and for a link that ran: t reached and cycle (`rst_info.py` on the newest rst), s/cycle
   (steady), last NON-CONVERGED count, max GPU memory, GPU-h (4 x elapsed, x2 interactive charge). Keep all run dirs
   (viper may ask for the rst).

Next for DeltaAI: TASK-2026-10-09-deltaai-prof.md (one short profiling job), after the interactive slot is free.
