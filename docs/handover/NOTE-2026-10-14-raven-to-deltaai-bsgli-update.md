# NOTE Raven -> DeltaAI: bsgli-val UPDATE (supersedes the hash in TASK-2026-10-14-deltaai-bsgli-val.md)

1. Build **8e018272** from the SNAPSHOT branch **fin4-linhao-bc-8e018272** (NOT the moving branch fin4-linhao-bc, now ahead at a diagnostic-only commit) instead of eb1e032c. Name `athena_hes_gpu_8e018272`.
   (Physics bitwise to 37798491/eb1e032c; only the L_in history diagnostic changed.) Run the same val.sh gate with it.
2. In your NOTE also report your **current ghx4 BATCH start estimate** for a 1 node x 4 GH200, 24 h job
   (`sbatch --test-only` or squeue --start), and the interactive situation.
3. BACKUP PRODUCTION (user-approved, conditional): ONLY IF the gate is BITWISE **and** your batch estimate is earlier than
   ~15 h from now, you may queue a backup production chain of the same run: same binary/deck/layout as the gate,
   fresh start, TLIM 4.96e6, 24 h batch links (expect ~13-16 h of compute at your 0.80 s/cycle x FIN4 +14 %),
   restart from the newest valid rst, prune keep 3, STOP on rc!=0 / FATAL / NaN / DIVERGED / NON-CONVERGED > 0.
   Do NOT use chained 1-h interactive links for this. If the estimate is later than ~15 h: do NOT queue it; just report.
4. Report job ids and the estimate in NOTE-2026-10-14-deltaai-bsgli-val.md.
