# TASK viper -> Caltech: STANDING TEST RUNNER (user 10-10: "ask Caltech to get them automatically")

Run a polling loop so viper's short test tasks start without the user pasting each one.

## Loop
Every ~15 min (`/loop 15m ...` or a background poll script; one loop only):
1. `git fetch <fork remote>` and look on branch **bsg-files-1009** for NEW files `docs/handover/TASK-2026-10-1*-caltech-*.md`
   that have no matching `NOTE-<same date>-caltech-<same name>.md` yet (answered tasks have a NOTE).
2. Skip this runner file itself and anything marked `WITHDRAWN` / `HOLD` in its first lines.
3. Process new tasks in order of commit time, one at a time per task (several jobs inside one task may run in parallel):
   build what the task names (incremental, as you do now), run it, push the NOTE it asks for on the same branch.
4. Push a short "picked up" line into the task's NOTE as soon as you start (so viper sees it was taken).

## Scope and limits
- SHORT tests only: smokes, gates, batteries, timing/scaling, profiling. Each job <= ~2 h, <= 2 nodes / 8 H200 total at a time;
  CPU batteries on free CPU nodes. Anything needing more: ask in a NOTE first.
- **Do NOT start or restart long production runs** (BSG, He giant, AG Car, accretor) from this loop: all long runs are on
  HOLD until the updated radiation scheme; a production restart will come as an explicit task AND the user's word.
- Keep a ledger `docs/handover/caltech-runner-1010/LEDGER.md` (task, jobs, GPU-h, CPU-h, result) pushed with each NOTE.
- Never force-push; never push code branches unless a task says so; do not touch other sites' directories.
- If a task fails to build or a job FATALs, report it in the NOTE and continue with the next task (no code fixes on Caltech
  unless the task asks for them).

## Already pending (process these first if not done)
- TASK-2026-10-10-caltech-steep-battery (done? its NOTE says complete -> skip)
- TASK-2026-10-10-caltech-vgdspeed (H200 timing; may be superseded by a vgdfuse-scaling task -> run the newer one)
- any TASK-2026-10-10-caltech-* pushed after this file.
