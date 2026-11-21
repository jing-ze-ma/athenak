# NOTE 2026-11-21jk orion -> caltech: graceful stop for v20.7 (operator note)

How to stop an AthenaK production run cleanly. Works in v20.6 989dfd32b and v20.7 af9e8d420
(code 5ba73ecab, your B_v207) and later. Old binaries (v20.5 6963db013, B_j3 a5257b7ee and
earlier) ignore all of this.

1. Stop now, from any shell: `touch <run dir>/STOP` (or `scancel --signal=USR1 <jobid>`).
   The run finishes the current cycle, writes `rst/<basename>.NNNNN.rst` through the normal
   restart output, prints
   `stop requested (STOP file|SIGUSR1|SIGTERM) at cycle N: restart written rst/...`
   and `Terminating on stop request`, and exits 0. Rank 0 removes the STOP file, so the next
   link restarts as usual. The restart is written even with `<job>/final_restart = false`.
2. Before the batch limit: add `#SBATCH --signal=USR1@900` (NO `B:`; `B:` signals only the
   batch shell, not the srun tasks). Pick the lead time > one cycle + one restart dump
   (~15 min is safe). Or pass `-t hh:mm:ss` to athena (wall time measured from the start;
   same stop path, message "wall time").
3. Do NOT stop with plain `scancel` / `--signal=KILL`: scancel sends TERM and then KILL after
   KillWait (~30 s), too short for a large restart dump. SIGTERM is handled like USR1, but
   rely on it only when the dump fits in KillWait. Wait for the job to exit by itself.
4. Ranks started by mpirun (not srun) live in the batch step: use
   `scancel --full --signal=USR1` (mpirun forwards USR1 to the ranks), or simply the STOP
   file, which works for every launcher.

The check runs once per accepted cycle (one MPI MAX reduction), so all ranks stop after the
same cycle; a signal or STOP file during a long cycle takes effect at its end.

-- Orion (LEAN6, for the USER)
