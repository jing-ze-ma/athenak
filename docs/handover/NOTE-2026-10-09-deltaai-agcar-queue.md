# NOTE DeltaAI: AG Car queue facts (answer to TASK-2026-10-09-deltaai-agcar, step 1)

Checked 2026-10-09 12:12 CDT (19:12 CEST). Account bivj-dtai-gh: 905 of 1001 GPU-h remaining.

| partition | MaxTime | MaxNodes | DefaultTime | 1 node x 4 GH200 --exclusive, start estimate (sbatch --test-only) |
|---|---|---|---|---|
| ghx4-interactive | 2 h | 2 | 30 min | --time 2 h: **in ~4 min** (2026-10-09 12:15 CDT); 20 min: same |
| ghx4 | 48 h | unlimited | 30 min | --time 48 h: **2026-10-14 11:41 CDT** (~5 days); --time 24 h: same |

- ghx4 (regular) has a deep queue: a 1-node job starts in ~5 days whatever its length (24 h or 48 h).
  Interactive starts in minutes but caps at 2 h per job, and the QOS allows 1 job per user (MaxSubmit 1), so no
  chained interactive links can wait in the queue.
- Interactive is charged at 2x (as seen 10-01).
- Smokes A and B (step 2) go next on ghx4-interactive (20 min each); results follow in this NOTE.

## Smokes (step 2): both PASS
ghx4-interactive, 1 node gh069, `agcar_dai_link.sh <D> 10` (fresh start, 10 cycles), binary 98835d99a16f md5 c9c6d167
(verified by the link script); module stack = the build stack (cray-mpich 9.0.1, cudatoolkit 25.5_12.9, PrgEnv-gnu 8.6.0).

| arm | job | rc | FATAL | NaN (run.log, hst) | NON-CONVERGED (last) | steady s/cycle | cycle 1 dt | max GPU mem per GH200 | job elapsed |
|---|---|---|---|---|---|---|---|---|---|
| smokeB | 3348878 | 0 | 0 | 0, 0 | 0 | **0.29** (cycles 3-9; 0.54 c1, 1.09 c10 = dump) | 241.9 s | 19205 MiB | 29 s |
| smokeA | 3349002 | 0 | 0 | 0, 0 | 0 | **0.31** (cycles 3-9; 0.53 c1, 1.07 c10 = dump) | -- | 21037 MiB | 29 s |

At the cycle-1 dt (242 s), B's tlim 1.3e6 s is ~5400 cycles, i.e. ~30 min of compute, unless dt shrinks.
