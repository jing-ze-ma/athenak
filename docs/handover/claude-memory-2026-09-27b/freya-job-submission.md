---
name: freya-job-submission
description: Preferred SLURM partition and submit settings for AthenaK runs on Freya
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 0df72474-f27b-49d2-9277-124d645c5193
  modified: 2026-08-05T04:23:47.002Z
---

Submit AthenaK jobs to the **p.shared** partition, not p.exclusive — it schedules faster
(shorter/no queue).

**Why:** user asked for it; shared nodes are usually idle while p.exclusive queues.
**How to apply:** in submit script use `#SBATCH --partition=p.shared` AND add
`#SBATCH --cpus-per-task=4` so SLURM reserves 4 cores per MPI task (16 tasks x 4 OMP threads
= 64 cores). On p.exclusive the whole node is allocated so cpus-per-task is unnecessary, but
on p.shared without it the OpenMP threads oversubscribe one core per task. p.shared nodes have
~224 cores / 512 GB, so 64 cores + 50 GB fits fine. Don't move an already-running job.

**User preference (firm, 2026-07-21):** after editing + building, JUST SUBMIT the real job.
Do NOT run pre-flight sanity checks first (no `-m` mesh dry-run, no short `mpirun nlim=N` smoke
test) — the user interrupted both and said "just run it without checking." Build, archive the
previous run's outputs, recreate an empty `bin/`, then `sbatch` directly; analyze partial output
from the real run instead. (Login-node pre-checks also hit spurious issues: a cosmetic segfault
in MeshBlockPack destructor on `-m` exit, and MPI-IO needing a `bin/` subdir.)

**EXCEPTION (2026-08-05):** p.shared is only faster when it is not already full of the user's
OWN jobs. A large `kt` job ARRAY (185145_[235-999], ~765 tasks pending, 4 running) had the
partition saturated, and a new p.shared job sat PD with reason "Priority". Moving it to
**p.exclusive** started it in ~2 seconds. So: run `squeue -u jinma` first when a job sits PD,
and switch partitions if the user's own array is hogging p.shared. When moving to
p.exclusive, drop `--mem` (whole node is allocated).

**DO NOT also drop `--cpus-per-task` when OMP_NUM_THREADS > 1.** This note used to say to
drop it; that is WRONG and cost an hour on 2026-08-17. With `--ntasks=16` and no
`--cpus-per-task`, SLURM still binds ONE core per task, so 16 tasks x 7 OpenMP threads fight
over 16 cores. The symptom is nasty because it looks like a physics or code problem: the run
burns CPU, reports no errors, keeps a healthy dt, and simply crawls. It read as a "140x
slowdown in the general EOS"; the real figure measured properly is 4.3x. Keep
`--cpus-per-task` on BOTH partitions whenever threads are used.

**AN ORION NODE HAS 112 PHYSICAL CORES, NOT 224.** `scontrol show node` reports
`CPUTot=224`, but the topology is `Sockets=2 CoresPerSocket=56 ThreadsPerCore=2` -- the 224
are HYPERTHREADS. So **16 ranks x 7 threads = 112 already saturates a node**; there is
nothing left idle. Do not read `NumCPUs=224` in `scontrol show job` as "half the node is
wasted" -- that mistake was made on 2026-08-17 and led to cancelling four healthy long runs
to "fix" it. Asking for 28 x 8 = 224 is REJECTED outright ("Requested node configuration is
not available") because SLURM allocates whole cores, and only 112 exist.

Sanity check before resizing anything: `sinfo -p <part> -o "%.10P %.8z %.6c"` prints
`S:C:T` -- read cores, not CPUs.

**p.exclusive is not automatically the escape hatch.** On 2026-08-17 p.shared had 4 IDLE
nodes while p.exclusive was 98/98 allocated, so jobs moved there sat PD longer. Check
`sinfo -p <part> -o "%.6D %.6t"` on both before switching. Also: p.shared nodes showing
`alloc` still have free cores (they are SHARED, ~224 each), so a smaller job
(`--ntasks=8 --cpus-per-task=4`) usually starts immediately alongside your own big ones --
that is how the 200-cycle cost tests got through while three full runs held three nodes.

Related: [[solar-convection-test]], [[freya-build-procedure]]. Incremental `make -j` (no build/
wipe) is correct after a source-only edit; only wipe build/ when CMake config changes.

## 16x7 MEASURED, not just assumed (2026-08-24)

The 16x7 convention above was never actually compared against anything. It has now been
scanned on one orion node: `deep_hot_jupiter_rt`, general EOS + correlated-k, mesh
64x64x128, meshblock 64x8x8 (128 blocks), 200 cycles, zone-cycles/s:

| ranks x threads | 1x112 | 2x56 | 4x28 | **8x14** | **16x7** | 28x4 | 56x2 | 112x1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| zone-cycles/s | 3.86e6 | 4.25e6 | 4.30e6 | **4.52e6** | 4.49e6 | 4.27e6 | 3.71e6 | 3.20e6 |

**Broad plateau over 4-28 ranks; 8x14 and 16x7 are tied to within 0.5 % (noise).** So
16x7 is right and there is nothing to gain by retuning it. The falloff is at the ends:
pure OpenMP (1x112) loses 14.5 % because one rank straddles both sockets (2 x 56), and
pure MPI (112x1) loses 29 % to boundary exchange.

**The decomposition, not the flags, is what wastes a node.** A meshblock of 64x64x64 on a
64x64x128 mesh gives just TWO meshblocks, which caps the run at 2 ranks = 56 cores no
matter what is requested. That is a GPU-shaped split (right for 2x MI300A on viper, wrong
for a CPU node). On CPU use meshblock 8x8 -> 128 blocks. Recorded in the
`<meshblock>` comment of `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`.

Note that input's mesh is now nx2 = **64**, not the 56 of the older CPU runs, so the
"legal nx2 are divisors of 56" trap is gone -- every power of two up to 64 works. Dropping
back to 56 would be a RESOLUTION change (14.3 % fewer cells), not a decomposition one.
