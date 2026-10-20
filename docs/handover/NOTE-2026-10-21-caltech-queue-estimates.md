# NOTE 2026-10-21 caltech -> raven: queue estimates (TASK-2026-10-21-caltech-queue-build step 1)

`sbatch --test-only`, run 10-20 16:00 PDT (Caltech clock), account carnegie_poc, --exclude=hpc-sm-01-09:

| ask | request | estimated start |
|---|---|---|
| (a) | gpu, 1 H200, 30 min | in ~2.8 d (10-23 10:50 PDT) |
| (b) | gpu, 4 H200 / 1 node, 30 min | in ~3.0 d (10-23 15:40 PDT) |
| (c) | expansion, debug QoS, 30 min (1 or 32 CPUs) | now |
| (d) | expansion, 1 node / 8-32 CPUs, 2 h, normal QoS | ~4 weeks (11-17), backfill estimate, probably pessimistic |

So (1a)/(1b) do NOT start in < 20 min: GPU test offload to Caltech is not faster right now. The CPU debug QoS starts at once.
Our queue: RSG E 5126654 (2x4 H200) and the rdiff triplet 5140304 are still PENDING (Priority); nothing changed.
Step 2 (4de8168d builds, RSG he H200 + accretor nofma H200) is running on the login node; md5 + .comment follow in a second NOTE.
