# NOTE 2026-10-26 caltech -> viper (USER): reply to 21fdedef: Caltech CAN run g25_full on 2 H200 nodes, but NOT soon (all our nodes are in production, queue is long); same blocker as DeltaAI (e93504ccef missing); nothing started

**Short answer.**
- **Hardware: yes.** 2 H200 nodes (8 H200, 256 cores, ~1.5 TB RAM per node) fit easily; 4 nodes are possible.
- **Blocker: the code.** Caltech has only the **10-24 bundle** plus our own local branches (idort-gpu-caltech-opt c663a54c5, nltree-g25-caltechfix). **e93504ccef is not here** (`git cat-file` fails in all three clones under arepo2_gpu/). So the incremental bundle ce74114d cannot be applied. We need either a full bundle of idort-gpu-nltree-2, or an incremental one whose prerequisite we have (6cef4d1d6 idort-gpu-s1 or 1108335a7 nltree-g25).
- **Queue: the real problem.** Our 6 held H200 nodes are all running the AthenaK productions: RSG on 3 nodes, HG on 2, the accretor on 1. The gpu partition has ~550 pending jobs. Our pending 1-node holder 5274196 has an estimated start in ~7 h. A new 2-node job would probably wait a day or more; `sbatch --test-only` on gpu returns no estimate here.
- **The only fast route** is to pause one production (e.g. HG for ~1 h on its 2 nodes). That is the user's call; we do NOT recommend it unless Raven and DeltaAI both fall through.
- **Nothing started.**

## Nodes and layout
- 2 nodes = 256 ranks (128 per node, 32 per GPU); or 4 nodes = 512 ranks if a holder frees up.
- Memory is no issue: Raven peaked at 135 GB host / 117 GB GPU per node (4 A100); the H200 nodes have ~1.5 TB host and 4 x 141 GB GPU.
- Inputs OK: g25_inputs.tgz md5 **2348f0db65cc644446240951bd518235** (matches); ics/g25_full.hdf5 md5 40c5ce81 OK.
- Build: ~15 min in a held allocation (10-24 recipe, local HYPRE). `--exclude=hpc-sm-01-09` (bad node) and OMPI_MCA_io=romio321 for multi-node jobs.

## Expected speed (estimate, not measured on g25_full multi-node)
- 10-24 Caltech, 1 node, 32 ranks: GPU 83.6 s vs CPU-IDORT 223 s vs pure CPU 235 s per top step.
- Raven 4 nodes (288 ranks, 16 A100): GPU 20.9 s, CPU 28.5 s per top step.
- 2 Caltech nodes (256 cores, 8 H200) should be about like the Raven 4-node job: **~20-30 s GPU, ~30-40 s CPU per top step**. 3 arms x (init + 4-5 top steps) **~30 min**, so a 1 h job.
- Cost: 2 nodes x 1 h = 8 H200-h (~$15).

## Recommendation
DeltaAI (interactive start ~35 min) is the faster site once the bundle with the prerequisite is copied there. Caltech is the backup. If the user picks Caltech, copy the same bundle to `/resnick/groups/carnegie_poc/jingze/arepo2_gpu/bundle_1026/` and say whether to pause a production or wait in the queue.
