# NOTE 2026-11-21fu orion -> caltech (USER rule via coordinator): **production runs NEVER share their node**

1. **Production nodes are exclusive:** BSG on hpc-sm-01-17, T9L6.2 on hpc-sm-02-11, m35mms on hpc-sm-01-11, and any future production node. **No test, build or analysis step may run there.** This includes holder steps doing python or seam/profile work, and srun into those nodes' holders. Login-node read-only monitoring (cp / tail / grep / ls / git) is fine.
2. **Tests only on 02-03, 02-04, 02-09, 02-10, 02-16**, with **one GPU job per node** unless we explicitly ask for packing. Builds and CPU analysis go on login nodes or one of these test nodes, never on a production node.
3. **Please confirm now, per production node (01-17, 02-11, 01-11):** `squeue -w <node>` (or the holder's step list) and `nvidia-smi` showing only that star's 4 ranks.
4. **Fix the 02-10 double-booking** that caused H2A_rsg's CUDA OOM (143 GB, NOTE ft). Make sure your runner launches the next arm on a node only after the previous step on that node has exited, e.g. by checking nvidia-smi for 0 MiB before launch. Tell us which arms overlapped on 02-10, so we can mark any affected timings (H2S_hg_a?) as invalid.

-- orion (coordinator, for the USER)
