# TASK 2026-11-19 orion -> Caltech (USER): can some staggered clean-M1 test runs move to Caltech (CPU or GPU)?

The user asks whether Caltech has machine time for part of the Orion staggered clean-M1 test work (local branches of this fork,
"stag2": the staggered clean scheme converged onto the old production physics). Our Orion / Raven / Delta queues are the bottleneck
tonight. CPU is fine; GPU (A100/H100/H200/GH200/MI300A) is better. Please reply with a NOTE on this branch.

What we would move (each item is self-contained; we ship a source tarball + input deck + restart file):
1. **HG / RSG evolved-restart arms** (the main load). Restarts `hegiant.00184.rst` (77 GB) and `rsg.01510.rst` (60 GB), 15-20 cycles each.
   Memory ~500 GB total per run (8 ranks x ~60 GB on CPU), so >= 2 CPU nodes of 256-512 GB, or 8 GPUs with >= 80 GB each.
   On 24 CPU ranks (12 Raven nodes) one run takes 25-40 min; on 8 GH200 a few minutes. About 46 runs (one-setting-off ablation, two
   code lines) = ~230 CPU node-hours, or ~22 GH200 node-hours.
2. **CPU gate rows** (small decks, 1-2 nodes, minutes each, ~70 rows per gate) — only if (1) is not possible.

Questions:
1. Which machine(s) and partitions could take these runs in the next ~24 h, and roughly how many node-hours / GPU-hours?
   Node memory (CPU) / GPU memory?
2. Can we get a login (or can you run them for us from a tarball + run script)? Is there a data path for ~140 GB of restart files
   (Globus, rsync from Orion/Raven, or regenerate from your own production restart of the same runs)?
3. Do you still have production restarts of the RSG / HG runs at or near rsg.01510 / hegiant.00184 on your side (that would avoid
   the transfer; we would also accept any later evolved restart of the same runs)?
4. Compilers / modules for AthenaK (Kokkos 4.6, CUDA 12.x or ROCm 6.x, MPI) on that machine.

Nothing urgent to change on your side; numbers and yes/no answers are enough.
