# NOTE 2026-11-27zo orion -> Caltech: queue priority + pack the GPUs (USER: use idle GPUs, priority by what holds up others)

Order for the open Orion NOTEs (finish higher first; run lower ones on any GPU left free):
1. 27zj restart-key check (merge of the trunk waits on it; short).
2. 27zl combined star check (starmix): gates the fix branches. Arms of 27zf that 27zl repeats on the mixed binary can be dropped if not yet run; keep the 27zf results already made.
3. 27zm strong-scaling breakdown and 27zn roofline profile (baselines for the ray-kernel and VET-solve work).
4. SCPR 27zh/27zk (PR #798 timing + emitting sphere), HLLRES W2/W1R2 (keep running where they are).
Packing: a node has 4 H200. Single-GPU arms that are NOT timing measurements (correctness, restart, bitwise, residual arms) should run 4 at a time per node (CUDA_VISIBLE_DEVICES per step), not one per node. Timing arms (27zm, 27zn, cost pairs) keep a node to themselves. Report idle GPUs in the node table as GPU counts.
