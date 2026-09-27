---
name: ck-store-0926
description: ck-store branch (8abca837, not merged) -- the storing pass runs WITHOUT the chain kernel (problem/ck_store_split default true) + lP 9->6 + scratch interleave + dup-store removal; with ck-jlin the total is ~35 % of ck (2 H200) / 38 % (1 H200); noise ratio 1.21x (at the edge); awaits the user's merge decision
metadata:
  type: project
---

Branch `ck-store`, worktree `wt_ckstore`, HEAD 8abca837 (4 commits on 5d6a68df). Work dir `/resnick/scratch/jingze/ckstore/`.

**Kept commits**
- 8fdeaef7: lP 9 -> 6 slots. Bitwise.
- 8414ebba: scratch interleaved by radial index. Bitwise.
- 5ad63742: drop pass 2's duplicate kappa-rho store. Bitwise. Those three together: -0.8 % / -1.1 %.
- 8abca837: `problem/ck_store_split = true`. `ck_coef` stores the operator with the chain kernel's own expressions (the operator is bitwise); the linear kernels produce Src/Fb. -7.0 % (2 H200), -8.5 % (1 H200). Round-off.

**Final same-job timing** (hpc-sm-02-15)

| setup | pre-jlin | jlin | final | cad-only |
|---|---|---|---|---|
| 2 H200 | 32.91 | 28.46 | 26.52 (-19.4 %) | 14.83 |
| 1 H200 | 59.34 | 49.78 | 45.55 (-23.2 %) | 23.39 |

- ck cost is defined as pre-jlin minus cad-only: 18.08 ms (2 GPUs) and 35.95 ms (1 GPU).
- By that definition: ck-jlin 24.6 %, ck-store +10.7 %, total 35 % (2 GPUs) / 38 % (1 GPU).

**Gates**
- CPU: split=false is bitwise; split on gives hst 1e-14.
- GPU 1 vs 2: bitwise.
- **Noise: 1.21x the mean of the noise pairs (range 1.11-1.25x), at the ~1.2 edge.** Mean T drift -3e-6 K.

**Dropped**: compile-time specialisation (more registers), ck_coef with a thread per (pair, cell), and the lin_build prologue split (no gain).

**Remaining ck cost per storing call (1 H200):**

| kernel | ms |
|---|---|
| ck_coef (throughput-bound) | 12.3 |
| ck_beam_tau | 22 |
| ck_lin_build | 22 |
| jlin | 14 |

Per lin pass: lin1p 15 ms, lin_sum 4 ms. HIP not gated.

**MERGED 09-26 (user: "merge it now pursue others"): 55762ada.** Next targets: ck_beam_tau, ck_lin_build, ck_coef, lin1p, lin_sum-into-rt_apply.
