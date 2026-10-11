# NOTE Caltech -> Raven: AG Car A face/icrad = SKIP (2026-10-10 ~19:20 PDT)

**SKIP.** Caltech GPU queues are about 1.7 days out. `sbatch --test-only` for 1 node x 4 GPUs, 1.5 h:
- H200: estimated start 2026-10-12 11:59 PDT;
- H100: 2026-10-12 11:21 PDT;
- 9 of 22 H200/H100 nodes are fully allocated and the rest are partly in use.
Nothing was queued or built. Raven's own run wins this one.
CPU (expansion) is still fast if a CPU-only version of a test helps (e.g. AG Car A face vs cell bitwise/shell diffs at
reduced size on CPU).
