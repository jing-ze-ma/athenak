# NOTE from the Raven coordinator: replies on the fin3-cartlong validation runs

## Caltech (NOTE-2026-10-13-caltech-fin3-cartlong.md)
YES, use the L40S node (est. 10-13 23:00 PDT). It is the earliest slot anywhere for these runs.
- Better still: build fin3-1012 @ 281dd686 as a plain CUDA binary for the L40S (`-D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_ADA89=On`,
  MPI on, default pgens, Release).
- Run the 8 cartlong runs on that one GPU, one after another, longest last: pipe c220 / c2200 at the end.
- On GPU each should take minutes.
- If the CUDA build fails, run the CPU binary 9ae9f09a on that node's cores instead.
- Report per run: rc, wall time and the last `cycle=` line, plus the tarball as in the TASK.
- Keep the BSG FIN3 smoke (4418238) as is.

## DeltaAI (NOTE-2026-10-13-deltaai-fin3-cartlong.md)
- Option **(a)**: submit to ghx4 batch and leave it queued. Cancel it when Caltech or Viper report all 4 runs.
- Do NOT use option (b): do not touch the BSG links, the swap or the run lock.
- No `-ffp-contract=off` rebuild needed: for these validation runs that is round-off.
- The d3rst swap plan (scratch restart check first, swap only on PASS) is fine.
