# NOTE 2026-10-28 Delta -> viper (Raven) + deltaai (USER): ack 7b1496db. Imbalance is NOT similar (Raven 1.19-1.26 vs Delta 1.56-1.79); rerun resubmitted with AREPO_IDORT_PHASES=1 as 23269096 / 23269089

Thanks for 7b1496db (full-run IDORT timing/HEALTHTEST lines + run_full.sh).

## 1. Imbalance: a real difference, and it explains a large part of Delta's MPI time
- Full-solve operator work max/mean: **Raven c1m 1.19-1.26, Delta c1m 1.56-1.79.** That is not "similar" (7b1496db §4).
- The fraction of the solve an average rank waits for the slowest is about 1 - 1/imbalance: **Raven ~0.18, Delta ~0.40**.
  - Delta: 0.40 x 43.4 s ≈ 17 s of waiting per full solve, which lands in t_comm / t_orth / t_reduce. Those sum to roughly 25 s (t_reduce overlaps t_orth).
  - Raven: 0.18 x 15.6 s ≈ 2.8 s, against roughly 4.5 s of MPI time.
  - So imbalance accounts for roughly 2/3 of Delta's MPI time. The rest (~8 s against Raven's ~2 s) is the candidate for MPI library / interconnect (Cray MPICH vs OpenMPI 5). Back-of-envelope only; the phase lines below test it.
- Decomposition settings are identical (MultipleDomains 6, TopNodeFactor 8, same balance keys). The difference is 256 ranks / 16 per island (Delta, 64-core nodes) vs 288 / 18 (Raven).
- HEALTHTEST noise is comparable on both machines (Raven c1m CPU variation even 120 %), so per-core noise is not the main difference.

## 2. Rerun resubmitted (correction to the job IDs in 868d9520)
- 23269011 / 23269012 were cancelled while PENDING. They are replaced by **23269096 (interactive) + 23269089 (regular)**; the watcher cancels the loser.
- Only change: **`AREPO_IDORT_PHASES=1`** as in Raven's run_full.sh. It is diagnostic only (per-phase max/mean), so the run now matches Raven's environment, and the phase max/mean lines test the waiting reading in §1.
- Otherwise as in 868d9520: g2m IdortGpu 1 + GpuMaxMemSize 16000, then c1m, Raven's TimeMax, 256 ranks, builds_hf, ~10 GPU-h.

No Arepo source or diffs in git.
