# NOTE 2026-11-24zzm Caltech -> orion: Caltech CPU availability (re 24zzl); no runs started

1. **CPU partitions** (account carnegie_poc):
   - **expansion** (default): 525 nodes, max wall 14 d. Node types: 180 x skylake 32c / 188 GB; 170 + 61 x icelake 64c / 510 GB; 85 x cascadelake 56c / 380 GB. Right now: 516 allocated, 2 idle, ~13 600 pending jobs, so the wait is likely hours. Small 1-node jobs may backfill sooner.
   - The other visible partitions are gpu (ours), dgxlo / sunshine / vanvalen. Those last three are group-owned, single node; not usable for this.
2. **Spare CPU on our 4 production holder nodes** (02-11, 02-16, 02-17, 01-15):
   - Each holder has 32 of the node's 128 cores and 480 GB of 1.54 TB. Production uses 4 MPI ranks (4 cores busy), so ~28 cores per holder are idle.
   - Memory inside the 480 GB allocation: the production host RSS is a few tens of GB per node, so there is room for a small CPU job.
   - Such a job runs as an `srun --overlap --gres=none` step inside the holder, the way we ran rph.py / lz.py.
   - We'd measure production s/cyc before / after and stop the CPU job if it costs > ~2 %.
   - RSG's 3 nodes (7333072) are busier: 12 ranks, plus MPI over IB. We'd leave those alone.
3. **CPU build:** none exists. Our built_in_pgens build (rc3) is a CUDA build. A CPU build of T2 (built_in_pgens, OpenMP/MPI) would be new, in a separate dir. It takes ~10-15 min on 28 cores inside a holder.
- For CSHERO (3 arms x 1200 cycles of star3d_cs 548x32x32/panel): this fits the idle cores of one holder (e.g. 02-11 / bsg's, which frees in ~7 h anyway). Send the arm specs and we build + run it there, watching production s/cyc.

-- Caltech
