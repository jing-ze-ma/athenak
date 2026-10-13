# NOTE DeltaAI: FIN3 final validation, 4 long Cartesian IDORT runs (TASK-2026-10-13-deltaai-fin3-cartlong) -- PICKED UP, but ghx4 batch is ~5 days out

**Picked up** 20:35 CDT 10-12.
- Build: fin3-1012 @ 281dd686, plain GPU target (`build_inc_deltaai.sh none_gpu`: built-in pgens, MPI ON, CUDA HOPPER90). Compiling now; md5 follows.
  NB: our DeltaAI builds (this one and all hes_gpu builds) do NOT pass `-ffp-contract=off`; nvcc/CC defaults. Say if you need it (new build dir, ~15 min).
- Job script ready (/work/nvme/bivj/jma20/cartlong_1013/cartlong.sub, dry run OK): ghx4 BATCH, 1 node, 1 h, one single-rank run per GH200
  (each run an `srun --overlap -n 1` step with CUDA_VISIBLE_DEVICES set inside, so the MPI build gets PMI; outputs as run_gpu.sh).
- **Problem: ghx4 batch start estimate for a 1-node 1-h job is 2026-10-18 08:17 CDT** (`sbatch --test-only`, 20:41 CDT; our 12 h ghx4 BSG job
  3367703 has had no estimate since noon). This cannot win the race vs Caltech CPU / Viper.
- Options (coordinator / user to choose):
  (a) submit to batch anyway (queued, free while pending; cancel when Caltech/Viper report);
  (b) run it inside the next BSG interactive link (3369777, est. ~07:35 CDT 10-13) BEFORE the BSG part: costs that link at most 30 min of GPU time
      (Raven CPU 45-120 min; on GPU expected minutes; runs not done by then are reported with their last cycle line);
  (c) SKIP.
  Default if no answer by the time the binary is ready: (a).
