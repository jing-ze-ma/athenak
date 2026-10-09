# NOTE Delta -> viper: prof batch2 closed without a Delta result (answered by Caltech); 1.18 GPU-h total

- 22778535 (default SRUN4 with --gpus-per-task=1 --gpu-bind=closest): hung in agc_plain, cancelled at 8:52,
  **1.18 GPU-h** - see NOTE-2026-10-09-delta-gpu-bind-hang.md (fb16dab9).
- 22779484 (resubmit, SRUN4="srun -n 4 --cpu-bind=cores"): cancelled while PENDING at 16:25 CDT (Slurm est. start
  ~21:14; fairshare 438 -> 405 after the hung run) when NOTE-2026-10-09-deltaai-prof.md appeared on bsg-files-1009.
  That NOTE says DeltaAI also cancelled because Caltech answered first (NOTE-2026-10-09-caltech-prof, 7a20aeec,
  H200 x2 job 4297456). 0 GPU-h.
- No Delta profiling numbers. The plain-srun binding is therefore still untested for prof_cuda.sh itself (it is the
  one spb2 prepA ran with). Setup kept on Delta (prof/agcar_files, prof/kt kernel timer, bin 98835d99 md5 b3f01de1).
- Batch total: 1.18 GPU-h (cap 4). Delta is idle; no batch pending.
