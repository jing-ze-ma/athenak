# NOTE DeltaAI: BSG comparison run STOPPED and cartlong CANCELLED (user decision, 10-12 21:37 CDT)

User: the DeltaAI BSG run is not the production setup (old method d3bf42d7; FIN3 runs on Raven/Viper/Caltech), so it stops here.
- **BSG (TASK-2026-10-12-deltaai-bsgprod, TASK-2026-10-13-deltaai-bsg-d3rst): STOPPED at t = 8.101625e5 s (9.38 d), cycle 9259**, after 2 clean links
  (0 FATAL / NaN / NOT ADMISSIBLE / NON-CONVERGED; convection developing: v_r -61..+19 km/s at 40.6 Rsun, lateral KE share 0.19).
  - Cancelled while PENDING (never ran): link 3 job 3369777 (would have run the d3rst restart check + swap) and the 12 h ghx4 job 3367703.
  - The d3rst swap never happened (no check run); 9ceb9cdc stays built (athena_hes_gpu_9ceb9cdc0ed4, md5 54c0a971) but unused.
  - Chain stopped. Data kept on DeltaAI: /work/nvme/bivj/jma20/bsgprod_1012/run (rst 00000-00011, bins every 0.5 d to 9.38 d, hst, logs).
    It can resume from rst/bsg3d.00011.rst on request (the d3rst binary reads it).
  - Page (https://claude.ai/artifact/MTSTfEYepRNC2kASacVBiY) shows this run to 9.38 d; it is free for the FIN3 production run.
  - GPU use: ~0.92 h + ~0.6 h wall x 4 GH200 (~6 GPU-h) plus the 25-cycle smoke.
- **cartlong (TASK-2026-10-13-deltaai-fin3-cartlong): CANCELLED** while pending (batch job 3370011, est. ~10-18); never ran.
  Binary athena_none_gpu_281dd686cb47 (md5 21427ecb) and job script kept in /work/nvme/bivj/jma20/cartlong_1013 if needed later.
- DeltaAI now has nothing running or queued.
