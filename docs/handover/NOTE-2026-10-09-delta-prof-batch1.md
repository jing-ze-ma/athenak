# NOTE Delta -> viper: prof batch1 cancelled while PENDING (withdrawn, moved to DeltaAI); 0 GPU-h

Job 22778067 (gpuA100x4-interactive, 1 node x 4 A100, 30 min; submitted 10-09 13:37 CDT, Slurm start estimate 18:52)
was cancelled while PENDING at 14:05 CDT, following NOTE-2026-10-09-viper-delta-prof-withdrawn (80564020). It never ran:
**0 GPU-h**. Setup is kept for reuse: prof/agcar_files (SETUP_OK, 4 x md5 OK), simple-kernel-timer built from the
existing kokkos-tools clone (1e98316, prof/kt -> bsg_mem/kt), bin/athena_he_gpu_98835d99 (md5 b3f01de1).
Delta queue note: priority 1616-1620 (fairshare 616 + interactive partition 1000); interactive start estimates swung
between 1.3 h and 5 h over the 28 min it was pending. Delta is idle; no batch pending.
