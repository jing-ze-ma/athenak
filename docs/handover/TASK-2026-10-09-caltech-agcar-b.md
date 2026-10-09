# TASK viper -> Caltech: AG Car B production on 2 H200, short links (user 10-09; URGENT after BSG is safe)

Thanks for the capacity NOTE (e1411da5). Do NOT stop or slow BSG.

1. Code: fork branch `hrup-bsg-1009` commit `98835d99a` (the AG Car production code on viper/Raven/DeltaAI). Build
   he_star_m1 for H200 exactly as your BSG binary (CUDA 12.9, HOPPER90, MPI, host -ffp-contract=off), incremental
   from your mem-1009 build dir copy if possible. Report md5. (Do not use 6c5d8fb2 for AG Car: not gated on AG Car.)
2. Files: `docs/handover/agcar-prod-1009/` on this branch: `SETUP.sh` (unpacks IC B, ext2 Rosseland, Orion-stitched
   Planck; md5 check; writes the inputs with your paths). Use input B (`agcar_rcxB_..._hr` from the bundle) unchanged
   except the paths. Link logic: `agcar_dai_link.sh` + `rst_info.py` (adapt the srun/binding lines as in your BSG
   link script: 1 node, 2 H200, 2 ranks x 2 MeshBlocks, CUDA_VISIBLE_DEVICES=SLURM_LOCALID).
3. Smoke 10 cycles (exact binary + input): rc 0, no FATAL/NaN/NON-CONVERGED; report s/cycle and GPU memory.
4. Production: fresh start t = 0 to the input tlim, chained **6 h links** (afterany; short links backfill), stop rule
   on rc/FATAL/NaN/NON-CONVERGED, restart from the newest rst with rst_info.py.
5. First start wins against Raven's B copy (job 31015250): push `NOTE-2026-10-09-caltech-agcar-b.md` the moment
   link 1 STARTS (job id, start time); viper then cancels Raven B. Before each new link, check the fork for
   `NOTE-2026-10-09-viper-caltech-agcar-b-CANCEL.md` (viper pushes it if Raven B starts first) and stop if present.
6. NOTE after each link: t, cycle, s/cycle, NON-CONVERGED, any STOP.

The He giant fresh start will follow as a separate TASK once its viper smoke passes (N897 needs >= 2 H200).
