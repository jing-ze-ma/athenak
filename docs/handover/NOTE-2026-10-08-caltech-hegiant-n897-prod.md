# NOTE 2026-10-08 Caltech: He giant N897 PRODUCTION started (10.42 d -> 30 d)

From: Caltech. To: viper, DeltaAI. Docs only.
- User OK 10-08. Chain 4242926 -> 4242927 -> 4242928 -> 4242929 (12 h links, 1 node x 2 H200, 2 ranks, 2 MeshBlocks/GPU),
  binary athena_gpu_he_star_m1_f3a66907 (md5 f5d85053), restart from our remapped N897 hegiant.00021.rst (t 900420.8),
  TLIM 2.592e6. First link with viper's first-link keys (cfl 0.3, restart_refill_ghosts=true, implicit_opac_newton_slope_max=3,
  output last_time keys); later links drop refill_ghosts (automatic, keyed on the newest rst).
- Link 1 STARTED 10-08 12:12 PDT. At 28 min: t 9.251e5, cycle 49200, dt 18.589, 0 FATAL/NaN/NON-CONVERGED, rst 00022 written.
  1.24 s/cycle -> ~31 h to 30 d (3 links + 1 spare).
- Duplicates: DeltaAI link 1 3337674 and viper hegiant897 (12119930-39) are the other copies of this leg. The user decides
  which continue; this note cancels nothing.
