# TASK 2026-10-19 raven -> caltech: H200 determinism + restart check of the RSG g25 build (before trusting link 1)

USER rule (10-19): every new build/configuration gets (a) a run-to-run determinism check and (b) a restart-equivalence check, both
BITWISE, on the target GPU and layout. The RSG g25 smoke on Raven (A100) and Viper (MI300A) passed both; the H200 build
`athena_gpu_he_star_m1_a2713be0_nofma` (md5 e31ea39e) has not been checked.

Please run on the H200 (any short queue, 1 node is fine; your production decomposition of 3 MeshBlocks/GPU if possible, otherwise the
reduced smoke layout from `docs/handover/rsg-g25-prod/` / Raven `rsg_g25_1019/jobs/smoke.sbatch`: cs N=48, 200 cycles):
1. C  = 200 cycles from t=0 with the production deck (rst at cycle 100, hst every cycle).
2. C2 = an identical rerun of C.                         -> C vs C2: hst + final bins bitwise (data words; bin headers carry the parameter dump).
3. R  = restart from C's cycle-100 rst to cycle 200.     -> R vs C: overlapping hst rows + final bins bitwise.
Report: IDENT/DIFF per comparison (first differing cycle if DIFF), s/cycle, memory per GPU. If DIFF: HOLD links 4993305-07 and tell us.
Do NOT cancel or modify the queued production links otherwise; this check can run while they wait (est. start ~16:47 CEST).

FYI (not Caltech-specific): we found that on Viper MI300A the <rad_m1> implicit_halo_mpi exchange gives wrong/random ghost cells with
>= 2 ranks; A100 (CUDA, same code) is correct and deterministic, so H200 is expected to be fine — this check confirms it.
