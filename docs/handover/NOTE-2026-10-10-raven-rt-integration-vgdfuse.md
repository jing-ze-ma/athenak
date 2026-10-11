# NOTE 2026-10-10 (raven): vgdfuse-1010 merged into rt-integration

viper has been down since about 07:20 CEST on 10-10. The Raven session now carries on with its work. It made this
merge at the user's request ("merge it now", 10-10). It did not wait for the Caltech H200 timing
(TASK-2026-10-10-caltech-vgdfuse-scaling on fork/bsg-files-1009). That task still stands; its BASE 5b304cf9 and NEW
712ace27 are unchanged.

## What merged
- fork/vgdfuse-1010 712ace27 is merged into rt-integration 5b304cf9 as 6ec86a56 (non-ff, no conflicts).
  rt-integration had not moved since the branch point, so `src/` after the merge equals 712ace27's `src/` exactly.
- Files: src/rad_m1/{rad_m1.hpp, rad_m1.cpp, rad_m1_vetgd.cpp, rad_m1_vetcol.cpp}.
- Default change, bitwise: `rad_m1/vet_gd_halo_exact` = 2 where it applies (MPI halo, compact halo, lists). The halo
  lists then hold only the band values the sweep reads.
- Opt-in, default off: `vet_gd_overlap` and `vet_gd_fuse_shells`. `vet_gd_fuse_shells` is NOT bitwise and failed
  viper's accuracy gate (thin atmosphere up to 6e-2). Do not use it in production.

## Gates (Raven gpudev, 1 node x 4 A100, 1 block per GPU; viper's runs, read from the job logs)
- BASE = /raven/ptmp/jinma/vgdspeed_1009/bin/athena_he_a100_89693b59_vgdmerge (md5 a1041576..., same src as 5b304cf9).
- NEW = /raven/ptmp/jinma/vgdfuse_1010/bin/athena_he_a100_712ace27_fuse12 (md5 0f9217f0...).
- Script: /raven/ptmp/jinma/vgdfuse_1010/raven_time.sh. Runs interleaved, dumps off, hst compared with the header
  stripped.

| case | jobs | BASE s/cycle (mean of reps) | NEW s/cycle | change | hst NEW vs BASE | run-to-run |
| --- | --- | --- | --- | --- | --- | --- |
| AG Car A hr, 30 cycles, 3 reps | 31041875 | 0.411 | 0.356 | -13 % | identical (all reps) | identical |
| BSG true-repro hr, reduced to 4 blocks, 30 cycles, 3 reps | 31041876 | 0.299 | 0.252 | -16 % | identical | identical |
| He giant fresh, 12 cycles, 2 reps | 31042231 | 1.479 | 1.241 | -16 % | identical | identical |

- NEW printed `vet_gd_halo_exact lists made` with fallbacks 0 and unsent-source reads 0 in every run.
- Job 31041877 (He giant, first try) ended with rc 1 and 4 FATAL in BOTH arms, BASE included. That was an input/run
  setup problem in the attempt, not the code. The retry, 31042231, passed.
- Kernel-timer runs: job 31041878 (AG Car and BSG, one run per arm).
- Limits: these gates compare hst only, not bin data. The CPU bitwise checks against 5b304cf9 (AG Car A at 4 and 2
  ranks, BSG hr, He giant) are from viper's TASK-2026-10-10-caltech-vgdfuse-scaling. Their logs are on
  /viper/ptmp2, which Raven cannot reach.

## Production
Nothing changes in the running productions. The held AG Car jobs on Raven (agcA_rc / agcB_rc) are untouched.
