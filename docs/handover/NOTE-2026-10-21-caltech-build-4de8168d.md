# NOTE 2026-10-21 caltech -> raven: 4de8168d built (TASK-2026-10-21-caltech-queue-build step 2)

Branch t2befb-ncfix-1cb @ 4de8168d2b23acf08d2c83fec22fe9d53b04fd6e, same recipe as the 3896d529 builds
(build_inc_nofma.sh, H200 nofma, login node, incremental: 45 objects recompiled each). Both BUILD_OK, no runs.

| binary | md5 |
|---|---|
| RSG he H200: builds_fin51/bin/athena_gpu_he_star_m1_4de8168d_nofma | 9a56b4a109c575b7aa87a9a94d0124eb |
| accretor nofma H200: builds_accr/bin/athena_gpu_ry_per_accretor_4de8168d_nofma | 0ff5fc545c9694302cfec8a3e6169af4 |

`readelf -p .comment` (both): `GCC: (GNU) 11.5.0 20240719 (Red Hat 11.5.0-2)` + `GCC: (Spack GCC) 13.2.0` (same as 3896d529).
Production untouched: RSG E 5126654 chain + rdiff 5140304 still PENDING; accretor stays cancelled until the v12lzc NOTE.
