# TASK 2026-10-24 caltech -> viper (USER): push the run12 deck to git, no rsync needed

Thanks for NOTE d880eaf7 / fd3a8e74. Caltech already has the 3 `in/` tables (acc3d_prod_1018/files: md5 5b813886, d81b7a46, 99d668f6,
identical to your bundle) and is building 958cac78 ry_per_accretor (H200 nofma). The ONLY missing file is the deck.

Please commit `acc3d_prod_D_B1v6_extrap2_caltech_abs.athinput` (md5 a93b8f05; paths -> /resnick/groups/carnegie_poc/jingze/acc3d_run12_move/in/)
plus `scripts/link_viper.sbatch` (for the phase/stop logic) to bsg-files-1009 under docs/handover/acc3d-run12-move/ and push.
Caltech then runs the static port check (nlim 10) on 1 node x 4 H200, starts production, and posts the STARTED NOTE (start rule as in your NOTE).
