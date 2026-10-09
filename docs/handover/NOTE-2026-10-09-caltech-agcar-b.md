# NOTE 2026-10-09 Caltech: AG Car B on HOLD acknowledged (answers NOTE-2026-10-09-viper-caltech-agcar-b-HOLD)

- Nothing ran. Pending smokes cancelled before they started: 4287959 (2 H200) and 4288657 (1 H200, added because the
  2-H200 estimate was 19:24 PDT). CANCEL file in agcarb_1009/prod; no production submitted.
- Binary kept: builds/bin/athena_gpu_he_star_m1_98835d99_nofma (hrup-bsg-1009 98835d99, CUDA 12.9, HOPPER90, MPI,
  host -ffp-contract=off), md5 bc6f4e9e11c106a3606246d2ac40d296. Files set up (SETUP.sh OK) in
  /resnick/groups/carnegie_poc/jingze/agcarb_1009/files; link script agcarb_1009/link_agcarb.sh (2 H200) and
  link_agcarb_1g.sh (1 H200, 4 blocks ~80 GB) ready if B or A ever comes back to Caltech.
- BSG continues unchanged (chain 4286643-45; link 1 at ~55 min: cycle ~2700, 0 FATAL/NON-CONVERGED, ~1.17 s/cycle).
