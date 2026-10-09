# NOTE 2026-10-09 Caltech: BSG true reproduction (answers TASK-2026-10-09-caltech-bsg-truerepro)

## PRODUCTION STARTED 2026-10-09T09:57:39 (PDT) on hpc-sm-02-11 — viper/Raven copies can be cancelled (first start wins)
- Chain 4286643 -> 4286644 -> 4286645 (12 h links, afterany), 1 node x 2 H200, 2 ranks x 8 MeshBlocks, CUDA_VISIBLE_DEVICES=SLURM_LOCALID
  wrapper, --cpu-bind=cores; fresh start (t = 0) from bsg3d_truerepro2_hr_lm.athinput (SETUP.sh, md5 2e213845 with Caltech paths);
  no command-line physics keys. Run dir /resnick/groups/carnegie_poc/jingze/bsg_1009/prod, link script bsg_1009/link_bsg.sh
  (viper rules: rst_info last_time keys on restarts, STOP on rc/FATAL/NaN/NON-CONVERGED, DONE at rst t >= 4.96e6).
- Binary: mem-1009 6c5d8fb2, PROBLEM=he_star_m1, CUDA 12.9 / gcc 13.2 / hpcx, Kokkos HOPPER90, host -ffp-contract=off:
  athena_gpu_he_star_m1_6c5d8fb2_nofma, md5 e9873137c319b5f29d85680b535a450c.
- Smoke 4286551 (nlim 10, same binary/input): rc 0, 0 FATAL/NaN/NON-CONVERGED, Picard mean 3.7 max 5, dt 87.37 s,
  GPU peak 66537 MiB (65.0 GiB, = DeltaAI 3348593), ~1.6 s/cycle after the 5.9 s setup (cpu time 21.7 s / 10 cycles)
  -> ~25 h for 57k cycles.
- He giant N897 stop record (TASK section 1, user 10-09): cancelled 09:38 PDT, jobs 4242927 (running link 2) + 4242928/4242929
  (pending); CANCEL file in n897. Resume point /resnick/groups/carnegie_poc/jingze/hegiant_1007/n897/rst/hegiant.00048.rst,
  t 2030400.0 s (23.50 d), cycle 112735 (live run was 65 cycles further); binary 7adb12d3 + hydro/sp_x2_periodic_image=true.
