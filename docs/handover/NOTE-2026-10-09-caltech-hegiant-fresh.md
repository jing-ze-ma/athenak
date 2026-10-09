# NOTE 2026-10-09 Caltech: He giant fresh N897 (answers TASK-2026-10-09-caltech-hegiant-fresh)

## PRODUCTION STARTED 2026-10-09T12:29:11 PDT on hpc-sm-02-11 — Raven chain can be cancelled (first start wins)
- Chain 4288960 4288961 4288962 4288963 4288964 4288965 (6 h links; link 1 afterok on the smoke, later links afterany), 1 node x 2 H200, 2 ranks x 2 MeshBlocks,
  CUDA_VISIBLE_DEVICES=SLURM_LOCALID wrapper, --cpu-bind=cores; fresh start t = 0; XKEYS time/cfl_number=0.3 only.
  Run dir /resnick/groups/carnegie_poc/jingze/hegiant_1009/fresh897, link script hegiant_1009/link_hegiant.sh
  (STOP on rc/FATAL/NaN/NON-CONVERGED, DONE at rst t >= 2.592e6; CANCEL note on the fork checked every 5 min by a login watcher).
- Binary athena_gpu_he_star_m1_6c5d8fb2_nofma, md5 e9873137c319b5f29d85680b535a450c (lowmem off: not in the input).
- Input: SETUP.sh hegiant_fresh897.athinput (Caltech paths) + he_ic_eint_from_t = true in <problem> (viper INPUT NOTE),
  md5 ea976e99a5b138ef99c4d120690903c7.
- Smoke 4288956 (nlim 10, same binary/input/XKEYS): passed the stop rule (afterok); GPU peak 70617 MiB; cpu time 1.929635e+01 s for 10
  cycles incl. setup; Picard mean=1.030000e+01 max=1.200000e+01.
