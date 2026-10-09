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

## STOPPED 2026-10-09 ~14:05 PDT (viper NOTE-2026-10-09-viper-caltech-hegiant-stop-for-tests, user 10-09)
- Chain 4288960-65 cancelled: link 1 4288960 after 1 h 46 min (live at t 1.127e5 s = 1.30 d, cycle 5600, dt 20.08 s,
  0 FATAL/NaN/NON-CONVERGED, ~1.11 s/cycle on 2 H200); links 2-6 cancelled while PENDING. CANCEL file in fresh897.
- Kept: run dir /resnick/groups/carnegie_poc/jingze/hegiant_1009/fresh897; newest rst rst/hegiant.00002.rst (t 86404.69 s =
  1.000 d, cycle 4291); bins every 0.25 d.
- Physics so far (for the restart on the fixed scheme): KE_int and its lateral share track our old N445 run within 0.3 %
  (0.69 d: 7.36e41 erg, lateral 0.10); 0.5-d dump rms v_r 1.0 km/s at 5 Rsun, 0.5 at 20, 0.2 at 30-40 Rsun (convection
  growing from the 1e-2 seed; scaffold weight w_mlt = 1).
- The freed 2 H200 are available for the xthinfix test batteries (TASK-2026-10-09-caltech-xthinfix-tests.md when it lands).
  Seeded BSG prod_seed link 1 4294781 stays queued.
