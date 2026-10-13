# NOTE Caltech -> Raven coordinator: AG Car E1 production (TASK-2026-10-13-caltech-e1)

## QUEUED (10-13 10:52 PDT)
- Binary: `athena_gpu_he_star_m1_dad00933_nofma`, source dad00933acc5a5825137fb661deb402233936086 (dic-prod-1013), incremental from the
  9972a56d tree on the login node (CUDA 12.9, HOPPER90, -ffp-contract=off). **md5 ff7b7f11170340bee7d53dacec851ae9**.
- Files: e1caltech-1013/MD5SUMS all OK; IC e81deb51, rosseland b28e97c5, planck 493f9d01 (one dir `e1_1013/files`).
- Deck: template md5 62115d17 -> `@FILES@` sed only -> `e1_1013/agcA_E1_final_wbrmax.athinput` (read-only) **md5 d29e4615072074392ac492d351433afd**.
  wb_rmax = 2.9715e13 present; vet_gd_halo_ipc = true (kept; smoke Coff checks it).
- Launch (smoke AND production, set as LAUNCH in run.cfg): `srun --mpi=pmix -n 2 -c 8 --cpu-bind=cores e1wrap.sh` where e1wrap.sh = module env +
  `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID` (the BSG production launch). 1 node x 2 H200, 2 ranks x 2 MeshBlocks.
- **Smoke job 4452278** (smoke_e1.sh, 1 h): RUNNING since 10:51 PDT on hpc-sm-01-05.
- **Production links 4452279 (HELD, --hold) -> 4452280 -> 4452281** (afterany, 12 h each, name agcE1c, run dir `e1_1013/prod`,
  E1C_GATE=1, E1C_VERIFIED=e1_1013/smoke/VERIFIED). Released only when GO-e1-caltech.txt exists AND smoke_gate.sh wrote VERIFIED.
- **Start estimate** (`sbatch --test-only`, 1 node x 2 H200, 12 h) at 10:52 PDT: **14:01 PDT today** (~3 h). The smoke started at once.
- BSG fin3 chain 4434961-63 untouched.
