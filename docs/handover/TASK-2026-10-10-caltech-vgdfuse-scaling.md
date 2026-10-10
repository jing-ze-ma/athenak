# TASK viper -> Caltech: H200 timing of vgdfuse-1010 + GPU scaling baseline (short jobs, user GO via main 10-10)

Supersedes the follow-ups of TASK-2026-10-10-caltech-vgdspeed. That task is DONE (its NOTE and results/4304128 are on
this branch). vgdfuse-1010 contains vgdspeed-1009, so its BASE here is the newer rt-integration 5b304cf9.

## What changed (fork branch `vgdfuse-1010`, tip 712ace27; only src/rad_m1/{rad_m1.hpp, rad_m1.cpp, rad_m1_vetgd.cpp, rad_m1_vetcol.cpp})

- **`vet_gd_halo_exact` (default 2 where it applies: MPI halo, compact halo, lists).** The per-(pass, shell) halo
  lists hold only the band values the sweep actually READS. Before, they held the whole mask: on A100 the reads are
  6-9 % of the mask. The lists are built from the sweep's own ray geometry once per direction set, with one request
  message per partner. There are no masks, no dense first exchange and no zero kernel.
  - Bitwise to 5b304cf9 on CPU: AG Car A at 4 and 2 ranks (incl. twin unfused / wall_interp off), BSG hr and He giant.
  - Bitwise on Raven 4 A100: AG Car A hr and BSG hr reduced, run-to-run included.
- **`vet_gd_overlap`, `vet_gd_fuse_shells`**: opt-in, default off.
  - `vet_gd_fuse_shells` = H, a lagged exchange every H shells. It is NOT bitwise and FAILED the accuracy gate on
    A100: interior up to 1e-5, thin atmosphere up to 6e-2, L_top up to 1e-3, against a round-off spread of about 1e-7.
    It was also not faster. It is timed here only because main asked.
- No new input keys are needed for the defaults.

A100 (Raven, 4 GPUs = 1 block per GPU, interleaved, s/cycle over cycles 5-30):

| case | base | new |
| --- | --- | --- |
| AG Car A hr | 0.399 / 0.405 | 0.362 / 0.361 (-10 %) |
| BSG hr reduced | 0.313 / 0.296 | 0.245 / 0.246 (-19 %) |

vet_gd halo kernel time: AG Car -51 %, BSG -68 %.

## Steps

1. **Build** two he_star_m1 H200 binaries exactly as the *_nofma ones (build_inc_nofma.sh: CUDA 12.9, HOPPER90, MPI,
   host -ffp-contract=off), incrementally from one build dir. Report both md5 values.
   - BASE = fork/rt-integration `5b304cf9`, giving athena_gpu_he_star_m1_5b304cf9_nofma.
   - NEW = fork/vgdfuse-1010 `712ace27`, giving athena_gpu_he_star_m1_712ace27_nofma.
2. **Inputs.** Use the files dirs set up by agcar-prod-1009/SETUP.sh and bsg-hrdet-1009/SETUP.sh. These are the ones
   the prof/vgdspeed tasks used, and they contain `agcar_rcxA_ge_accel_st_mg_hr.athinput` (AG Car **A** hr) and
   `bsg_hr_dc5.athinput`. Copy `docs/handover/caltech-vgdfuse-1010/` to `/resnick/groups/carnegie_poc/jingze/vgdfuse_1010/`.
   Check the CONFIG block of both sbatch files (binary names, input paths, gpuwrap, kernel-timer .so).
3. **Smoke.** Run `bash vgdfuse_run.sh <dir> smoke <NEW> <AG Car A input> 2 3 ""` inside a 2-GPU allocation of ≤ 10 min.
   It must cycle with rc 0 and no FATAL.
4. **Job 1**, `sbatch j1_1node.sbatch` (1 node x 4 H200, ≤ 1 h):
   - (1) timing: BASE / NEW / NEW+H2 / NEW+H4, interleaved, 2 repeats, AG Car A hr and BSG hr reduced at 1 block per
     GPU, plus one kernel-timer run of BASE and NEW per case;
   - (2) AG Car A strong scaling on 2 and 1 H200;
   - (3) BSG hr (16 blocks) on 4 and on 2 H200, the latter also with `vet_gd_twin_lowmem` (OOM allowed and reported);
   - (4) weak scaling: BSG reduced on 1 H200 (4 blocks/GPU) vs BSG hr on 4 H200 (4 blocks/GPU).
5. **Job 2**, `sbatch j2_2node.sbatch` (2 nodes x 4 H200, ≤ 45 min): BSG hr on 8 H200, BASE + NEW + NEW repeat.
   It may run in parallel with job 1. If 2 nodes do not backfill within ~2 h, report that and skip it.
6. **NOTE** back as `docs/handover/NOTE-2026-10-10-caltech-vgdfuse-scaling.md` with RAW tables, plus results/<jobid>/
   containing results.txt, the j*.out files and the kp_reader tables:
   - every ARM line (s/cycle mean / median / steady without the rotation cycles, peak MiB, nranks);
   - bitwise: diff the `*.hst` data rows (lines not starting with `#`):
     - NEW vs BASE in every arm pair;
     - NEW rep 2 vs rep 1 (expected identical; H2/H4 are expected to differ from BASE);
   - the strong-scaling table:
     - AG Car A on 1/2/4 H200;
     - BSG hr on 2/4/8;
     - parallel efficiency E(N) = T(N0)·N0 / (T(N)·N), using the steady column, for BASE and NEW;
   - the weak-scaling ratio T(BSG hr, 4 GPU) / T(BSG reduced, 1 GPU), for BASE and NEW;
   - kernel groups of the kt runs (prof_group.py, with `m1_vgd_hx_` and `m1_vgd_hl_` in the vet_gd halo group).

## Paste line for main

`Read docs/handover/TASK-2026-10-10-caltech-vgdfuse-scaling.md on fork/bsg-files-1009 and run it (2 short H200 jobs: 1 node x 4 ≤ 1 h, 2 nodes x 4 ≤ 45 min).`
