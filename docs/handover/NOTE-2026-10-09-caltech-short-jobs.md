# NOTE 2026-10-09 Caltech -> viper: Caltech can take SHORT jobs now (user 10-09)

User: "let viper know caltech can take short jobs". Snapshot 13:20 PDT:
- Running on Caltech: He giant fresh N897 (2 H200, chain 4288960-65), accretor RHD (1 H200 + queued smokes). Queued: seeded
  BSG true repro (2 H200, link 1 4294781, est. start ~18:25 PDT).
- GPU partition is busy: 1 H200 free, ~200 jobs pending; our fair share is 0.25 (above share). sbatch --test-only for a NEW
  6 h job: 1 x H200 ~10-10 17:10 PDT, 2 x H200 ~10-10 22:10 PDT. In practice short jobs backfill much sooner (today's
  10-30 min smokes started within 0-2 h).
- So: send SHORT tasks (smokes, gates, profiling, <= ~2-6 h, preferably 1 H200, 2 at most). Long 2-GPU productions would
  start ~a day out unless they replace/follow a current run.
- Ready on Caltech without a build: athena_gpu_he_star_m1_6c5d8fb2_nofma (md5 e9873137, mem-1009, host -ffp-contract=off),
  athena_gpu_he_star_m1_98835d99_nofma (md5 bc6f4e9e), athena_gpu_he_star_m1_7adb12d3 (f3a66907 + seam fix); incremental
  builds of other commits take ~5-15 min on a CPU node.
- Put the task in a TASK-2026-10-xx-caltech-*.md on any fork branch as before (Caltech fetches and checks the fork when asked).
