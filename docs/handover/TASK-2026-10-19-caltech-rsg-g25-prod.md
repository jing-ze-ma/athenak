# TASK 2026-10-19 raven -> caltech: QUEUE NOW the RSG g25 production on H200 (race with Raven)

USER (10-19): RSG g25 (the star Orion runs with Arepo2 + IDORT) on our method. The user wants the production queued NOW so the
job builds priority. Files follow later. Race rule: the first site to START wins. The other site is cancelled only after the
winner has started. Raven chain: jobs 31285583 -> 31285584 -> 31285585 (6 nodes x 4 A100, 24 h links). Their start estimate is
~10-20 04:45 CEST.

## The run
- Star: MESA r15140 golden25 profile438, Orion's model (see `NOTE-2026-10-19-orion-rsg-g25-setup.md`).
  - M 20.54 Msun, L 1.79e5 Lsun, R_ph 843.8 Rsun, Teff 4088 K, X 0.6538, Z 0.0198.
- Pgen he_star_m1. Same production method as the He giant 14p5M production deck (`hg145_prod.athinput`), plus:
  - mlt_closure adaptive;
  - **cubed sphere, full sphere, N = 128 per face**;
  - nx1 591 (r 0.03-3 R_ph);
  - **24 MeshBlocks 591 x 64 x 64, about 58M cells**.
- Size: ~7.6 KB/cell measured on GH200, so ~440 GB in total. Suggested layouts:
  - **2 nodes x 4 H200 = 8 GPUs, 3 blocks per GPU**: 7.3M cells, ~55 GB per GPU (recommended);
  - 6 H200 with 4 blocks each: ~74 GB per GPU, if you prefer fewer GPUs.
  - Choose and report it. Number of ranks = number of GPUs; the blocks are spread evenly.
- Cost estimate: dt ~230 s; ~0.1 h wall per simulated day on 24 A100. Expect about the same on 8 H200.
- tlim 1e9 (long-tlim rule). Use your longest walltime links, chained afterany (3 links), the usual STOP rule:
  - stop only on the end-of-run summary NON-CONVERGED=[1-9], DIVERGED, FATAL, NaN or rc != 0;
  - never on per-pass NC lines;
  - prune rst, keep 3.

## Do NOW
1. **Build** the current tip of `origin/fin51-r-fix`, which is a2713be0 now (cs second order + AM conservation, gated).
   - Use your production GPU recipe (H200, -ffp-contract=off as the fin51 builds), PROBLEM he_star_m1.
   - **If fin51-r-fix moves before READY, rebuild from the new tip.** We will say which tip the deck was gated on.
2. **Queue** the chain (link 1 + 2 afterany links) with a link script that reads a `run.cfg` AT JOB START.
   - If `READY=1` is missing (READY=0 or no run.cfg): exit 0 at once. Raven does the same and does NOT resubmit, because a
     resubmitted job loses its accrued priority.
   - If READY=1: run as your usual production link (FRESH from the deck if no rst, else the newest valid rst).
   - Reply with the job ids and the start estimate in a NOTE.
3. **Stage**: the files will be put in `docs/handover/rsg-g25-prod/` on bsg-files-1009 when the Raven smoke passes:
   - the deck (Raven layout; you change only the meshblock/rank layout);
   - IC column `ic_rsg_g25.txt`;
   - the Rosseland (MESA-kap) and Planck tables;
   - the run.cfg values.
   Fix the paths in the deck to your copies. Set READY=1 only when we post "READY" in a NOTE.

## Report (NOTE on bsg-files-1009)
- build commit and md5;
- job ids, layout (GPUs, blocks per GPU) and start estimate;
- later: start time, the first link's s/cycle, memory per GPU, Picard summary.
