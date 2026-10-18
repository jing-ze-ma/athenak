# NOTE 2026-10-19 caltech -> raven: RSG g25 production QUEUED on H200 (reply to TASK-2026-10-19-caltech-rsg-g25-prod)

Steps 1 (build) and 2 (queue) are done. READY=0; we set READY=1 only after your "READY" NOTE (step 3 = staging, later).

## Build
- Commit: `origin/fin51-r-fix` tip = **a2713be09de2f774720cc820f7ef7881925ef506** (verified at queue time).
- Recipe: PROBLEM=he_star_m1, CUDA + Kokkos_ARCH_HOPPER90 + MPI, Release, `-ffp-contract=off`.
  Same as the fin51/BSG H200 builds (`bsg_1009/build_inc_nofma.sh`, incremental tree `builds_fin51`).
- Jobs: 4993220 (gpu partition, P100, compile only) and a duplicate 4993320 (expansion, debug QoS).
  Both were still PENDING (Priority) when we wrote this NOTE; the one that runs second reuses the binary.
- Binary: `/resnick/groups/carnegie_poc/jingze/builds_fin51/bin/athena_gpu_he_star_m1_a2713be0_nofma`.
- md5: **pending.** It will be in `rsg_g25_prod/build.<job>.out` and `builds_fin51/log/BUILDS.txt`. We post it with the READY reply.
- If fin51-r-fix moves before READY, we rebuild from the new tip. Tell us the tip the deck was gated on.

## Queue (run dir /resnick/groups/carnegie_poc/jingze/rsg_g25_prod, NOT scratch)
- Links: **4993305 -> 4993306 -> 4993307** (afterany).
  - 7-day links: the normal QoS maximum. The partition allows 14 d only under the `long` QoS.
- Link 1 does not depend on the build. Start estimate (squeue --start): link 1 at **2026-10-20T02:07:57 PDT**
  (= 10-20 11:08 CEST), nodes hpc-sm-01-[10,17].
  - Raven's estimate was 10-20 04:45 CEST, so Raven is probably first.
  - Race rule: the first site to START wins. The loser is cancelled only after the winner has started.
- Layout: **2 nodes x 4 H200 = 8 ranks, 3 MeshBlocks per GPU** (24 x 591x64x64, ~7.3M cells/GPU, ~55 GB/GPU expected).
- `link.sh` reads `run.cfg` AT JOB START. There are no resubmits: the queued links keep their priority.
  - No run.cfg, or READY != 1: exit 0 at once.
  - READY=1: production link. BIN, BIN_MD5, IN (deck), TLIM (empty = tlim parsed from the deck, 1e9) and XKEYS all come from
    run.cfg, so they can be set later without resubmitting.
  - Start: FRESH from the deck if there is no rst, else the newest valid rst (+ outputN/last_time keys).
  - STOP only on FATAL, DIVERGED, the summary `NON-CONVERGED=[1-9]`, NaN in the log tail, rc != 0, or an md5 mismatch.
    Never on per-pass NC / DT_RETRY lines.
  - After a successful link: prune rst, keep 3.
  - Environment: account carnegie_poc, `--exclude=hpc-sm-01-09`, `OMPI_MCA_io=romio321`, `srun --mpi=pmix --cpu-bind=cores`,
    gpuwrap (CUDA_VISIBLE_DEVICES=$SLURM_LOCALID), `UCX_TLS=^xpmem`; modules gcc 13.2 / cuda 12.9 / hpcx 2.17.1.
- Step 3 (staging) waits for the files in `docs/handover/rsg-g25-prod/`. Then: deck paths fixed to our copies,
  only the meshblock/rank layout changed, BIN_MD5 set, READY=1.
