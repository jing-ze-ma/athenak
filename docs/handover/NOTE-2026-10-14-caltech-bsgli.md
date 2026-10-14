# NOTE Caltech: bsgli QUEUED (TASK-2026-10-14-caltech-bsgli)

**No smoke/gate: user decision 10-14, production submitted directly.** run.cfg BLC_GATE=0 (bsgli_hook.sh line 21 checks
VERIFIED only when BLC_GATE=1). No other change from the TASK.

- Binary /resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_8e018272_nofma, md5 4a92728d016609effe6d4586fb2e55bf,
  git rev-parse HEAD 8e0182723da23d8c864c8b079bd24e75689f2d28 (= origin/fin4-linhao-bc-8e018272); -ffp-contract=off, HOPPER90, MPI.
- Deck bsgli_1014/files/bsgLI_bc_caltech.athinput (read-only), md5 d77f9a102c2fcc4db73f97422cee7d56 (template 532fc7af, @FILES@ sed only).
- linhao_profile_full.txt md5 41edc225; ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805; all MD5SUMS pass.
- Layout 1 node x 2 H200, 2 ranks, LAUNCH = e1_1013/e1launch.sh, --exclude=hpc-sm-01-09.
- Run dir /resnick/groups/carnegie_poc/jingze/bsgli_1014/prod (fresh). TLIM 4.96e6, STOP_NC=1.
- Links (-J bsgLIc, 12 h): 4487415 (no dependency), 4487610, 4487611, 4487612, 4487613 (afterany chain).
- Start estimate: sbatch --test-only 2026-10-13T23:49 (hpc-sm-01-05); squeue estimate for 4487415 2026-10-14T11:40.

Next: STARTED when link 1 runs, then per-link reports as in TASK section 6.
