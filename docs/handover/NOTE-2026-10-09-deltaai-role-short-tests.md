# NOTE DeltaAI -> viper: DeltaAI = fast short-test machine (user 10-09)

The user agrees with this proposal: route **short tests** to DeltaAI (gates, smokes, memory/bitwise checks),
keep long **production** on Caltech / viper / Raven. Please write DeltaAI TASKs to fit the limits below.

## Evidence (10-09)
| task | layout | waited | ran | GPU-h |
|---|---|---|---|---|
| BSG half-range gate (3347471) | 2 nodes x 4 GH200, 1 h limit | ~15 min | 4 min | 0.48 |
| BSG mem round 1 (3347781) | 1 node, 30 min | ~30 min | 1.5 min | 0.09 |
| BSG mem round 2 (3348593) | 1 node, 30 min | ~2 min | 2 min | 0.12 (+0.05 failed try) |
| AG Car smokes B/A (3348878/3349002) | 1 node, 20 min | ~33 min / <1 min | 29 s each | ~0.03 |
| AG Car B link 1 (3349057) | 1 node, 2 h (cut to 1 h) | first estimate ~1.75 h; cutting to 1 h moved it up | -- | -- |

DeltaAI answered all three BSG tasks first (Delta 0 GPU-h, NOTE-2026-10-09-delta-no-bsg-gates).

## What fits DeltaAI
- **<= 1 h, <= 2 nodes (8 GH200)** on ghx4-interactive: usually starts within ~5-35 min.
- GH200: 96 GB HBM per GPU (BSG with lowmem: 8 blocks per GPU fit at 65 GiB); 4 GPUs + 288 Grace cores per node.
- Builds: `build_inc_deltaai.sh` (incremental, 4-10 min per commit, worktree per target); md5s recorded.
- Cheap: today's gate + mem tasks were < 1 GPU-h total. Balance ~905 GPU-h (bivj-dtai-gh).

## Limits (please design around them)
- **Interactive: 1 job per user at a time (QOS MaxSubmit 1)** -> arms run one after another. Put several arms in ONE
  job script (as gate_hrdet.sh / mem_cuda*.sh did) rather than separate jobs.
- **Interactive: MaxNodes 2, MaxTime 2 h.** Ask for 1 h or less: 2 h jobs wait 1-2 h.
- **Regular ghx4: ~5 days to start (1 node, any length)** -> no multi-day production here.
- Interactive is charged 2x. Fairshare drops with use (took 3-4 days to recover on 10-01).
- 4-node layouts (16 GPUs as 4x4) are not practical; 2 nodes x 8 ranks (2 per GPU) worked for the BSG gate.
- AthenaK accepts a command-line key only if it is already in the input file: put new keys in the input
  (the mem_cuda2.sh lesson).

## Still running here
AG Car B then A on interactive links (GO 4fd88b70), as agreed; NOTEs on bsg-files-1009.
