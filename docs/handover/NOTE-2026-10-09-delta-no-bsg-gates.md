# NOTE Delta -> viper: please stop sending BSG gates to Delta (user 10-09)

**Request (user-approved): route BSG gates / "first answer wins" BSG tasks to DeltaAI (and viper/Raven), not to Delta.
Keep Delta for 1-node A100 test arms (spb2 / beam / csfloor style).**

Why: on 10-09, DeltaAI answered all three BSG tasks before Delta could start, so Delta built, set up and queued for nothing:

| Delta batch (rad-beam-1008) | Delta | DeltaAI |
|---|---|---|
| beam batch3 (hrdet gate, 4 nodes x 4 A100) | job 22772932 PENDING, Slurm est. start ~1.5 d out; cancelled (NOTE edd27b30) | 2 x GH200 nodes, 3.6 min, bit+det PASS |
| mem batch1 (1 node) | job 22774481 PENDING ~30 min, est. start ~1 h out; cancelled (NOTE a0bd7c31) | answered ~30 min after the TASK |
| mem batch2 (1 node) | not submitted: DeltaAI answered while the binary built (NOTE 40e3c2d9) | answered ~20 min after the TASK |

- Delta queue at our fairshare (~0.05): gpuA100x4-interactive starts 1-node jobs in ~1 h at best; 4-node jobs ~1.5 d;
  regular gpuA100x4 tests to start in ~17 days for 4 nodes. DeltaAI's ghx4-interactive starts in minutes.
- A100 40 GB is also the wrong card for BSG memory questions (mem_c4 is 72.7 GB/GPU; 8 blocks/GPU impossible).
- Delta cost of the three: 0 GPU-h (all cancelled while PENDING). Account ~98 of 100 GPU-h.
- If a Delta A100 data point is ever needed, everything is kept: bin/athena_he_gpu_{98835d99,d385de40,af5147cc,6c5d8fb2},
  bsg_hrdet/ (files + gate.sbatch), bsg_mem/ (mem.sbatch, mem2.sbatch, Kokkos Tools memory-events). Ask with a TASK
  that is NOT also sent to DeltaAI.

Delta stays on the runner task (TASK-2026-10-08-delta-test-runner.md) and keeps polling sp-blend2-1008,
rad-beam-1008 and cs-floor-diag-1008.
