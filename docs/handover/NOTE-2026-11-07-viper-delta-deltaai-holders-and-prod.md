# NOTE 2026-11-07 viper -> Delta + DeltaAI (USER: "can you also do those on delta or deltaai? can you ask them if they can queue allocations like what we did on caltech")

**Answer with a NOTE on this branch, facts only. Start nothing and queue nothing until viper sends a TASK.** Delta: thanks for 82eca4c38 (queue, 935 GPU-h, CUDA 12.9 build); this extends it.

## What "those" are
The Caltech productions now continue off Caltech (HANDOVER-2026-11-07-caltech-to-raven.md), queued on Raven AND viper, first machine to start wins (claim files): RSG E, HG, R136a1, plus BSG test arms and the default-candidate tests (NOTE-2026-11-07-viper-delta-deltaai-can-you-run-defaults-tests.md). USER asks whether Delta / DeltaAI can be a third/fourth place for them.

## Questions (both)
1. **Holder allocations, Caltech-style:** can you queue long allocations (e.g. 24 h or your max wall) that sit in the queue now, and when they start run whatever task viper puts in a queue directory (inside the allocation, srun steps, one after another; exit after ~60 min idle)? Allowed by your site policy / QoS? Max wall, max nodes per job, max running and pending jobs per user, and does an idle allocation bill fully?
2. Queue estimate now (`sbatch --test-only`) for: 1 node / 24 h; 2 nodes / 24 h; 3 nodes / 24 h; 6 nodes / 24 h.
3. Budget left (GPU-h or SU) and the billing factor per node-hour.
4. prod-1107v 698b82413 build status (he_star_m1, CUDA, MPI, `-ffp-contract=off` + device `--fmad=false`); compile-only check is enough for now.
5. Free fast scratch (RSG rst 45.4 GB, HG 97.5 GB, R136a1 12.1 GB).

## Sizing per run (Caltech per-block memory incl. pools: RSG ~19 GB, HG ~25 GB; R136a1 est. 18-30 GB)
| run | blocks | Delta (4x A100 40 GB/node) | DeltaAI (4x GH200 96 GB/node) | Caltech speed |
|---|---|---|---|---|
| RSG E | 24 | 6 nodes, 24 x 1 | 2 nodes, 8 x 3 (Caltech layout, ~57 GB/GPU) | 1.35 s/cycle on 8 H200 |
| HG | 24 | 6 nodes, 24 x 1 | 3 nodes, 12 x 2 (~49 GB/GPU) or 2 nodes 8 x 3 (~75 GB, untested) | 0.97 s/cycle on 12 H200 |
| R136a1 | 4 | 1 node, 4 x 1 | 1 node, 4 x 1 (or 2 x 2) | 1.15 s/cycle on 2 H200 |
| BSG test arm | 16 | 1 node, 4 x 4 | 1 node, 4 x 4 | ~15 h / 57.4 d on 4 H200 |
DeltaAI is the natural fit for RSG/HG (GH200 memory = Caltech layouts, 2-3 nodes instead of 6). Delta at 6 nodes x 24 h = 576 GPU-h per day per run: only 1-node holders (R136a1, BSG arms, tests) look affordable there.

## Data route
Restart files only via the USER (Caltech `$J` or Raven `/raven/ptmp/jinma/caltech_1107/`, `raven_gap_1105/{rsgE,hg}/rst/`). The parallel-stream copy (pcp.sh, NOTE sec 7 of the handover-received note) gets ~5-10x a single rsync. Say which source and how fast inbound copies are.

## Rules you would run under (same as Raven/viper)
Preflight every command line before a job can start (`preflight.py` logic: fresh starts FATAL on command-line keys absent from the deck; restarts need the block in the rst input); stop on FATAL / summary NON-CONVERGED=[1-9] / NaN / rc != 0; kill exact PIDs/jobs only; keys and binaries only as viper's TASK says; claim per run so only one machine continues a production (viper mirrors claims; you would report your claim by NOTE before srun and wait for viper's ack — slower, so production on Delta/DeltaAI only if viper hands the run over explicitly).
