# TASK 2026-10-20 raven -> deltaai: capacity check for long global runs (quick, no jobs submitted)

USER asks whether DeltaAI can take a long production run. Candidates: the fresh global He giant (cs, 24 MeshBlocks 698x64x64, ~33 GB
each -> 3 per GH200 = 2 ghx4 nodes / 8 GPUs) or the fresh global RSG (cs N=128, 24 blocks, smaller per block). Viper is ~44 h out, Raven ~20 h.
Please reply in a NOTE deltaai->raven within your next poll:
1. `sbatch --test-only` start estimates on ghx4 for 2 nodes (8 GH200) and 1 node (4 GH200), 24 h and 48 h walltime (your usual account/flags).
2. Remaining SU balance on the allocation, and the SU cost per node-hour on ghx4.
3. Max walltime on ghx4, and whether the nvme/scratch quota fits a long run (~100-300 GB of rst+dumps; prune to 3 rst as the chain scripts do).
Nothing to build or submit yet.
