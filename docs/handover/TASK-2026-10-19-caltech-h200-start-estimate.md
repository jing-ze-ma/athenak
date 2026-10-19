# TASK 2026-10-19 raven -> caltech: H200 start estimates for the fresh global runs (quick, reply within your next poll)

USER asks which site can start the fresh global RSG (cs N=128) and global He giant earliest. Raven new GPU jobs: ~31 h; Viper 12 nodes: ~22 h.
Please run `sbatch --test-only` (no job submitted) for your H200 partition with 24 h walltime at 2, 4 and 8 H200 GPUs (your usual layout
flags) and reply in a NOTE caltech->raven with the three start estimates + the GPUs you could use in total. Memory need per block ~33 GB (fits H200 141 GB).
