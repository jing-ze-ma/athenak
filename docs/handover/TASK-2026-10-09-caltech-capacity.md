# TASK viper -> Caltech: GPU capacity question (user 10-09; answer by NOTE, no runs)

BSG (chain 4286643-45) keeps the 2 H200 of hpc-sm-02-11. The user wants to run more on Caltech: AG Car B
(4 MeshBlocks, ~20 GB GPU per block, 1 node x 2-4 GPUs) and the He giant fresh start (being prepared on viper).
Please answer in NOTE-2026-10-09-caltech-capacity.md on this branch:
1. Which GPU partitions/nodes can the account use besides hpc-sm-02-11 (sinfo; GPUs per node and type; MaxTime;
   your group's share)? Start estimate (`sbatch --test-only`) for 1 node x 2 and x 4 GPUs at 12 h and 24 h.
2. Is there any second H200/H100/A100 node available to us now?
3. Do NOT stop or slow BSG. Do not submit anything; viper will send a task once the answer is in.
