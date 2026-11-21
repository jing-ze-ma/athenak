# NOTE 2026-11-21 viper -> orion (USER): Raven pool 32638311 (ravg1142, 4 A100) is free for you now

The Raven pool holder started 11-21 01:18 PST and first picked up a leftover viper task, dhj W189 va1e7. On the USER's word
("give to orion") that run is paused: step killed, rst kept, `wasp189b_mhd100_va1e7_1116/STOP`. The holder (holder.sbatch)
stays alive until wall - 10 min and only cycles a no-GPU keepalive task, so all 4 GPUs on ravg1142 are yours.
Use `srun --jobid=32638311 --overlap ...` steps, as on viper. Never cancel the holder. Also do not drop tasks into
`/raven/ptmp/jinma/pool_1107/queue`: any of our holders takes them.

-- viper (for the USER)
