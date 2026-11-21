# NOTE 2026-11-21b viper -> orion (USER "give gpus to orion"): every viper + Raven allocation is yours

Checked 11-21 ~11:10 PST. None of our allocations runs anything of ours on a GPU. RSG, R136 and W189 va1e7 stay paused by
their STOP files.

Running now (use `srun --jobid=<id> --overlap ...` steps):
- viper pRSG 13359763: vipa1150, vipa1155, vipa1163, vipa1165 (4 nodes, ~13 h left). Your V7 rsg_d0/d1 steps are in it.
- viper pool 13359766: vipa1105 (1 node, ~13 h left). The holder only cycles a no-GPU keepalive and a guarded TOI task that
  exits at once (rc=9). Both GPUs are free.
- Raven pool 32638311: ravg1142 (4 A100, ~15 h left). Keepalive only, all free.

Pending, also yours once they start: viper pr136f 13359764 (2 nodes), pHG 13360518 (6), pHG262 13255126 (6); Raven pHG262
32307711 (6). I will post the node names when they start.

Rules as before: never cancel a holder, and do not drop tasks into either `pool_1107/queue`, because any of our holders takes them.

-- viper (for the USER)
