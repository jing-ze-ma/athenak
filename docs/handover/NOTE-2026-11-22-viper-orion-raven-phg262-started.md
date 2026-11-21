# NOTE 2026-11-22 viper -> orion: Raven pHG262 32307711 started (6 nodes, yours); successors queued

Raven pHG262 32307711 started 11-21 15:19 PST on ravg[1027-1028,1031,1033,1098,1100] (24 A100, 24 h wall). Its batch
(`hg262_1113/jobs/link.sbatch`) is in HOLD mode because `hg262_1113/STOP` is present, so it runs nothing of ours. You are
already running the a100hr steps in it. Do not write `hg262_1113/FIXED` or `CANCEL`.

Successor holders, all yours once they start (USER 11-21):
- viper pRSG 13514310 (4 nodes, plain holder.sbatch), afterany 13359763 (~23:30 PST)
- viper pool 13514311 (1 node), afterany 13359766
- Raven pool 32747999 (1 node), afterany 32638311 (~01:15 PST)

Still pending: viper pr136f 13359764, pHG 13360518, pHG262 13255126.

-- viper (for the USER)
