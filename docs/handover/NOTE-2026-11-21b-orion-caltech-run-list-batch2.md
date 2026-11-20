# NOTE 2026-11-21b orion -> caltech (USER): run list batch 2 (10 arms + 1 optional), after batch 1

Follows 111be2f40 (NOTE-2026-11-21-orion-caltech-accept-h200-holders). Everything is set up as in batch 1: `armc.sh`, 4 ranks on
1 node per arm, rsg.01295 / hegiant.00139, results to `docs/handover/caltech_results/<TAG>/`. Please start these on nodes as
batch 1 frees them. The arms are in `docs/handover/caltech_runs_1121/arms_batch2.txt`; new keys files are in the same directory
(`rsg_thv7`, `hg_thv7`, `rsg_v19`, `hg_v19`, `rsg_fidp1`, `hg_fidp1`; `@J@` as before).

New builds (same cmake line as batch 1):
| build name | branch | commit |
|---|---|---|
| B_thinf | caltech-runs/thinf | 2aaaaa5d19b7c6673ec19831ecce1ef837cfa693 |
| B_geom-ca | caltech-runs/geom-ca | 25ed67139b5be083439633a27b8ef0f3dcdd04ea |
| B_geom-fm2 | caltech-runs/geom-fm2 | ae8ead4512713b5f586b8dc332fdcb523c51ddd5 |
| B_speedhead | caltech-runs/speedhead | fe6e94f03f8e5992c5e749b5630c84478388f9b9 |
| B_realf | caltech-runs/realf | 4d1904a7dc8c270276d94ff8890ec9debcd5e7d6 |

Priority: P6 thinf (THF300_rsg, THF100_hg, THF_hg_rt_v7) > P7 geometry (GCA300_rsg, FMP300_hg) > P8 v19 head (HGV40_12, URF_rsg /
UR8_rsg) > P9 real F (RF300_rsg_f1, RF100_hg_f1) > P10 optional nsys profile.

Several of these arms are single runs whose reference was a Raven run from a newer restart. **The reference arms on the same
restart (rsg.01295 / hegiant.00139) will come in a batch-3 NOTE.** Please run batch 2 as listed; the references can run later
on the same nodes.

Expected times: 300-cycle RSG about 1.5-2 h; 300-cycle HG about 3 h; 100-cycle HG about 1 h; the rest under 30 min.

-- orion (for the USER)
