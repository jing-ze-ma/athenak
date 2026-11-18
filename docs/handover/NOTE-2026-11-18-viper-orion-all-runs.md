# NOTE 2026-11-18 viper -> Orion (USER): every run viper has, running or queued (snapshot 11-17 18:53 PST)

The USER asked viper to list all its runs for you. Times are PST. Supersedes the run list in NOTE 46b0c0780.

## viper: running
| job / step | nodes | what | ends |
|---|---|---|---|
| pRSG 13199029 step r136_C2 (and following r136_*) | vipa1087 (2 GPU) | cc test: R136 r136.00095 gate of clean-cc-physfix2-viper (60-cycle base vs final, BE-only arm, restart + forced retry); physfix2_1117 | before ~19:45 |
| pRSG 13199029 step r136_D2 (and following r136_*) | vipa1164 (2 GPU) | cc test: R136 defaults vs clean_apply_store=float + clean_tol_opac keys; ccdef2_1117 | before ~19:45 |
| pRSG 13199029 step ccf_C1 (and following ccf_*) | vipa1162 (2 GPU) | cc test: cubed-sphere gates (clean_consistency on/off) + R136 smoke; ccfull_1117 | before ~19:45 |
| pRSG 13199029 | vipa1108 | nothing of viper's now (W189 100 G filler, cancelled by you, rst dhj.00156) | - |
| pool 13199033 step r136_XC64 | vipa1089 | last speed8 memory step (timeout 240 s, started 18:52); viper starts no new steps on 13199033 | ~18:57 |
| pr136f 13199030 step memmon_wrap_viper.sh | vipa1151-1152 | **R136 production** (old prod, link chain) | holder end ~19:50; successor pr136f 13255124 continues it |
| pHG 13199031 | 6 nodes | nothing of viper's (HG production stopped by you 18:04, STOP set; your athena_rd6 runs there) | - |

## viper: stopped productions / fillers (to re-home later; the USER decides where)
| run | newest rst | stopped by |
|---|---|---|
| HG (He Giant Global) | hegiant.00186 (17:28) | Orion 18:04 (USER), prod_1107/hg/STOP |
| HG262 | hegiant.00046 (16:31) | Orion 16:33 (USER), hg262_1113/run_viper/STOP |
| M35MMS | m35mms.00098 (18:44) | viper 18:46 (USER: stop runs for cc tests) |
| M35ZAMS | m35zams.00026 (18:07) | viper 18:48 (USER) |
| W189 30 G | dhj.00175 (18:33) | viper 18:48 (USER) |
| W189 100 G | dhj.00156 (18:42) | Orion ~18:47 |
| W76 / TOI-1518b / W178 (dhj) | 00400 / 00372 / 00193 | paused for tests 11-17 (USER) |

## viper: queued
- Successors (afterany): pRSG 13255123 (4 nodes), pr136f 13255124 (2), pHG 13255125 (6), pHG262 13255126 (6, waits for the
  requeued 13199032, PENDING Priority), pool 13255127 (1). Split between Orion tests / viper cc tests / fillers: USER's call,
  viper posts it when decided.
- apudev: lzq_n5 13255938 (cc rotating-sp L_z test, realisability work).
- Held by the USER: pw121h 12810139-41.

## Raven
- Running: nothing of viper's. The orion_cpuhold1-4 jobs are yours.
- RSG E: STOPPED since 11-17 10:13 (GPU OOM in pHG262 32307709). The watcher rsg_in_phg_1116.sh (raven02 pid 60939) starts RSG
  from rsg.01658 in pHG262 32307710 when it starts (PENDING Priority), then 32307711 (afterany).
- Queued pool holders (1 node, tests): 32435128 (W189 va1e7 twin task), 32490420, 32490421.
- Held by the USER: pw121h 32037237-38.

## viper agents working now (all cc)
physfix2 (floor-root rule + R136 gate), defaults2 gating, speed8 halo lists (report pending), the cs audit + full-problem test
prep (R136 long, RSG E and HG restarts with clean-cc; nothing launched), realisability L_z (apudev).

-- viper
