# NOTE 2026-11-24zzzzc8 Caltech -> orion: t9l62, cells at the velocity ceiling (zzzzc4)

Read-only, CPU. Source: hydro_w bins (float32 primitives after c2p, real cells only); no rst was read (all needed times have bins). Wedge values (theta pi/4-3pi/4, phi 0-pi), not scaled, cgs. Full per-file tables and histograms: caltech_results/T9_VCEIL/*.txt; counters.txt; ceil.py, logs.py.

| file (prod_*/t9l62/bin) | run, binary | t [s], cycle | N >= 0.5 / 0.8 / 0.9 / 0.99 / 0.999 | \|v\|max [cm/s] (/vceil) at r [cm], rho |
|---|---|---|---|---|
| rc4v_1124 00023 (quiet) | T2, 262affd6f | 247248, c2553 | 0 / 0 / 0 / 0 / 0 | 4.65e6 (0.0155), 4.89e12, 4.98e-9 |
| rc3_1124 00043 | rc3, af836f5a3 | 464458, c4967 | 0 / 0 / 0 / 0 / 0 | 6.27e7 (0.209), 2.50e12, 2.49e-9 |
| rc3_1124 00045 | rc3 (restarted 11-24 07:16) | 486035, c5783 | 0 / 0 / 0 / 0 / 0 | 1.28e8 (0.427), 3.36e12, 3.43e-10 |
| rc3_1124 00047 | rc3 | 507605, c6501 | 5 / 1 / 1 / 0 / 0 | 2.87e8 (0.956), 2.38e12, 1.42e-13 |
| rc3_1124 00049 | rc3 | 529235, c7651 | 0 / 0 / 0 / 0 / 0 | 1.25e8 (0.415), 4.44e12, 6.0e-11 |
| rc3_1124 00051 | rc3 | 550845, c8060 | 0 / 0 / 0 / 0 / 0 | 1.25e8 (0.415), 4.54e12, 1.07e-10 |
| rc3_1124 00052 | rc3 | 561631, c8246 | 0 / 0 / 0 / 0 / 0 | 1.18e8 (0.394), 4.76e12, 6.2e-11 |
| fresh_1121b 00148 (**late, old scheme**, j3 a5257b7ee) | old | 1.59846e6 | 0 / 0 / 0 / 0 / 0 | 5.16e7 (0.172), 5.61e13, 1.16e-13 |
| fresh_1121b 00156 (**late, old scheme**, j3) | old | 1.67629e6 | 0 / 0 / 0 / 0 / 0 | 6.34e7 (0.211), 4.29e12, 1.58e-12 |

- **Clips (|v| within 1e-6 of vceil): 0 in every file. Above vceil: 0. Smooth tail 0.99-1: empty in every file.** Float32 primitives carry ~6e-8, so the 1e-6 test is meaningful.
- 00047: the 5 cells >= 0.5 vceil hold mass 7.9e19 g, KE 9.4e35 erg: 3.0e-12 / 4.0e-9 of M / KE (r > 2e12), 3.4e-13 / 4.0e-9 of M / KE (rho > 1e-10). The single cell >= 0.9: r 2.38e12, rho 1.42e-13, T 2.99e5 K, tau_cell 7.9e-4, c_g/|v| 0.027, v_r -8.5e6, v_theta -2.80e8, v_phi 6.1e7. (N = 1, so no distribution.)
- |v| histogram, rho < 1e-11, bins 1e7-3e8: in the .txt files. In 00047: 2 cells (6.3e30 cm3) in 1.906e8-2.391e8, 1 cell (3.1e30 cm3) in 2.391e8-3.0e8, none in 1.519e8-1.906e8; 00049 and 00156 are empty above 3.9e7 and 7.7e7, the quiet file above 1e7. At rho < 1e-14 the separate floor-gas cap |v| <= 3e7 applies (deck floorcap_v / floorcap_rho).
- **What the dumps can't show.** They are snapshots every 10.8 ks. The ceiling acts inside c2p (source, af836f5a3 general_c2p_hyd.hpp: with floor_from_ekin on, E is fixed and the KE scaled, its L_z part kept; with it off, the clipped KE is taken out of E and tallied in vceil_de. vceil_de is 0 in rc3, T2 and S5G, and > 0 in T9HOT2 K, T9ARMS and fresh 1121b), so a cell clipped inside a step and slowed or heated by the dump time leaves no mark, and these files cannot date or count clips.

**Item 4, counters.**
- `skip(dfloor/vceil)` is one counter (DC_SKIP): it **mixes** dfloor and vceil cells and can't be split from the logs. I do not split it. Totals per run and per-cycle maxima are in counters.txt (cumulative since process start; differenced).
- The event counters in t9l62.log are separate: eos_vceil (c2p calls that scaled the momentum) and eos_dfloor, per ~1e4 s interval. Per run, eos_vceil events: rc3 652 in 5 of 73 intervals (c5623-7462); T2 rc4v 0; c4 start 0 (to c104); S5Gon 68, S5Goff 105 (c5279-5578, t 4.75-4.86e5); T9HOT2 B 90, K 37, P1 0; T9ARMS L 11713, S 20478; fresh 1121b 1629 in 33 of 229 intervals. F has no t9l62.log rows.
- eos_dfloor is 0 after the first interval in rc3 and T2, while skip stays at about 1e6 per audit line there ([refresh], quiet phase) and 2e2 (rc3 outburst restart). The two counters don't count the same cells, so I don't form a vceil fraction from them. The eos_dfloor column is not 0 in T9ARMS L (5.5e9) and fresh 1121b (8.8e11).

**Files not used:** rsts (88 GB each in rc3/rc4v, 28 GB in 1121b); the cells of the t9l62 c2553 stop rst are represented by dump 00023 (same cycle). The old-scheme rst 00037-00039 were not opened; their time could not be confirmed from a bin.

-- Caltech
