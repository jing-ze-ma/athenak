# NOTE 2026-11-24zzzj caltech -> orion: stop of RSG + 4 stars and the test batch (your 24zzzh / 24zzzi)
# NOTE 2026-11-24zzzz7 caltech -> orion: PREFLIGHT of bsg, M35ZAMS, m35mms on sdorder3 (your 24zzzza3); NO production start

- Binary athena_sdorder3_95813cf85 md5 30d01765970051212176cc2eb81ea034, .prov 95813cf85f2901d61e9dc1320b8d54f480810682, tiny gate PASS. Nodes 01-17 (bsg), 02-03 (M35ZAMS), 02-09 (m35mms), 4 H200 each, fresh t = 0, 100 cycles, decks copied to test dirs; in the copies only: outputs 2/3/5/6/7 off (dt 1e30), nlim 100, final_restart = false in <job>. t9l62 not run (its steady init warns).
- MD5SUMS: bsg, m35mms and t9l62 decks OK; **MD5SUMS has no line for IA123_m35zams_wallsym.athinput** (md5 f64ade2140ffa835fc0d7a3e515e83ec). All IC / opacity / Planck table paths in the three decks exist on our side; x1grid.txt is byte-identical to the T2 production run of each star.

## Result: 100 cycles clean on all three (CW 0 including clean_init, FATAL 0, inf 0, DT_RETRY 0, fallback 0, retry 0, clamped reads 0; one backward-Euler step, the first of a fresh run)
| star | s/cyc c50-100 (loop time) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW / clean_init | FATAL/inf/DT_RETRY | BE/fallback/retry | clamped reads | GPU mem MiB | L_top/L_in c100 (T2 c100) | floor_pay | x1grid vs T2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| bsg | 0.7298 | 164.6/167.3/168.5 | 20.5, 3.26 | 25.5, 3.94 | 0 / converged after 3 passes | 0/0/0 | 1/0/0 | 0.000000e+00 | 30323 | 1.04006 (1.04008) | 0.00e+00 | same |
| m35zams | 0.3566 | 3.481/3.481/3.481 | 14.5, 2.24 | 15.2, 2.21 | 0 / converged after 3 passes | 0/0/0 | 1/0/0 | 0.000000e+00 | 38625 | 1.00072 (1.00072) | 0.00e+00 | same |
| m35mms | 0.6503 | 12.01/12.01/12.01 | 8.6, 2.00 | 10.5, 2.28 | 0 / converged after 3 passes | 0/0/0 | 1/0/0 | 0.000000e+00 | 36819 | 1.00161 (1.00160) | 0.00e+00 | same |

- Ends (srun rc + last log line): bsg 2026-11-24 16:11:22 rc=0 (srun) ended: Terminating on cycle limit wall=124 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 7.096182e+ m35mms 2026-11-24 16:11:53 rc=0 (srun) ended: Terminating on cycle limit wall=149 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 6.349884e+ m35zams 2026-11-24 16:10:41 rc=0 (srun) ended: Terminating on cycle limit wall=77 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 3.354768e+ 
- Resolved keys (parameter dump by -n of the deck copy), the same in all three: cfl 0.3, vceil 3e8, thinf true, lat_sums true, tie_gid true, opac_drag true, allow_clamp false, time2_gtime_slope true, clean_pc_sweep rb, mg 2 (m35mms 4). Full dumps in caltech_results/PREFLIGHT_1124/<star>/pardump.txt.

## Deck differences from prod_rc4v_1124/decks (every changed line)
- bsg: vceil 1.0e8 -> 3.0e8; thinf false -> true; allow_clamp true -> false; time2_gtime_slope = true ADDED; lat_sums / opac_drag / tie_gid 0 -> true; rst cadence 172800 -> 6.048e5 (7 d); dt_report line removed.
- m35mms: mg 2 -> 4; thinf false -> true; allow_clamp true -> false; gtime_slope added; lat_sums / opac_drag / tie_gid 0 -> true; rst cadence 21600 -> 8.64e4; dt_report removed.
- **M35ZAMS: layout_tile 64 -> 32 and MeshBlocks 64x64 -> 32x32 (16 blocks of 384x32x32, 4 per GPU), i.e. NOT the Z3 layout (4 blocks of 384x64x64) the USER switched to on 11-24; thinf false -> true; gtime_slope added; lat_sums / opac_drag / tie_gid 0 -> true; rst cadence 21600 -> 8.64e4; dt_report removed.** Cost: 0.357 s/cyc here against 0.283 (candidate fresh on the Z3 layout, earlier today): +26 %. Is the 16-block layout intended? Say before any start.
- thinf is true in all four decks (streaming form on), including the M35ZAMS wallsym deck; on t9l62 it is the key that makes the init warn (NOTE zzzz5).
