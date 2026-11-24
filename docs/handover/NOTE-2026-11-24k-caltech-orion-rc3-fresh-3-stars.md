# NOTE 2026-11-24k caltech -> orion: IA12.3-rc3 af836f5a3 build + gate, GPU key tests, FRESH bsg / m35mms / t9l62 production (your 24e/24h/24j)

## Build (11-24 00:25)
```
VET_BUILD_OK 00:23:11 md5=09f50abee4ffa6cd399ce15879626353 TUs=53 wall=194 s (incremental, from 8197e0c4f)
BIP_BUILD_OK 00:24:52 md5=8dc36df51c211522cb96f26f1c25edcb wall=295 s (cold, shared FS)
GATE FAIL (seam+wall sym deck copy): athena_ia123rc3_af836f5a3 md5 09f50abe node hpc-sm-02-17 pass1=
GATE PASS (seam+wall sym deck copy): athena_ia123rc3_af836f5a3 md5 09f50abe node hpc-sm-02-17 pass1=3.686e-11
<rad_vet> clean_init steady pass 0: max(|dE|, |dF|/c)/max E = 9.534e-02
<rad_vet> clean_init steady pass 1: max(|dE|, |dF|/c)/max E = 3.686e-11
<rad_vet> clean_init steady pass 2: max(|dE|, |dF|/c)/max E = 9.994e-15
```
- Recipe: the usual in-place incremental build (v203/build_n, PROBLEM=he_star_vet; from 8197e0c4f) + a separate build dir v203/build_bip for PROBLEM=built_in_pgens (cold, fresh cmake, shared FS). The 24i /dev/shm + ninja recipe is NOT used yet (setup cost), so there are no 24i times; I can set it up between arms.
- Gate deck: bsg_wngl_smoke copy + vet_gd_theta_seam_sym + vet_gd_wall_sym (as 23zz24).

## Fresh production starts (11-24 00:26)
- bsg: `2026-11-24 00:25:57 chunk 0 STAR=bsg BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_bsg.athinput md5=9f4884524a7d230edd2fc32eb63525fd start=-i /resnick/gro`
- t9l62: `2026-11-24 00:25:57 chunk 0 STAR=t9l62 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_t9l62.athinput md5=194a66017c67d98275cee4c9f9af0815 start=-i /resnick`
- m35mms: `2026-11-24 00:25:57 chunk 0 STAR=m35mms BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_m35mms.athinput md5=ae635eee87acbd5436e2c61f058453e6 start=-i /resni`
- m35zams: `2026-11-24 00:25:57 chunk 0 STAR=m35zams BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_m35zams.athinput md5=f29145b5c76d3158f42f68a5709f9009 start=-i /res`
- Launcher prod_rc3_1124/launch_rc3.sh -> prodg.sh (1 node, 4 ranks, romio321, pmix; KEEPRST 3 only in prod_rc3_1124/<star>/rst). END = holder end - 20 min (7333072: 11-26 05:01; 7333073: 11-26 17:34; 7333036: 11-26 23:07). m35zams (USER OK, 24l): deck ia123rc3_decks/IA123_m35zams_wallsym.athinput as is (md5 f29145b5), on hpc-sm-01-15 (7333036). Tile decks ia123rc3_decks/tile (md5 OK), outputs as in the decks (m35mms hydro_w / m1 1800 s).

## m35zams (11-24 00:27) 
- start: `2026-11-24 00:25:57 chunk 0 STAR=m35zams BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_m35zams.athinput md5=f29145b5c76d3158f42f68a5709f9009 start=-i /resnick/groups/carnegie_poc/jingz`
- FATAL 0, CLEAN WARNING 0, stage fallbacks 0, DT_RETRY 0, created inf 0; never-read: ### mesh/cs_seam_blockend = extend;### mesh/cs_vertex_3cell = true;### mesh/cs_seam_rho_guard = 4;### time/dt_retry_max = 8;### hydro/wb_cache_every = 0;### rad_vet/clean_stag_thinf_held = false;### rad_vet/time2_stage_local_be = false;### rad_vet/vet_gd_rotate_align = true;### rad_vet/vet_rotate_av
- s/cyc 0.3397 (c10-240); row 200 (cycle 200): t 6.961014e+02 dt 3.48051 L_top/L_in 1.000463 L_bot/L_in 0.999890
- md5 of the first 201 non-# rows (cycles 0-200): user.hst 389b393268498febe19f29b63adeefd1, hydro.hst 80043f8487b2a2ba33dea9a71eb7d49a (rows in caltech_results/PROD_RC3_1124/m35zams/user_hst_first201rows.txt)

## m35mms (11-24 00:28) 
- start: `2026-11-24 00:25:57 chunk 0 STAR=m35mms BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_m35mms.athinput md5=ae635eee87acbd5436e2c61f058453e6 start=-i /resnick/groups/carnegie_poc/jingze/`
- FATAL 0, CLEAN WARNING 0, stage fallbacks 0, DT_RETRY 0, created inf 0; never-read: ### mesh/cs_seam_blockend = extend;### mesh/cs_vertex_3cell = true;### mesh/cs_seam_rho_guard = 4;### time/dt_retry_max = 8;### hydro/wb_cache_every = 0;### rad_vet/clean_stag_thinf_held = false;### rad_vet/time2_stage_local_be = false;### rad_vet/vet_gd_rotate_align = true;### rad_vet/vet_rotate_av
- s/cyc 0.6174 (c50-200); row 200 (cycle 200): t 2.401526e+03 dt 12.0076 L_top/L_in 1.001077 L_bot/L_in 0.999344
- md5 of the first 201 non-# rows (cycles 0-200): user.hst b1d5579c37ff27ec25fdb174fa27baf2, hydro.hst bfdd1f3ae0e7fd6066faff5d44df98e3 (rows in caltech_results/PROD_RC3_1124/m35mms/user_hst_first201rows.txt)
- (COORD stopped the IA12.2 B m35mms production on 02-03 at 00:30 after this pass.)

## bsg (11-24 00:28) 
- start: `2026-11-24 00:25:57 chunk 0 STAR=bsg BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_bsg.athinput md5=9f4884524a7d230edd2fc32eb63525fd start=-i /resnick/groups/carnegie_poc/jingze/prod_r`
- FATAL 0, CLEAN WARNING 0, stage fallbacks 0, DT_RETRY 0, created inf 0; never-read: ### mesh/cs_seam_blockend = extend;### mesh/cs_vertex_3cell = true;### mesh/cs_seam_rho_guard = 4;### time/dt_retry_max = 8;### hydro/merge_lat_mom = rigid;### hydro/wb_cache_every = 0;### rad_vet/clean_stag_thinf_held = false;### rad_vet/time2_stage_local_be = false;### rad_vet/vet_gd_rotate_align 
- s/cyc 0.6644 (c50-200); row 200 (cycle 200): t 3.359495e+04 dt 168.721 L_top/L_in 1.054684 L_bot/L_in 1.534098
- md5 of the first 201 non-# rows (cycles 0-200): user.hst f907df1f36dc5a8d11439e43ab157a5c, hydro.hst 446533081340295d2a239eb200762c9b (rows in caltech_results/PROD_RC3_1124/bsg/user_hst_first201rows.txt)

## GPU key tests, 1 H200 (hpc-sm-02-17), built_in_pgens binary (11-24 00:29)
- BIP_BUILD_OK 00:24:52 md5=8dc36df51c211522cb96f26f1c25edcb wall=295 s (cold, shared FS)
```
clean VET test suite: 8 tests, full, 4 at a time x 1 threads, binary /resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3bip_af836f5a3
PASS restart_sp           47 s  files 4 maxdiff 0
PASS dt_retry             59 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS layout_tile         107 s  maxdiff_1x2 0 maxrel_1x2 0 maxdiff_2x2 7.9e-12 maxrel_2x2 3.11e-10 maxdiff_dfl 0 dE_dfl 2.22e-16 dM_dfl 1.8e-16 dLz_dfl 1.96e-16
PASS conservation_lz     110 s  sp_dE 2.22e-16 sp_dM 8.99e-16 sp_dLz 1.96e-16 sp_dLz_gr 3.91e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.24e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS floor_pay            54 s  dE 2.22e-16 dE0 2.22e-16 paid 0.286 unpaid 0 warnings 0 worst_lines 196 worst_bitwise True ekin_dE 1.11e-16 ekin_dM 7.19e-16 ekin_dLz 1.96e-16
PASS inf_reject            5 s  dE 2.22e-16 nonfinite 0 rejected 1
PASS wedge_sym            72 s  fphi [1.22e-09 1.18e-10 8.25e-10] warnings 0 wedge_fatal 1 fth 8.06e-15 theta_fatal 1
PASS ceiling_metric       14 s  sp_dvphi 0 sp_dcap 0 sp_vtheta 0 sp_scale_ok 1 sp_capped 1408 cs_dvphi 6.31e-16 cs_dcap 6.31e-16 cs_vtheta 1.5e-16 cs_scale_ok 1 cs_capped 2928
8 passed, 0 failed, 121 s
```
- Full output + record: caltech_results/RC3_KEYTESTS/.

## m35zams thin-layer |F_phi|/|F_r| (11-24 00:33)
- No vet bin at cycle 200 (t 696; bins every 1800 s). These are the nearest dumps: t=0 (00000) and t=1800.02 / cycle 601 (00001). Float32 bins, vet_f1/f2/f3 cell moments; thin layer = cells with rho 1e-16..4e-13 (the band of the 23zt v_phi growth, r ~5.53-5.55e11).
```
time 0 cycle 0 (vet bin 00000; float32 output)
thin layer rho 1e-16..4e-13 (cells): n=196608  max|F_phi|/|F_r| 1.833e-02  p99 1.443e-02  median 6.698e-03 | max|F_theta|/|F_r| 2.782e-02
rho 1e-14..4e-13 (cells): n=98304  max|F_phi|/|F_r| 1.833e-02  p99 1.509e-02  median 8.129e-03 | max|F_theta|/|F_r| 2.782e-02
shell ratio sum|F_phi|/sum|F_r| over the thin-layer shells i=320..331 (r 5.5292e+11..5.5590e+11): max 8.176e-03

time 1800.02 cycle 601 (vet bin 00001; float32 output)
thin layer rho 1e-16..4e-13 (cells): n=131072  max|F_phi|/|F_r| 1.996e-02  p99 1.470e-02  median 7.282e-03 | max|F_theta|/|F_r| 2.474e-02
rho 1e-14..4e-13 (cells): n=81920  max|F_phi|/|F_r| 1.996e-02  p99 1.506e-02  median 8.058e-03 | max|F_theta|/|F_r| 2.474e-02
shell ratio sum|F_phi|/sum|F_r| over the thin-layer shells i=321..328 (r 5.5319e+11..5.5509e+11): max 8.093e-03
```
- For reference, IA12.1 smoke 23zt (S22 2x2 64^2 layout, non-wallsym deck) at t 3.6e3 already had a coherent v_phi of 3.9e6 cm/s there. Script: caltech_results/PROD_RC3_1124/m35zams/fphi.py.

## t9l62 (11-24 00:35) 
- start: `2026-11-24 00:25:57 chunk 0 STAR=t9l62 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_t9l62.athinput md5=194a66017c67d98275cee4c9f9af0815 start=-i /resnick/groups/carnegie_poc/jingze/pr`
- FATAL 0, CLEAN WARNING 0, stage fallbacks 0, DT_RETRY 0, created inf 0; never-read: ### mesh/cs_seam_blockend = extend;### mesh/cs_vertex_3cell = true;### mesh/cs_seam_rho_guard = 4;### time/dt_retry_max = 8;### hydro/wb_cache_every = 0;### rad_vet/clean_stag_thinf_held = false;### rad_vet/time2_stage_local_be = false;### rad_vet/vet_gd_rotate_align = true;### rad_vet/vet_rotate_av
- s/cyc 2.3070 (c50-200); row 200 (cycle 200): t 1.937038e+04 dt 96.8505 L_top/L_in 1.002203 L_bot/L_in 0.938786
- md5 of the first 201 non-# rows (cycles 0-200): user.hst c6f22bb0fab6105f8be69389d0cd09ee, hydro.hst 694e3bcfadb96b67ba8f2af07fadbcd5 (rows in caltech_results/PROD_RC3_1124/t9l62/user_hst_first201rows.txt)

## m35zams thin layer, signed sums (orion 24o) at 11-24 00:42: 4 vet dumps
```
 dump         t  cycle | sumFphiV/|Fr|V sumFthV/|Fr|V |     <v_phi>_m    max|v_phi| |       <v_r>_m     <|v_r|>_m
00000         0      0 |    1.0897e-11   -3.4790e-14 |    0.0000e+00    0.0000e+00 |    0.0000e+00    0.0000e+00
00001      1800    601 |    1.1295e-04   -4.0297e-08 |   -4.6006e+02    3.0000e+07 |    2.5043e+06    2.5048e+06
00002      3600   1320 |   -1.5021e-04   -1.9715e-07 |    3.7450e+03    2.2801e+05 |    1.9148e+06    2.4344e+06
00003      5401   2041 |   -1.5355e-03   -2.2347e-06 |   -1.9721e+03    6.4105e+05 |   -3.9430e+06    4.0677e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
- Volume-weighted over shells i=320..331 (the 23zt thin layer). F_phi / F_theta are SIGNED sums over the |F_r| sum; v_phi and v_r are mass-weighted means; max|v_phi| is over the cells. Script: caltech_results/PROD_RC3_1124/m35zams/thin.py (srun CPU step on 01-15).

- (11-24 00:42) **bsg STOPPED (USR1 = graceful rst + exit) on:  CW created-inf**, step 7333072.82, log.0:
```
<hydro> floor_audit cycle 777 pay: paid 2.820596039e+39 unpaid 0 given_to_rad 2.396834386e+23 cells 34730701 | fofc_retest: retested 0 unflagged 0
CLEAN cycle 777 t2s 2 it 0 (sum 7) pos 0 mark 0 all 0 neg 379 r1 379 r2 1241 first 379 1241 guard 0 solves 6 cwin 0
CLEAN WARNING XT cycle 777 t2s 2 write-back closure inf
CLEAN WARNING cycle 777 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 1)
CLEAN cycle 777 t2s 1 it 1 (sum 7) pos 0 mark 0 all 0 neg 0 r1 0 r2 187 first 0 187 guard 0 solves 5 cwin 0
CLEAN FALLBACK cycle 777: pos_gas 8 (from E 4.300e+31) e_floor 0 (from gas 0.000e+00, column 0.000e+00) ffct_left 0 auto_be 0 floorcap v 0 T 0 c2p_dfloor 91488 (mass 1.024e+07) floor_pay 3.961e+38 (unpaid 0.000e+00)
<hydro> floor_audit cycle 778 [update] cell-c2p 2.61053481e+10 e_c<=0 599876 kedom 1481599 sub_tfloor 615829 sub_edge 0 dual 0 dual_de 0 dual_sub_tfloor 0 floored 693417 (kedom 659702) floor_de 4.568888467e+38 <r> 3.978187927e+12 <r>_de 3.989411779e+
--
<hydro> floor_audit cycle 779 pay: paid 3.67640457e+39 unpaid 0 given_to_rad 2.40745216e+23 cells 34781317 | fofc_retest: retested 0 unflagged 0
CLEAN cycle 779 t2s 2 it 0 (sum 10) pos 0 mark 0 all 0 neg 615 r1 615 r2 1944 first 615 1944 guard 0 solves 6 cwin 0
CLEAN WARNING XT cycle 779 t2s 2 write-back closure inf
CLEAN WARNING cycle 779 t2s 2: implicit_pos_floor created 1.459516e+42 of E (neither the gas nor the column covered the raise)
CLEAN WARNING cycle 779 t2s 2: stage solve failed, the step is redone with backward Euler (stage fallback 2)
CLEAN WARNING cycle 779 t2s 2: the implicit_pos_floor creation 1.459516e+42 of the failed attempt is discarded with it
```

## s/cyc before / after FOFC (orion 24q; COORD relaunched all 4 with hydro/fofc=true hydro/fofc_retest=true; log.0 = no FOFC, log.1 = FOFC)
| star | log.0 (no FOFC) | log.1 (FOFC) | change |
|---|---|---|---|
| bsg | 0.673 (c50-150), 0.759 (c1050-1300, toward the bad 00002; dt fell to 8.9e-6) | 0.657 (c600-700, first after the restart from 00001), 0.532 (c1900-2150) | early window -2 %; late windows not comparable (log.0 was collapsing) |
| t9l62 | 2.303 (c50-150), 2.344 (c150-400) | 2.386 (c500-600), 2.412 (c550-800) | +3 % |
| m35mms | 0.605 (c50-150), 0.757 (c1600-1850) | 0.595 (c1950-2050), 0.682 (c3150-3400) | about -2 % at the start of each log; different phases later |
| m35zams | 0.340 (c3410-3460) | 0.370 (c3480-3500), 0.366 (c6370-6420) | +8 % |

- From the 50-cycle (m35zams: 10-cycle) elapsed prints; windows as listed. The runs are on different cycles, so this is a rough comparison. dt is unchanged (m35mms 12.0076; t9l62 96.85; m35zams 2.49).
- Watch re-armed (prod_rc3_1124/prodwatch2.sh, 5 min): it posts every CLEAN WARNING with context (no stop). On FATAL / 'created inf' it posts and COORD asks the USER; only m35mms 'created inf' writes STOP at once (USER pre-approved). Status every 2 h.

## m35zams thin layer, signed sums (orion 24o) at 11-24 01:13: 11 vet dumps
```
00005      9000   3485 |    7.1752e-04    1.9220e-07 |    3.1269e+04    2.1059e+06 |   -1.5987e+06    1.6334e+06
00006  1.08e+04   4203 |    4.8903e-04   -3.7440e-07 |    9.4914e+03    5.0767e+06 |   -9.8358e+05    9.8358e+05
00007  1.26e+04   4927 |    2.6021e-04   -9.2387e-09 |   -6.1489e+03    5.5981e+06 |   -5.6087e+05    5.6183e+05
00008  1.44e+04   5653 |    6.7599e-04   -2.5037e-07 |   -2.2011e+04    4.8252e+06 |   -4.9114e+05    4.9206e+05
00009  1.62e+04   6378 |    4.8321e-04   -1.7067e-07 |   -1.9673e+04    4.1902e+06 |   -2.6207e+05    2.6863e+05
00010   1.8e+04   7104 |    7.1289e-04   -4.7976e-08 |   -3.9166e+04    2.7038e+06 |    1.1445e+05    1.5089e+05
00011  1.98e+04   7830 |    1.6817e-04   -1.5113e-07 |   -4.8355e+04    4.8276e+06 |    2.7961e+05    2.8401e+05
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams thin layer, signed sums (orion 24o) at 11-24 01:45: 16 vet dumps
```
00009  1.62e+04   6378 |    4.8321e-04   -1.7067e-07 |   -1.9673e+04    4.1902e+06 |   -2.6207e+05    2.6863e+05
00010   1.8e+04   7104 |    7.1289e-04   -4.7976e-08 |   -3.9166e+04    2.7038e+06 |    1.1445e+05    1.5089e+05
00011  1.98e+04   7830 |    1.6817e-04   -1.5113e-07 |   -4.8355e+04    4.8276e+06 |    2.7961e+05    2.8401e+05
00012  2.16e+04   8555 |    4.4789e-04    5.8762e-08 |   -5.9545e+04    1.3885e+07 |    3.8050e+05    3.8290e+05
00013  2.34e+04   9281 |    5.3270e-04    1.3909e-07 |   -7.2101e+04    9.4772e+06 |    1.5443e+05    1.6704e+05
00014  2.52e+04  10008 |    4.8077e-04    1.0076e-07 |   -7.5365e+04    7.7836e+06 |    4.6922e+04    8.8922e+04
00015 2.619e+04  10408 |    6.7870e-04   -1.3452e-07 |   -8.2293e+04    2.9810e+06 |   -1.9709e+04    9.2487e+04
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams_z3 (11-24 01:57) 
- start: `2026-11-24 01:55:47 chunk 0 STAR=m35zams_z3 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/layopt/decks/IA123_m35zams_lay.athinput md5=6811bf5735a0d56d6ebf8c2860bccec1 start=-i /resnick/groups/carnegie_`
- FATAL 0, CLEAN WARNING 0, stage fallbacks 0, DT_RETRY 0, created inf 0; never-read: ### mesh/cs_seam_blockend = extend;### mesh/cs_vertex_3cell = true;### mesh/cs_seam_rho_guard = 4;### time/dt_retry_max = 8;### hydro/wb_cache_every = 0;### rad_vet/clean_stag_thinf_held = false;### rad_vet/time2_stage_local_be = false;### rad_vet/vet_gd_rotate_align = true;### rad_vet/vet_rotate_av
- s/cyc 0.2661 (c10-280); row 200 (cycle 200): t 6.961014e+02 dt 3.48051 L_top/L_in 1.000460 L_bot/L_in 0.999890
- md5 of the first 201 non-# rows (cycles 0-200): user.hst 1e0c57b57e944bf7f534fe471424f9d2, hydro.hst 7e8b70e0f8f80341fba8ebbf18fad8e5 (rows in caltech_results/PROD_RC3_1124/m35zams_z3/user_hst_first201rows.txt)

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 01:59: 0 vet dumps
```
 dump         t  cycle | sumFphiV/|Fr|V sumFthV/|Fr|V |     <v_phi>_m    max|v_phi| |       <v_r>_m     <|v_r|>_m
00000         0      0 |    1.1650e-11    1.5311e-14 |    0.0000e+00    0.0000e+00 |    0.0000e+00    0.0000e+00
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
- Volume-weighted over shells i=320..331 (the 23zt thin layer). F_phi / F_theta are SIGNED sums over the |F_r| sum; v_phi and v_r are mass-weighted means; max|v_phi| is over the cells. Script: caltech_results/PROD_RC3_1124/m35zams/thin.py (srun CPU step on 01-15).

## m35zams_z3 (11-24 01:59) 
- start: `2026-11-24 01:58:05 chunk 0 STAR=m35zams_z3 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/layopt/decks/IA123_m35zams_lay_lz.athinput md5=5a4db7fb77013d8f546781bc76bd605d start=-i /resnick/groups/carneg`
- FATAL 0, CLEAN WARNING 0, stage fallbacks 0, DT_RETRY 0, created inf 0; never-read: ### mesh/cs_seam_blockend = extend;### mesh/cs_vertex_3cell = true;### mesh/cs_seam_rho_guard = 4;### time/dt_retry_max = 8;### hydro/wb_cache_every = 0;### rad_vet/clean_stag_thinf_held = false;### rad_vet/time2_stage_local_be = false;### rad_vet/vet_gd_rotate_align = true;### rad_vet/vet_rotate_av
- s/cyc 0.2647 (c10-310); row 200 (cycle 200): t 6.961014e+02 dt 3.48051 L_top/L_in 1.000460 L_bot/L_in 0.999890
- md5 of the first 201 non-# rows (cycles 0-200): user.hst 1e0c57b57e944bf7f534fe471424f9d2, hydro.hst 7e8b70e0f8f80341fba8ebbf18fad8e5 (rows in caltech_results/PROD_RC3_1124/m35zams_z3/user_hst_first201rows.txt)

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 02:02: 2 vet dumps
```
 dump         t  cycle | sumFphiV/|Fr|V sumFthV/|Fr|V |     <v_phi>_m    max|v_phi| |       <v_r>_m     <|v_r|>_m
00000         0      0 |    1.1650e-11    1.5311e-14 |    0.0000e+00    0.0000e+00 |    0.0000e+00    0.0000e+00
00001      1802    525 |    1.2257e-04   -4.0754e-08 |   -4.6103e+02    3.0732e+05 |    2.4903e+06    2.4905e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
- Volume-weighted over shells i=320..331 (the 23zt thin layer). F_phi / F_theta are SIGNED sums over the |F_r| sum; v_phi and v_r are mass-weighted means; max|v_phi| is over the cells. Script: caltech_results/PROD_RC3_1124/m35zams/thin.py (srun CPU step on 01-15).
- (11-24 02:10) m35zams_z3 dt per 1000 cycles: c1000 t 3.020044e+03 dt 2.503296e+00 (wall 2.803689e+02 s);c2000 t 5.521967e+03 dt 2.498548e+00 (wall 5.472784e+02 s);
- (11-24 02:41) m35zams_z3 dt per 1000 cycles: c3000 t 8.014082e+03 dt 2.491558e+00 (wall 8.226977e+02 s);c4000 t 1.050364e+04 dt 2.487109e+00 (wall 1.101079e+03 s);c5000 t 1.281244e+04 dt 2.029738e+00 (wall 1.375524e+03 s);c6000 t 1.492077e+04 dt 1.751397e+00 (wall 1.647705e+03 s);c7000 t 1.717329e+04 dt 2.479212e+00 (wall 1.926822e+03 s);c8000 t 1.965230e+04 dt 2.480016e+00 (wall 2.217261e+03 s);c9000 t 2.213235e+04 dt 2.479691e+00 (wall 2.524838e+03 s);

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 03:00: 18 vet dumps
```
00007  1.26e+04   4908 |    2.7728e-04   -1.1194e-08 |   -5.7812e+03    6.4861e+06 |   -5.4839e+05    5.4990e+05
00008  1.44e+04   5746 |    7.0642e-04   -2.5115e-07 |   -2.1968e+04    6.0379e+06 |   -4.8675e+05    4.8788e+05
00009  1.62e+04   6603 |    5.0484e-04   -1.7026e-07 |   -1.9859e+04    4.8381e+06 |   -2.6539e+05    2.7162e+05
00010   1.8e+04   7334 |    7.4581e-04   -3.9750e-08 |   -3.9508e+04    2.9522e+06 |    1.0864e+05    1.4717e+05
00011  1.98e+04   8060 |    2.1889e-04   -1.3988e-07 |   -4.8738e+04    5.5133e+06 |    2.7454e+05    2.7915e+05
00012  2.16e+04   8786 |    5.1058e-04    6.5663e-08 |   -5.9933e+04    1.4293e+07 |    3.7576e+05    3.7826e+05
00013  2.34e+04   9513 |    6.0035e-04    1.3899e-07 |   -7.2437e+04    1.4793e+07 |    1.5352e+05    1.6627e+05
00014  2.52e+04  10244 |    5.6557e-04    8.8849e-08 |   -7.5332e+04    3.1623e+07 |    4.9082e+04    9.0248e+04
00015   2.7e+04  10975 |    5.6671e-04    8.1218e-08 |   -9.2537e+04    7.3910e+06 |   -1.3842e+05    1.4653e+05
00016  2.88e+04  11713 |    1.1395e-03   -9.2259e-08 |   -1.0714e+05    1.6875e+07 |   -1.1054e+05    1.3759e+05
00017  3.06e+04  12440 |    1.3947e-03   -1.2365e-07 |   -1.2783e+05    6.1348e+06 |   -1.1113e+05    1.3168e+05
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 03:04: 19 vet dumps
```
00008  1.44e+04   5746 |    7.0642e-04   -2.5115e-07 |   -2.1968e+04    6.0379e+06 |   -4.8675e+05    4.8788e+05
00009  1.62e+04   6603 |    5.0484e-04   -1.7026e-07 |   -1.9859e+04    4.8381e+06 |   -2.6539e+05    2.7162e+05
00010   1.8e+04   7334 |    7.4581e-04   -3.9750e-08 |   -3.9508e+04    2.9522e+06 |    1.0864e+05    1.4717e+05
00011  1.98e+04   8060 |    2.1889e-04   -1.3988e-07 |   -4.8738e+04    5.5133e+06 |    2.7454e+05    2.7915e+05
00012  2.16e+04   8786 |    5.1058e-04    6.5663e-08 |   -5.9933e+04    1.4293e+07 |    3.7576e+05    3.7826e+05
00013  2.34e+04   9513 |    6.0035e-04    1.3899e-07 |   -7.2437e+04    1.4793e+07 |    1.5352e+05    1.6627e+05
00014  2.52e+04  10244 |    5.6557e-04    8.8849e-08 |   -7.5332e+04    3.1623e+07 |    4.9082e+04    9.0248e+04
00015   2.7e+04  10975 |    5.6671e-04    8.1218e-08 |   -9.2537e+04    7.3910e+06 |   -1.3842e+05    1.4653e+05
00016  2.88e+04  11713 |    1.1395e-03   -9.2259e-08 |   -1.0714e+05    1.6875e+07 |   -1.1054e+05    1.3759e+05
00017  3.06e+04  12440 |    1.3947e-03   -1.2365e-07 |   -1.2783e+05    6.1348e+06 |   -1.1113e+05    1.3168e+05
00018  3.24e+04  13167 |    1.9584e-03   -2.0321e-07 |   -1.5985e+05    7.8233e+06 |   -5.1163e+04    9.3972e+04
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
- (11-24 03:11) m35zams_z3 dt per 1000 cycles: c10000 t 2.460643e+04 dt 2.478854e+00 (wall 2.824096e+03 s);c11000 t 2.706197e+04 dt 2.369481e+00 (wall 3.127529e+03 s);c12000 t 2.951167e+04 dt 2.478231e+00 (wall 3.429198e+03 s);c13000 t 3.198850e+04 dt 2.478397e+00 (wall 3.726071e+03 s);c14000 t 3.446340e+04 dt 2.478232e+00 (wall 4.021253e+03 s);c15000 t 3.694098e+04 dt 2.477430e+00 (wall 4.310705e+03 s);
- (11-24 03:41) m35zams_z3 dt per 1000 cycles: c16000 t 3.941773e+04 dt 2.477433e+00 (wall 4.599349e+03 s);c17000 t 4.183726e+04 dt 2.477388e+00 (wall 4.886294e+03 s);c18000 t 4.431469e+04 dt 2.477510e+00 (wall 5.178268e+03 s);c19000 t 4.674285e+04 dt 2.478028e+00 (wall 5.483689e+03 s);c20000 t 4.883177e+04 dt 2.031900e+00 (wall 5.762350e+03 s);c21000 t 5.079662e+04 dt 1.988411e+00 (wall 6.042611e+03 s);

- (11-24 03:46) RSG status: log.1 cycle=776600 time=7.267617e+08 dt=9.166667e+02, CW 0, DT_RETRY 0

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 04:02: 34 vet dumps
```
00023  4.14e+04  16824 |   -1.3506e-04   -9.7258e-07 |   -3.4826e+05    8.3834e+06 |    1.8269e+02    6.8996e+04
00024  4.32e+04  17551 |   -1.9640e-03   -8.1356e-06 |   -4.0281e+05    1.1119e+07 |   -1.0651e+04    6.8357e+04
00025   4.5e+04  18277 |   -4.9829e-03   -1.5081e-05 |   -4.6216e+05    2.3433e+07 |   -1.6582e+04    6.6731e+04
00026  4.68e+04  19024 |   -9.8300e-03   -2.4698e-05 |   -5.3198e+05    1.6315e+07 |   -5.3704e+03    6.4789e+04
00027  4.86e+04  19883 |   -1.7411e-02   -3.7062e-05 |   -5.9712e+05    2.9909e+07 |   -8.6471e+03    6.2472e+04
00028  5.04e+04  20795 |   -2.9517e-02   -5.5786e-05 |   -6.8223e+05    3.6876e+07 |   -5.3587e+03    6.1985e+04
00029  5.22e+04  21711 |   -4.5665e-02   -7.4805e-05 |   -7.7048e+05    3.8301e+07 |    1.1087e+03    6.2355e+04
00030   5.4e+04  22635 |   -6.4908e-02   -7.1257e-05 |   -8.6075e+05    3.7824e+07 |    5.4585e+03    6.2311e+04
00031  5.58e+04  23585 |   -8.5154e-02   -6.2605e-05 |   -9.5757e+05    3.1832e+07 |   -1.1319e+02    6.3573e+04
00032  5.76e+04  24523 |   -1.0561e-01   -4.7191e-05 |   -1.0471e+06    3.2316e+07 |   -3.4366e+03    6.2300e+04
00033  5.94e+04  25468 |   -1.2477e-01   -2.7321e-05 |   -1.1158e+06    2.9919e+07 |   -2.0833e+03    6.3493e+04
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 04:05: 34 vet dumps
```
00024  4.32e+04  17551 |   -1.9640e-03   -8.1356e-06 |   -4.0281e+05    1.1119e+07 |   -1.0651e+04    6.8357e+04
00025   4.5e+04  18277 |   -4.9829e-03   -1.5081e-05 |   -4.6216e+05    2.3433e+07 |   -1.6582e+04    6.6731e+04
00026  4.68e+04  19024 |   -9.8300e-03   -2.4698e-05 |   -5.3198e+05    1.6315e+07 |   -5.3704e+03    6.4789e+04
00027  4.86e+04  19883 |   -1.7411e-02   -3.7062e-05 |   -5.9712e+05    2.9909e+07 |   -8.6471e+03    6.2472e+04
00028  5.04e+04  20795 |   -2.9517e-02   -5.5786e-05 |   -6.8223e+05    3.6876e+07 |   -5.3587e+03    6.1985e+04
00029  5.22e+04  21711 |   -4.5665e-02   -7.4805e-05 |   -7.7048e+05    3.8301e+07 |    1.1087e+03    6.2355e+04
00030   5.4e+04  22635 |   -6.4908e-02   -7.1257e-05 |   -8.6075e+05    3.7824e+07 |    5.4585e+03    6.2311e+04
00031  5.58e+04  23585 |   -8.5154e-02   -6.2605e-05 |   -9.5757e+05    3.1832e+07 |   -1.1319e+02    6.3573e+04
00032  5.76e+04  24523 |   -1.0561e-01   -4.7191e-05 |   -1.0471e+06    3.2316e+07 |   -3.4366e+03    6.2300e+04
00033  5.94e+04  25468 |   -1.2477e-01   -2.7321e-05 |   -1.1158e+06    2.9919e+07 |   -2.0833e+03    6.3493e+04
00034  6.12e+04  26423 |   -1.3973e-01   -4.7284e-05 |   -1.1944e+06    3.1248e+07 |   -1.2250e+04    6.1401e+04
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
- (11-24 04:11) m35zams_z3 dt per 1000 cycles: c22000 t 5.276485e+04 dt 1.845216e+00 (wall 6.321371e+03 s);c23000 t 5.469825e+04 dt 1.805268e+00 (wall 6.600154e+03 s);c24000 t 5.659522e+04 dt 1.714710e+00 (wall 6.879727e+03 s);c25000 t 5.851048e+04 dt 1.900824e+00 (wall 7.162713e+03 s);c26000 t 6.038583e+04 dt 1.833869e+00 (wall 7.448624e+03 s);c27000 t 6.232854e+04 dt 1.945029e+00 (wall 7.738806e+03 s);

- (11-24 04:23) status: bsg [log.2 cycle=21700 time=1.364828e+06 dt=6.009534e+01, CW 0, DT_RETRY 0]; t9l62 [log.2 cycle=4300 time=4.162843e+05 dt=9.090845e+01, CW 0, DT_RETRY 0]; m35mms [log.4 cycle=15200 time=1.330905e+05 dt=1.111823e+01, CW 0, DT_RETRY 0]; m35zams_z3 [log.0 cycle=30300 time=6.842044e+04 dt=1.873972e+00, CW 0, DT_RETRY 0];
- (11-24 04:41) m35zams_z3 dt per 1000 cycles: c28000 t 6.424542e+04 dt 1.752589e+00 (wall 8.030806e+03 s);c29000 t 6.615237e+04 dt 1.945076e+00 (wall 8.334468e+03 s);c30000 t 6.789950e+04 dt 1.885244e+00 (wall 8.623248e+03 s);c31000 t 6.966278e+04 dt 2.001032e+00 (wall 8.907828e+03 s);c32000 t 7.130907e+04 dt 1.441911e+00 (wall 9.189618e+03 s);c33000 t 7.274837e+04 dt 1.413751e+00 (wall 9.458798e+03 s);c34000 t 7.415260e+04 dt 1.422293e+00 (wall 9.734291e+03 s);

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 05:03: 45 vet dumps
```
00035   6.3e+04  27343 |   -1.4394e-01   -7.1680e-05 |   -1.2693e+06    3.8108e+07 |   -7.4715e+03    6.6387e+04
00036  6.48e+04  28289 |   -1.3707e-01   -3.1745e-05 |   -1.2810e+06    2.3825e+07 |    2.1117e+03    6.4637e+04
00037  6.66e+04  29258 |   -1.1881e-01   -2.4660e-05 |   -1.1987e+06    3.0615e+07 |    3.5707e+03    6.8156e+04
00038  6.84e+04  30289 |   -9.5952e-02   -2.4945e-05 |   -9.9492e+05    3.6148e+07 |    1.2888e+04    7.0722e+04
00039  7.02e+04  31307 |   -6.1124e-02   -6.8923e-06 |   -6.6442e+05    4.2480e+07 |   -5.8321e+03    7.4005e+04
00040   7.2e+04  32453 |   -2.0502e-02   -2.8044e-05 |   -3.6937e+05    6.4121e+07 |    7.9024e+04    1.2341e+05
00041  7.38e+04  33742 |    2.1431e-02   -3.0137e-04 |   -2.8893e+05    5.7772e+07 |   -3.3434e+04    1.5222e+05
00042  7.56e+04  35084 |    7.7238e-02   -4.7680e-04 |   -4.0957e+05    4.7021e+07 |   -1.1147e+05    3.0652e+05
00043  7.74e+04  36404 |    1.0243e-01   -1.7389e-04 |   -6.1845e+05    4.4741e+07 |   -3.3259e+04    3.9639e+05
00044  7.92e+04  37771 |   -1.7989e-02   -6.2940e-06 |   -7.1947e+05    1.3832e+07 |   -1.2318e+05    4.3943e+05
00045   8.1e+04  38867 |   -1.0011e-02    5.8863e-04 |   -7.9167e+05    2.0377e+07 |    2.5814e+05    7.9796e+05
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 05:06: 46 vet dumps
```
00035   6.3e+04  27343 |   -1.4394e-01   -7.1680e-05 |   -1.2693e+06    3.8108e+07 |   -7.4715e+03    6.6387e+04
00036  6.48e+04  28289 |   -1.3707e-01   -3.1745e-05 |   -1.2810e+06    2.3825e+07 |    2.1117e+03    6.4637e+04
00037  6.66e+04  29258 |   -1.1881e-01   -2.4660e-05 |   -1.1987e+06    3.0615e+07 |    3.5707e+03    6.8156e+04
00038  6.84e+04  30289 |   -9.5952e-02   -2.4945e-05 |   -9.9492e+05    3.6148e+07 |    1.2888e+04    7.0722e+04
00039  7.02e+04  31307 |   -6.1124e-02   -6.8923e-06 |   -6.6442e+05    4.2480e+07 |   -5.8321e+03    7.4005e+04
00040   7.2e+04  32453 |   -2.0502e-02   -2.8044e-05 |   -3.6937e+05    6.4121e+07 |    7.9024e+04    1.2341e+05
00041  7.38e+04  33742 |    2.1431e-02   -3.0137e-04 |   -2.8893e+05    5.7772e+07 |   -3.3434e+04    1.5222e+05
00042  7.56e+04  35084 |    7.7238e-02   -4.7680e-04 |   -4.0957e+05    4.7021e+07 |   -1.1147e+05    3.0652e+05
00043  7.74e+04  36404 |    1.0243e-01   -1.7389e-04 |   -6.1845e+05    4.4741e+07 |   -3.3259e+04    3.9639e+05
00044  7.92e+04  37771 |   -1.7989e-02   -6.2940e-06 |   -7.1947e+05    1.3832e+07 |   -1.2318e+05    4.3943e+05
00045   8.1e+04  38867 |   -1.0011e-02    5.8863e-04 |   -7.9167e+05    2.0377e+07 |    2.5814e+05    7.9796e+05
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

- (11-24 05:08) **t9l62 FATAL** (log.3); clean STOP written + rsts hard-linked to rst_keep (COORD permission, USER asleep):
```
[LOG_CAT_SBGP] Failed to dlopen libnuma.so. Fallback to GROUP_BY_SOCKET manual.
[LOG_CAT_SBGP] libnuma.so: cannot open shared object file: No such file or directory
[LOG_CAT_SBGP] Failed to dlopen libnuma.so. Fallback to GROUP_BY_SOCKET manual.
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/src/parameter_input.cpp at line 460
Parameter 'fofc_retest' in block 'hydro' on command line not found
Kokkos::Cuda ERROR: Failed to call Kokkos::Cuda::finalize()
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/src/parameter_input.cpp at line 460
Parameter 'fofc_retest' in block 'hydro' on command line not found
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/v203/src/parameter_input.cpp at line 460
Parameter 'fofc_retest' in block 'hydro' on command line not found
```
- (11-24 05:12) m35zams_z3 dt per 1000 cycles: c35000 t 7.549501e+04 dt 1.430958e+00 (wall 1.001073e+04 s);c36000 t 7.684122e+04 dt 1.278883e+00 (wall 1.028138e+04 s);c37000 t 7.824394e+04 dt 1.476283e+00 (wall 1.055695e+04 s);c38000 t 7.957407e+04 dt 1.748991e+00 (wall 1.082441e+04 s);c39000 t 8.122855e+04 dt 1.804176e+00 (wall 1.109264e+04 s);c40000 t 8.301159e+04 dt 1.862168e+00 (wall 1.137189e+04 s);
- (11-24 05:42) m35zams_z3 dt per 1000 cycles: c41000 t 8.512593e+04 dt 1.322463e+00 (wall 1.166540e+04 s);c42000 t 8.661554e+04 dt 1.645723e+00 (wall 1.195257e+04 s);c43000 t 8.809211e+04 dt 1.276568e+00 (wall 1.224082e+04 s);c44000 t 8.971747e+04 dt 1.805507e+00 (wall 1.253254e+04 s);c45000 t 9.102972e+04 dt 1.899379e+00 (wall 1.281890e+04 s);c46000 t 9.232048e+04 dt 1.190016e+00 (wall 1.309962e+04 s);c47000 t 9.373723e+04 dt 1.446071e+00 (wall 1.339352e+04 s);

- (11-24 05:47) RSG status: log.1 cycle=781000 time=7.308515e+08 dt=9.565217e+02, CW 0, DT_RETRY 0

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 06:04: 56 vet dumps
```
00046  8.28e+04  39872 |   -2.2421e-02    8.8406e-04 |   -1.9393e+06    1.5090e+07 |   -1.0661e+05    7.7886e+05
00047  8.46e+04  40699 |    8.8149e-03    1.8554e-03 |   -1.7741e+05    1.0989e+07 |    1.1918e+06    1.4808e+06
00048  8.64e+04  41870 |    5.2545e-02    2.6327e-05 |    5.9916e+05    9.1008e+06 |   -1.4969e+06    2.5110e+06
00049  8.82e+04  43089 |    2.6582e-03   -1.1619e-03 |   -1.0153e+06    9.9664e+06 |   -1.4665e+06    1.8209e+06
00050     9e+04  44190 |    1.7188e-02   -1.2273e-03 |   -4.7230e+05    1.2707e+07 |   -1.6077e+05    1.1891e+06
00051  9.18e+04  45622 |   -1.8034e-02    4.2613e-04 |   -6.5111e+05    1.2133e+07 |    8.8820e+05    1.8574e+06
00052  9.36e+04  46902 |   -1.7948e-02   -1.0762e-03 |   -9.5310e+05    1.0789e+07 |    6.2964e+05    1.6839e+06
00053  9.54e+04  48266 |   -1.0553e-02   -4.3803e-04 |   -9.1428e+05    1.0323e+07 |    2.9847e+05    1.4853e+06
00054  9.72e+04  49438 |    1.5288e-02   -2.0951e-03 |    1.3181e+04    1.3456e+07 |    7.0854e+05    1.5516e+06
00055   9.9e+04  50641 |   -3.1458e-02   -4.2084e-04 |    2.6609e+04    9.7496e+06 |   -2.0216e+05    2.2075e+06
00056 1.008e+05  51732 |    1.3359e-02    1.7033e-04 |   -9.8890e+05    9.9247e+06 |   -2.0806e+05    1.5081e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 06:08: 57 vet dumps
```
00046  8.28e+04  39872 |   -2.2421e-02    8.8406e-04 |   -1.9393e+06    1.5090e+07 |   -1.0661e+05    7.7886e+05
00047  8.46e+04  40699 |    8.8149e-03    1.8554e-03 |   -1.7741e+05    1.0989e+07 |    1.1918e+06    1.4808e+06
00048  8.64e+04  41870 |    5.2545e-02    2.6327e-05 |    5.9916e+05    9.1008e+06 |   -1.4969e+06    2.5110e+06
00049  8.82e+04  43089 |    2.6582e-03   -1.1619e-03 |   -1.0153e+06    9.9664e+06 |   -1.4665e+06    1.8209e+06
00050     9e+04  44190 |    1.7188e-02   -1.2273e-03 |   -4.7230e+05    1.2707e+07 |   -1.6077e+05    1.1891e+06
00051  9.18e+04  45622 |   -1.8034e-02    4.2613e-04 |   -6.5111e+05    1.2133e+07 |    8.8820e+05    1.8574e+06
00052  9.36e+04  46902 |   -1.7948e-02   -1.0762e-03 |   -9.5310e+05    1.0789e+07 |    6.2964e+05    1.6839e+06
00053  9.54e+04  48266 |   -1.0553e-02   -4.3803e-04 |   -9.1428e+05    1.0323e+07 |    2.9847e+05    1.4853e+06
00054  9.72e+04  49438 |    1.5288e-02   -2.0951e-03 |    1.3181e+04    1.3456e+07 |    7.0854e+05    1.5516e+06
00055   9.9e+04  50641 |   -3.1458e-02   -4.2084e-04 |    2.6609e+04    9.7496e+06 |   -2.0216e+05    2.2075e+06
00056 1.008e+05  51732 |    1.3359e-02    1.7033e-04 |   -9.8890e+05    9.9247e+06 |   -2.0806e+05    1.5081e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
- (11-24 06:12) m35zams_z3 dt per 1000 cycles: c48000 t 9.510491e+04 dt 9.443427e-01 (wall 1.367600e+04 s);c49000 t 9.654218e+04 dt 1.512507e+00 (wall 1.396141e+04 s);c50000 t 9.792555e+04 dt 1.436974e+00 (wall 1.424888e+04 s);c51000 t 9.962515e+04 dt 1.610053e+00 (wall 1.453304e+04 s);c52000 t 1.012748e+05 dt 1.693847e+00 (wall 1.482363e+04 s);c53000 t 1.027746e+05 dt 1.495586e+00 (wall 1.510895e+04 s);
- (11-24 06:42) m35zams_z3 dt per 1000 cycles: c54000 t 1.041553e+05 dt 1.534326e+00 (wall 1.538963e+04 s);c55000 t 1.056979e+05 dt 1.856368e+00 (wall 1.567789e+04 s);c56000 t 1.070214e+05 dt 1.589882e+00 (wall 1.595505e+04 s);c57000 t 1.085416e+05 dt 1.640182e+00 (wall 1.625370e+04 s);c58000 t 1.101009e+05 dt 1.695892e+00 (wall 1.653417e+04 s);c59000 t 1.114595e+05 dt 9.591549e-01 (wall 1.680752e+04 s);

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 07:06: 66 vet dumps
```
00055   9.9e+04  50641 |   -3.1458e-02   -4.2084e-04 |    2.6609e+04    9.7496e+06 |   -2.0216e+05    2.2075e+06
00056 1.008e+05  51732 |    1.3359e-02    1.7033e-04 |   -9.8890e+05    9.9247e+06 |   -2.0806e+05    1.5081e+06
00057 1.026e+05  52877 |   -1.4030e-02    2.4582e-04 |   -5.0228e+05    1.1461e+07 |   -1.3880e+06    1.8466e+06
00058 1.044e+05  54167 |   -1.5626e-02   -2.9914e-04 |   -3.8919e+05    9.9746e+06 |   -1.0189e+03    1.1504e+06
00059 1.062e+05  55288 |   -3.9251e-03    1.4375e-04 |   -1.3109e+06    1.2013e+07 |   -6.6584e+05    1.3023e+06
00060  1.08e+05  56662 |    2.4675e-02    3.6483e-04 |   -1.3377e+06    1.5177e+07 |   -4.2925e+05    1.7398e+06
00061 1.098e+05  57808 |   -2.5602e-02   -5.5409e-04 |   -6.0754e+05    1.3758e+07 |   -5.4347e+05    1.5175e+06
00062 1.116e+05  59163 |    8.7139e-03   -3.6906e-04 |   -5.3741e+05    1.2655e+07 |   -9.3157e+05    1.9500e+06
00063 1.134e+05  61478 |   -2.9009e-02   -1.7391e-03 |   -8.4898e+05    1.0218e+07 |   -2.3161e+05    1.8683e+06
00064 1.152e+05  62993 |    5.3770e-02   -5.2738e-04 |    4.4741e+05    9.9865e+06 |   -1.4230e+05    1.2481e+06
00065  1.17e+05  64295 |   -3.8670e-02   -3.8995e-04 |   -4.8497e+04    1.1280e+07 |   -1.1743e+05    1.3751e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 07:09: 66 vet dumps
```
00056 1.008e+05  51732 |    1.3359e-02    1.7033e-04 |   -9.8890e+05    9.9247e+06 |   -2.0806e+05    1.5081e+06
00057 1.026e+05  52877 |   -1.4030e-02    2.4582e-04 |   -5.0228e+05    1.1461e+07 |   -1.3880e+06    1.8466e+06
00058 1.044e+05  54167 |   -1.5626e-02   -2.9914e-04 |   -3.8919e+05    9.9746e+06 |   -1.0189e+03    1.1504e+06
00059 1.062e+05  55288 |   -3.9251e-03    1.4375e-04 |   -1.3109e+06    1.2013e+07 |   -6.6584e+05    1.3023e+06
00060  1.08e+05  56662 |    2.4675e-02    3.6483e-04 |   -1.3377e+06    1.5177e+07 |   -4.2925e+05    1.7398e+06
00061 1.098e+05  57808 |   -2.5602e-02   -5.5409e-04 |   -6.0754e+05    1.3758e+07 |   -5.4347e+05    1.5175e+06
00062 1.116e+05  59163 |    8.7139e-03   -3.6906e-04 |   -5.3741e+05    1.2655e+07 |   -9.3157e+05    1.9500e+06
00063 1.134e+05  61478 |   -2.9009e-02   -1.7391e-03 |   -8.4898e+05    1.0218e+07 |   -2.3161e+05    1.8683e+06
00064 1.152e+05  62993 |    5.3770e-02   -5.2738e-04 |    4.4741e+05    9.9865e+06 |   -1.4230e+05    1.2481e+06
00065  1.17e+05  64295 |   -3.8670e-02   -3.8995e-04 |   -4.8497e+04    1.1280e+07 |   -1.1743e+05    1.3751e+06
00066 1.188e+05  65713 |    1.3193e-02   -4.7166e-04 |    1.0310e+06    1.2344e+07 |   -4.0784e+05    1.5105e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
- (11-24 07:12) m35zams_z3 dt per 1000 cycles: c60000 t 1.122316e+05 dt 6.918412e-01 (wall 1.706984e+04 s);c61000 t 1.129352e+05 dt 8.614701e-01 (wall 1.732959e+04 s);c62000 t 1.139244e+05 dt 1.173690e+00 (wall 1.759469e+04 s);c63000 t 1.152088e+05 dt 1.130959e+00 (wall 1.786604e+04 s);c64000 t 1.165886e+05 dt 1.418600e+00 (wall 1.814086e+04 s);c65000 t 1.179651e+05 dt 8.827438e-01 (wall 1.841417e+04 s);c66000 t 1.191647e+05 dt 1.259213e+00 (wall 1.868509e+04 s);

- (11-24 07:13) status: bsg [log.2 cycle=37550 time=2.201836e+06 dt=3.342286e+01, CW 0, DT_RETRY 0]; t9l62 [log.4   , CW 0, DT_RETRY 0]; m35mms [log.4 cycle=25300 time=2.525947e+05 dt=1.200577e+01, CW 0, DT_RETRY 0]; m35zams_z3 [log.0 cycle=66760 time=1.199515e+05 dt=9.432424e-01, CW 0, DT_RETRY 0];
- (11-24 07:43) m35zams_z3 dt per 1000 cycles: c71000 t 1.273835e+05 dt 2.490599e+00 (wall 4.211671e+02 s);

- (11-24 07:47) RSG status: log.1 cycle=785350 time=7.350219e+08 dt=9.565217e+02, CW 0, DT_RETRY 0

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 08:08: 79 vet dumps
```
00068 1.224e+05  68434 |   -5.1818e-03   -3.4711e-04 |    6.9182e+05    1.2849e+07 |   -1.7294e+05    1.2642e+06
00069 1.242e+05  69488 |   -3.4919e-02   -7.0825e-05 |    8.1837e+04    9.3326e+06 |    4.9143e+05    1.5696e+06
00070  1.26e+05  70445 |   -1.2955e-02    1.6361e-04 |   -1.7635e+05    1.0591e+07 |    1.5608e+05    1.2376e+06
00071 1.278e+05  71168 |   -3.1100e-02   -1.0735e-04 |   -4.0064e+05    9.2809e+06 |    1.4999e+04    1.6218e+06
00072 1.296e+05  71891 |   -4.1281e-03   -1.7200e-04 |   -2.2608e+05    9.2019e+06 |    4.8625e+04    1.1916e+06
00073 1.314e+05  72614 |   -7.0327e-03   -7.5935e-05 |   -2.0773e+05    1.0515e+07 |   -1.6329e+05    1.2638e+06
00074 1.332e+05  73337 |   -1.4237e-02    1.5365e-04 |    9.3604e+04    9.9615e+06 |   -2.9606e+05    1.2079e+06
00075  1.35e+05  74061 |   -4.8311e-03    1.1060e-04 |   -1.5419e+04    1.0582e+07 |    5.2890e+05    1.1118e+06
00076 1.368e+05  74784 |   -3.6742e-04   -2.6356e-04 |    3.3147e+05    1.0603e+07 |    2.8535e+05    1.2104e+06
00077 1.386e+05  75508 |    2.8798e-02   -3.4859e-04 |   -1.3187e+05    9.5536e+06 |   -1.6329e+04    1.1986e+06
00078 1.404e+05  76233 |    7.6284e-03   -3.2475e-04 |    7.8709e+04    8.3006e+06 |    3.4336e+05    1.5519e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 08:10: 80 vet dumps
```
00069 1.242e+05  69488 |   -3.4919e-02   -7.0825e-05 |    8.1837e+04    9.3326e+06 |    4.9143e+05    1.5696e+06
00070  1.26e+05  70445 |   -1.2955e-02    1.6361e-04 |   -1.7635e+05    1.0591e+07 |    1.5608e+05    1.2376e+06
00071 1.278e+05  71168 |   -3.1100e-02   -1.0735e-04 |   -4.0064e+05    9.2809e+06 |    1.4999e+04    1.6218e+06
00072 1.296e+05  71891 |   -4.1281e-03   -1.7200e-04 |   -2.2608e+05    9.2019e+06 |    4.8625e+04    1.1916e+06
00073 1.314e+05  72614 |   -7.0327e-03   -7.5935e-05 |   -2.0773e+05    1.0515e+07 |   -1.6329e+05    1.2638e+06
00074 1.332e+05  73337 |   -1.4237e-02    1.5365e-04 |    9.3604e+04    9.9615e+06 |   -2.9606e+05    1.2079e+06
00075  1.35e+05  74061 |   -4.8311e-03    1.1060e-04 |   -1.5419e+04    1.0582e+07 |    5.2890e+05    1.1118e+06
00076 1.368e+05  74784 |   -3.6742e-04   -2.6356e-04 |    3.3147e+05    1.0603e+07 |    2.8535e+05    1.2104e+06
00077 1.386e+05  75508 |    2.8798e-02   -3.4859e-04 |   -1.3187e+05    9.5536e+06 |   -1.6329e+04    1.1986e+06
00078 1.404e+05  76233 |    7.6284e-03   -3.2475e-04 |    7.8709e+04    8.3006e+06 |    3.4336e+05    1.5519e+06
00079 1.422e+05  76957 |    6.9725e-02   -5.6202e-05 |   -2.3860e+05    9.4145e+06 |    1.2412e+05    1.2170e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
- (11-24 08:13) m35zams_z3 dt per 1000 cycles: c72000 t 1.298735e+05 dt 2.489743e+00 (wall 7.071221e+02 s);c73000 t 1.323616e+05 dt 2.488904e+00 (wall 9.839953e+02 s);c74000 t 1.348502e+05 dt 2.487631e+00 (wall 1.260586e+03 s);c75000 t 1.373373e+05 dt 2.487620e+00 (wall 1.540544e+03 s);c76000 t 1.398229e+05 dt 2.486076e+00 (wall 1.820501e+03 s);c77000 t 1.423089e+05 dt 2.486266e+00 (wall 2.099140e+03 s);c78000 t 1.447955e+05 dt 2.485524e+00 (wall 2.375828e+03 s);
- (11-24 08:43) m35zams_z3 dt per 1000 cycles: c79000 t 1.472806e+05 dt 2.485558e+00 (wall 2.652398e+03 s);c80000 t 1.497652e+05 dt 2.483509e+00 (wall 2.930283e+03 s);c81000 t 1.522498e+05 dt 2.483339e+00 (wall 3.218661e+03 s);c82000 t 1.547341e+05 dt 2.484810e+00 (wall 3.498421e+03 s);c83000 t 1.572186e+05 dt 2.484137e+00 (wall 3.772581e+03 s);c84000 t 1.597025e+05 dt 2.483679e+00 (wall 4.053137e+03 s);

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 09:09: 97 vet dumps
```
00087 1.566e+05  82751 |    9.6236e-02    5.2446e-05 |   -2.4079e+05    9.5885e+06 |    3.1737e+03    1.0696e+06
00088 1.584e+05  83476 |   -2.0883e-02   -2.1669e-04 |   -1.8940e+05    1.0093e+07 |   -2.2004e+05    9.8952e+05
00089 1.602e+05  84201 |    2.0214e-03    2.3688e-04 |   -3.3401e+05    7.0129e+06 |   -5.1959e+05    7.9213e+05
00090  1.62e+05  84926 |   -9.9913e-02   -3.5458e-04 |   -2.0162e+05    8.6740e+06 |    1.7133e+05    7.5303e+05
00091 1.638e+05  85650 |   -2.7081e-02   -2.9595e-04 |   -2.2548e+05    1.0173e+07 |   -1.7525e+05    7.8941e+05
00092 1.656e+05  86375 |   -6.5476e-02   -7.1287e-04 |   -2.4914e+05    1.2681e+07 |    1.2730e+05    1.0942e+06
00093 1.674e+05  87099 |    1.2481e-02   -3.2919e-04 |   -3.1271e+05    1.0357e+07 |    2.3847e+05    9.7670e+05
00094 1.692e+05  87824 |   -3.6205e-02   -2.7818e-04 |    4.0652e+05    1.1722e+07 |    8.4120e+05    1.7884e+06
00095  1.71e+05  88548 |    6.9022e-02   -4.6077e-04 |    1.6492e+04    1.0213e+07 |   -5.1045e+04    1.5201e+06
00096 1.728e+05  89271 |    2.3347e-02   -2.9118e-04 |    1.1354e+05    8.4023e+06 |   -2.4930e+05    1.1731e+06
00097 1.746e+05  89994 |    1.7748e-02    4.6021e-05 |    4.5492e+04    8.8661e+06 |   -6.6955e+05    1.3069e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```

## m35zams_z3 (Z3 layout) thin layer, signed sums (orion 24o/24zb) at 11-24 09:12: 98 vet dumps
```
00087 1.566e+05  82751 |    9.6236e-02    5.2446e-05 |   -2.4079e+05    9.5885e+06 |    3.1737e+03    1.0696e+06
00088 1.584e+05  83476 |   -2.0883e-02   -2.1669e-04 |   -1.8940e+05    1.0093e+07 |   -2.2004e+05    9.8952e+05
00089 1.602e+05  84201 |    2.0214e-03    2.3688e-04 |   -3.3401e+05    7.0129e+06 |   -5.1959e+05    7.9213e+05
00090  1.62e+05  84926 |   -9.9913e-02   -3.5458e-04 |   -2.0162e+05    8.6740e+06 |    1.7133e+05    7.5303e+05
00091 1.638e+05  85650 |   -2.7081e-02   -2.9595e-04 |   -2.2548e+05    1.0173e+07 |   -1.7525e+05    7.8941e+05
00092 1.656e+05  86375 |   -6.5476e-02   -7.1287e-04 |   -2.4914e+05    1.2681e+07 |    1.2730e+05    1.0942e+06
00093 1.674e+05  87099 |    1.2481e-02   -3.2919e-04 |   -3.1271e+05    1.0357e+07 |    2.3847e+05    9.7670e+05
00094 1.692e+05  87824 |   -3.6205e-02   -2.7818e-04 |    4.0652e+05    1.1722e+07 |    8.4120e+05    1.7884e+06
00095  1.71e+05  88548 |    6.9022e-02   -4.6077e-04 |    1.6492e+04    1.0213e+07 |   -5.1045e+04    1.5201e+06
00096 1.728e+05  89271 |    2.3347e-02   -2.9118e-04 |    1.1354e+05    8.4023e+06 |   -2.4930e+05    1.1731e+06
00097 1.746e+05  89994 |    1.7748e-02    4.6021e-05 |    4.5492e+04    8.8661e+06 |   -6.6955e+05    1.3069e+06
(float32 bins; vet_f* = cell radiation flux moments; shells i=320..331 = r 5.5279e+11..5.5603e+11; v in cm/s)
```
