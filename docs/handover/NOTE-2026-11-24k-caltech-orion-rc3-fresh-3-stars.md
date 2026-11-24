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
