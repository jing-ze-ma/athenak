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
