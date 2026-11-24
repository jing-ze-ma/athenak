# NOTE 2026-11-24d caltech -> orion: release-branch GPU smoke B (IA12.2 465aff485) vs R (caltech-release-38f3d3c95), your 24b

- R: separate worktree wt_orion1121/rel38f + fresh build dir build_rel (USER/COORD: keep the rc1 / rc1-wd builds intact), CUDA 12.9.41 + g++ 13.2, HOPPER90, -ffp-contract=off, **PROBLEM=he_star_m1** (alias -> he_star_vet: `-- Including user-specified problem generator file: he_star_vet`).
  **BUILD_OK 23:35:09 md5=ad08ae9ba54268f65ac3c5ea76c0b18f TUs=258**; B = athena_ia122_465aff485 md5 cfa3dc77.
- Identical command lines with the old key names (rad_m1/...); hst only, rst/bin/log outputs off, no rst kept. bsg fresh (tile deck, nlim 30); m35mms from prod_ia122_1123/m35mms/rst_keep/m35mms.00017.rst (hard link, cycle 38385) **with the notile deck** (= the production layout, block; a restart cannot change the layout mode) + seam/wall sym, +20 cycles; t9l62 from rst_keep/t9l62.00012 + IA122_t9l62_restart.args (md5 6f8743fc), +10 cycles.
- bsg and m35mms: B and R side by side on one node (2+2 H200); t9l62: B then R on 4 H200 (2+2 does not fit at ~104 GB/GPU).

| arm | file | B md5 (non-# rows) | R md5 (non-# rows) | equal? | rc B/R | FATAL B/R | CW B/R |
|---|---|---|---|---|---|---|---|
| bsg | bsg3d.hydro.hst | 54c7c6ff3de2 | 54c7c6ff3de2 | **yes** (31 rows) | rc=0/rc=0 | 0/0 | 0/0 |
| bsg | bsg3d.user.hst | 971c7fd4c59d | 971c7fd4c59d | **yes** (31 rows) | rc=0/rc=0 | 0/0 | 0/0 |
| m35mms | m35mms.hydro.hst | bc9528cea87a | bc9528cea87a | **yes** (20 rows) | rc=0/rc=0 | 0/0 | 0/0 |
| m35mms | m35mms.user.hst | 91b682bff5eb | 91b682bff5eb | **yes** (20 rows) | rc=0/rc=0 | 0/0 | 0/0 |
| t9l62 | t9l62.hydro.hst | c6dc1d998c4e | c6dc1d998c4e | **yes** (10 rows) | rc=0/rc=0 | 0/0 | 0/0 |
| t9l62 | t9l62.user.hst | 07273eaa8f46 | 07273eaa8f46 | **yes** (10 rows) | rc=0/rc=0 | 0/0 | 0/0 |

- R log checks: "Block name 'rad_m1' ... not found": 0. Never-read keys in R but not in B:
  - bsg: ### rad_vet/clean_stag_thinf_held = false;### rad_vet/cs_cross_mono = true;### rad_vet/time2_stage_local_be = false;### rad_vet/vet_gd_rotate_align = true;### rad_vet/vet_rotate_avg_ema_change_form = linear; | 'no effect' line R: ### WARNING (scheme = clean): no effect: rad_vet/implicit_pos_floor_solve rad_vet/cs_cross_mono
  - m35mms: ### rad_vet/clean_stag_thinf_held = false;### rad_vet/cs_cross_mono = true;### rad_vet/implicit_pos_floor_solve = true;### rad_vet/time2_stage_local_be = false;### rad_vet/vet_gd_rotate_align = true;### rad_vet/vet_rotate_avg_ema_change_form = linear; | 'no effect' line R: ### WARNING (scheme = clean): no effect: rad_vet/implicit_pos_floor_solve rad_vet/cs_cross_mono
  - t9l62:  | 'no effect' line R: ### WARNING (scheme = clean): no effect: rad_vet/implicit_pos_floor_solve rad_vet/cs_cross_mono
- Files: caltech_results/RELSMOKE_1124/<arm>_{B,R}/ (hst, run.log, log.0.gz).
- Correction to the never-read list above: it compared raw names, but B prints the keys as rad_m1/... and R as rad_vet/.... After mapping rad_m1 -> rad_vet, R's extra never-read keys are only **rad_vet/cs_cross_mono and rad_vet/implicit_pos_floor_solve** (m35mms; bsg has cs_cross_mono). These are exactly the expected removed keys. So R reads every key that B reads.
