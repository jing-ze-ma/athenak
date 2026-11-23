# NOTE 2026-11-23f orion -> Caltech: GO - launch fresh bsg, m35mms, t9l62 on IA9 with the attached decks

USER-approved (overnight, via COORD): launch the three FRESH production runs now.

## Binary
- GitHub branch `caltech-ia9-77551e823` (unchanged since 23c). Build CUDA 12.9 + g++13; your smoke as before (post md5).

## Decks (final, validated)
- `docs/handover/caltech_ia9/IA9_{bsg,m35mms,t9l62}.athinput` on this branch (md5 in MD5SUMS).
  They are your PROD1121b J4 decks (original ICs, NOT i1_/ic_star) with the IA9 key set INSIDE <rad_m1> (per 23d), plus:
  - ray quadrature GL 7x18 staggered, ray rotation OFF (vet_rotate_avg 0, ema 0, vet_gd_rotate_dt 0, clean_rot_ring off) - USER choice;
  - clean_stag_flim / thinf / thinf_held = false; clean_stag_fs_line = pcr; clean_maxit 200; clean_cs_lz_tel false (sp decks);
  - time scheme sdirk2 + xdonor_freeze 4 + time2_k1_stiff 0 + time2_cfl_post true; faces w_fse / hctau / jse on, jse_cell off,
    bc_area on, upwind_central (TVD) on, vet_seg_plateau off.
- Before launching: replace the file paths (he_ic_file and the opacity/table files) with YOUR copies of the same files (our decks
  point at viper paths). Change nothing else.
- cfl 0.3 is in <time> (time/cfl_number = 0.3). Fresh start: give NO rad_m1 keys on the command line (they FATAL on a fresh start).
- Disk: a fresh run writes an rst at t = 0 (bsg ~22 GB, m35mms ~29 GB, t9l62 ~88 GB).

## Our smoke of exactly these decks (viper MI300A, IA9 HIP build, 55 cycles, fresh)
| run | ranks | s/cycle | dt (= your prod) | CW / fallbacks / dt_retry | L_top/L_in end |
|---|---|---|---|---|---|
| bsg | 2 | 1.41 | 165.46 | 0 / 0 / 0 | 1.0556 |
| m35mms | 4 | 0.83 | 12.008 | 0 / 0 / 0 | 1.0021 |
| t9l62 | 8 | 1.19 | 96.85 | 0 / 0 / 0 | 1.0023 |
Cycle-0 BE converged in 8-30 iterations (maxit 200); t9l62 prints some MCAP lines in the first ~17 cycles (start transient).

## Please post
Build md5 + smoke; command lines, node counts, holders; then the usual PROD table every ~20 min (CLEAN WARNING, fallbacks,
dt_retry, dt, s/cycle, L_top/L_in, L_bot/L_in). If anything FATALs or warns at startup, stop and post before retrying.

-- orion (coordinator)
