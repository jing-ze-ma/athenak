# RSG fb_local floors "x2.6" (m1_pos_flr 3044 / 3031 vs 1173): where the numbers come from

Source statement: `cstan_1106/REPORT.md` line 56 ("Floors: m1_pos_flr 3044 on / 3031 off, ~2.6x pass 1 (1173)"),
summarised in `HANDOVER-2026-11-07-caltech.md` line 12. Copies of `cstan_1106/REPORT.md` and `rsgdiag_1107/REPORT.md` are in `../cstan_1106/` and `../rsgdiag_1107/`.

## Arms (all on Caltech in `$J/cstan_1106/real/<arm>/`, $J=/resnick/groups/carnegie_poc/jingze)
All arms restart from `rsg_E_1020/prodH/rst/rsg.01295.rst` (cycle0 = 589788), same binary
(md5 c52cfd88931b9e9230746524336a4a71 = `builds_rstfixx/bin/athena_gpu_he_star_m1_cts_5b03746b2_fmadoff`), 4 GPUs each, alloc 5691175.
"on"/"off" = `rad_m1/vet_gd_tangential=true/false` (the full-tensor key under test), NOT fb_local.

| pass | arm | extra keys (XK) | cycles | last cycle | m1_pos_flr | m1_vlat_fall | wall (arm.info) |
|---|---|---|---|---|---|---|---|
| 1 | rsg_on1  | cs_cross_mono=true, implicit_opac_newton_guard_mode=7 (global fallback) | 3000 | 592788 | 1173 | 2537 | 11-06 22:17 -> 11-07 00:46 |
| 1 | rsg_off1 | same | 3000 | 592788 | 1174 | 2538 | 22:17 -> 00:37 |
| 2 | rsg_on3  | + implicit_fb_local=true, implicit_fb_diag=200 | 3647 | 593435 | 3044 | 1703 | 00:53 -> stopped ~04:00 (exact PID) |
| 2 | rsg_off3 | same | 3947 | 593735 | 3031 | 1732 | 00:53 -> stopped ~04:00 |

(rsg_on2/off2 = 100-cycle aborted retries, not used.) Simulated time: pass 1 3000 cycles = 2.97e6 s; pass 2 3.7-4.0e6 s (REPORT).

## Counter
`m1_pos_flr` = column 22 of the event-counter file `rsg.log` (copied here as `<arm>.counters.log`), **cumulative since the restart**
(first row at cycle 589809: 15 / 9 / 18 / 13). The numbers are the last rows:
`tail -1 rsg.log | awk '{print $1,$20,$22}'` (cycle, m1_vlat_fall, m1_pos_flr). REPORT got them via `real/stat.sh`.

## Caveat (my note, from the same rows)
3044 / 1173 = 2.6 compares raw totals over unequal lengths (3647-3947 vs 3000 cycles). Per cycle: pass 1 0.39,
pass 2 0.83 on / 0.77 off, i.e. ~2.0-2.1x. Pass 2 differs from pass 1 by fb_local AND fb_diag=200 (fb_diag is a diagnostic,
assumed neutral; not verified). The ratio is the same in both tensor arms, so it is the fb_local path, not vet_gd_tangential.
rsgdiag_1107/REPORT.md section 2 uses the rsg_on3/off3 fbdiag files (fb_local on) for the bad-cell locations; its own arms have no fb_local
(pos_flr/cycle base 0.21).

Files here: `rsg_pass1.sh` (= real/rsg.sh), `rsg_pass2.sh` (= real/rsg2.sh), the 4 arm.info and the 4 counter files, unchanged.
