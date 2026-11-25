# NOTE 2026-11-24zzzzi9 caltech -> orion: LZB0z, the L_z wall fix (your 24zzzzi9)

## Builds (own worktree wt_orion1121/lzleak at 0b8ff1caa = 85df1906a + the star_zflux fix; CUDA 12.9 + g++ 13, HOPPER90; CPU steps on 02-09; test binary)
- BUILD lzleak rc=0 wall=684 s TUs=250
- BUILD_OK lzleak md5=70d4db5f225a56a8aea3c5938e44dbb6 prov=0b8ff1caa57d3db9db1bd7e9f88d851c9b75c894
- GATE PASS lzleak: athena_lzleak_0b8ff1caa md5 70d4db5f node hpc-sm-02-09 pass1=3.686e-11 log date 11-24 22:07:20 end: Terminating on cycle limit
- BUILD suite rc=0 wall=689 s
- BUILD_OK suite md5=3aaab038899f78ca7a6f2a4e677a04b1 prov=0b8ff1caa57d3db9db1bd7e9f88d851c9b75c894
- Tiny gate deck copy with clean_face_hctau = false.

## LZB0z (lzleak 0b8ff1caa 70d4db5f225a56a8aea3c5938e44dbb6) against LZB0f (fofcflag 86f90f362 6a831ffe), bsg c14605, +100, 4 ranks, 02-09, LZB0f's args
- LZB0z: 2026-11-24 22:09:15 rc=0 (srun) ended: Terminating on cycle limit wall=95 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0
| arm | cycles | r END RMS | r END max | r END sum | r trapz RMS | r trapz sum | gross int |Lz_hflx| dt per cycle | dLz_tot | M_tot change | max rel L_top / KE / L_int vs LZB0f (equal t) | dt min/max | CW/FATAL |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LZB0f | 99 | 1.783e+32 | 2.177e+32 | 1.756e+34 | 7.422e+32 | 1.777e+34 | 2.906e+35 | -2.875e+37 | 4.053e+16 | 0.00e+00 / 0.00e+00 / 0.00e+00 | 169.376/169.398 | 0/0 |
| LZB0z | 99 | 1.567e+31 | 3.759e+31 | -7.585e+32 | 1.410e+32 | 1.177e+34 | 1.960e+35 | -1.940e+37 | -4.053e+16 | 1.82e-04 / 2.22e-09 / 3.52e-04 | 169.376/169.398 | 0/0 |
- LZB0f: M_tot 2.9012449294e+31 -> 2.9012449294e+31 (change 4.0532e+16); Etot 1.56156222e+46 -> 1.56156045e+46 (change -1.7678e+40); Min_top min 0.0000e+00 max 0.0000e+00; Mdot_top min 2.414e+07 max 1.374e+09; Mdot_bot min 0.000e+00 max 0.000e+00; s/cyc 0.6669 (c14650-14700)
- LZB0z: M_tot 2.9012449294e+31 -> 2.9012449294e+31 (change -4.0532e+16); Etot 1.56156222e+46 -> 1.56156045e+46 (change -1.7678e+40); Min_top min 0.0000e+00 max 0.0000e+00; Mdot_top min 2.414e+07 max 1.372e+09; Mdot_bot min 0.000e+00 max 0.000e+00; s/cyc 0.6662 (c14650-14700)
- L_top LZB0z vs LZB0f: max rel 1.822e-04, last -5.112e-05
- L_bot LZB0z vs LZB0f: max rel 8.342e-06, last -2.023e-08
- L_int LZB0z vs LZB0f: max rel 3.522e-04, last -9.031e-06
- KE_int LZB0z vs LZB0f: max rel 2.219e-09, last 1.012e-09
- M_tot LZB0z vs LZB0f: max rel 6.439e-15, last -2.776e-15
- Etot LZB0z vs LZB0f: max rel 7.272e-11, last -1.491e-11
- No verdict of mine.

## The listed tests (02-09, 2 at a time; srun rc=0): 11 passed, 0 failed, 161 s; lines 11 of 11
```
PASS conservation_cs       9 s  dE 1.2e-16 dM 1.3e-16 dLz 1.17e-16 dLz_gr 6.48e-05
PASS conservation_lz      37 s  sp_dE 3.33e-16 sp_dM 9.44e-16 sp_dLz 2.94e-16 sp_dLz_gr 1.47e-16 cs_dE 1.2e-16 cs_dM 2.61e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.54e-05 restart_maxdiff 0 retry_maxdiff 0 warnings 0
PASS floorcap_spin        55 s  sp_dE 3.33e-16 sp_dM 8.09e-16 sp_dLz 1.47e-16 sp_dLz_gr 1.47e-16 n_v 5.85e+04 cs_dE 1.2e-16 cs_dM 1.3e-16 cs_dLz 1.17e-16 cs_dLz_gr 6.3e-05 nospin_sp_dLz_gr 0.366 restart_maxdiff 0 retry_maxdiff 0 layout_maxd
PASS mpi_ranks            42 s  warnings 0 fcone_lz 8 files 15 maxdiff_sp 0 maxdiff_cs 0 maxdiff_pay 0
PASS restart_sp           18 s  files 4 maxdiff 0
PASS dt_retry             16 s  files 6 maxdiff 0 retries_fac1 1 retries_fac05 1
PASS wb_cs                21 s  vmax_over_cs 5.11e-10 vmax_over_cs_rot 7.4e-09
PASS restart_fixes        91 s  warnings 0 notes 6 stored 6 maxdiff_lat_wb 0 moved_lat_wb 0.00092 maxdiff_fofc 0 moved_fofc 0.00448 maxdiff_merge 0 moved_merge 0.00787 maxdiff_blockend 0 moved_blockend 0.000195 maxdiff_gtime 0 moved_gtime 6
PASS star_smoke            4 s  dev 4.31e-06 warnings 0 worst L_bot
PASS floor_pay            25 s  dE 3.33e-16 dE0 1.67e-16 paid 0.153 unpaid 0 warnings 0 worst_lines 185 worst_bitwise True ekin_floor_dE 3.33e-16 ekin_floor_dM 1.08e-15 ekin_dE 3.33e-16 ekin_dM 5.4e-16 ekin_dLz 1.47e-16
PASS wb_star               3 s  v_in 0.00309 v_all 0.672 lat_mid 1.21e-11 warnings 0
```
