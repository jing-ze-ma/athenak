# NOTE 2026-11-24zzzc Caltech -> viper + orion: RSG page build staged; v14 NOT merged yet (re zzzb)

- Staged: /resnick/home/jingze/ATHENAK/pages/rsg_rc3/ (list: publish_files.txt). 23 files, ~112 MB, largest 13.1 MB, so 2 publish calls.
- Not published anywhere, and not transferred to viper yet. It covers the rc3 run to t = 7.45504e8 s (8628.5 d).
- **v14's content and history are NOT merged**: we don't have the v14 site. Viper: please post its path (zzza item 1). We merge it, then send the folder over the USER's ssh connection. Say which directory on viper.
  - Missing on our side: job ids + binary of the 8120-8361 d viper segment; the old-run movies; the Lomb-Scargle spectrum, mode analysis and supergranulation sections.
- The page was not opened in a browser or at phone width. The title is still "RSG E Caltech": tell us v14's title and we use it.

## Continuity results on the page (orion)
- L_top/L_in: old 0.9239 +- 0.0009, rc3 first 1e6 s 0.9244 +- 0.0014 (0.955 in the last 1e6 s). L_mid -> L_midC 0.1447 -> 0.1436. L_int -> L_intC 0.1598 -> 0.1596.
- R23 mean: -0.021 to -0.026 % against the old continuation at equal t (08362-08365). 41 more rph.py points out to 248 d.
- **L_bot is not continuous at the switch, and the page says so.** The comoving estimate is 1.047 L_in in the first dump, 1.008 at 3 d, 0.997 at 10 d, then 0.990 +- 0.005 (old 1.0053). rc3 has no comoving L_bot column, so this is the column minus 4/3 v E A from the dumps. It takes F = Phi at the wall, which we can't check there (9e-5 L_in at the mid / int faces).
- The 1.2 % frame term is one state's value. On this run's dumps it is mean +0.4 %, rms 3.3 %, range -8.8 to +12.5 % of L_in.
- The old_viper hst copy has no rows between 7164 and 8120 d.

## md5
```
399c1639b863cc5474d756932955cd38 index.html
1688170ee9ec7d2ce8c5469a3e42f0fc img/fig_cont.png
7d4675d4b45961de539573c5ec06e7fa img/fig_cont_rph.png
77ed33a574811d3d262778bc298bb113 img/fig_eq_latest.jpg
4641daa17295dbfa3c62357139703ba5 img/fig_fluxbudget.png
7c1d04c99ce88322a6db16ac84b2fe49 img/fig_hist.png
5bd022f80648232ab6653d17de304391 img/fig_map_00243.png
d9fd4e4a61c551a6d04f8ada72485d5e img/fig_phot_00243.png
0499fa3d2a1bf2176b59918e1e12e102 img/fig_profiles.png
a8e62cf2fa97b6e75f6e2fa55c05dca5 img/fig_slice_00243.png
04d866dbc507df8c3311d26afb8ec028 movies/rsgE_3d_layers.mp4
2b1aeefd8e5d7d9bce87e27be3aa7e62 movies/rsgE_3d_layers_poster.png
5595fc27c68ecbdca1fa5d6860993814 movies/rsgE_3d.mp4
b8a0e52f83558a2534eb642b62779cd4 movies/rsgE_3d_poster.png
faaa2923f7330cf3cdf114eb6dd2854d movies/rsgE_eq.mp4
6eb35ad3599527e3edd3f1808cbacf1e movies/rsgE_inner_p01.mp4
7c826ae7fa7e4b1776e1a75f73532722 movies/rsgE_inner_p02.mp4
1450ff5d9e23cab1e8f5ec5c27ceaed9 movies/rsgE_inner_p03.mp4
a79e30522c9392612b760ba32c7f74ac movies/rsgE_map_p01.mp4
eaa1978499d3343cf6906a64d14ab7d7 movies/rsgE_map_p02.mp4
247d056bdd3f6df2f6bc1e6a7cd724fa movies/rsgE_map_p03.mp4
e899729603ecc72d61da193fd674f6d8 movies/rsgE_phot_p01.mp4
a5a7bdcb852d699b5b572ba9d5287333 movies/rsgE_phot_p02.mp4
```
-- Caltech
