# NOTE 2026-11-24zzzzj6 caltech -> orion: luminosity continuity of the rel-f restart tests (your 24zzzzj6)

## L continuity and smoothness, rel-f (b09a7b3c1) restart, five stars, from existing histories (read-only; R50 = rel-f, C50 = running binary, both 50 steps from the same rst; production rows before and after; (11-24 23:06)
- Columns: the value at the rst (last production row before it, also as L/L_in); the star's own noise = mean and rms of the per-step relative change over the 50 production rows before the rst; the jump is the first restarted row against the last production row before the rst. (jump - mean)/rms is in units of that noise. L_bot is in the table as its own row.


### m35zams: rst m35zams_m35zams.00003.rst, cycle 50816, t 1.768656e+05; last production row before it t 1.768656e+05 (dt 3.4805); production rows after it: 28488; R50 first row t 1.768690e+05, C50 1.768690e+05
| L | value at the rst (L/L_in) | own step-to-step noise: mean, rms dev (rel per step, 50 rows before) | jump at the 1st step: rel-f | running binary | production after | (jump - mean) / rms: rel-f, running | steps 1-10 max |rel step change|: rel-f, running | row 50 vs the rst: rel-f, running | mean rel change per step rows 11-50: rel-f, running |
|---|---|---|---|---|---|---|---|---|---|
| L_top | 2.1601e+35 (1.00014) | -6.32e-07, 4.12e-07 | -5.61e-07 | -5.61e-07 | -5.61e-07 | 0.2, 0.2 | 7.51e-07, 7.52e-07 | -2.13e-05, -2.22e-05 | -3.94e-07, -4.17e-07 |
| L_bot | 2.1599e+35 (1.00001) | -3.77e-09, 1.57e-08 | -1.09e-08 | -1.10e-08 | -1.10e-08 | -0.5, -0.5 | 2.31e-08, 2.66e-08 | -2.12e-07, -2.10e-07 | -4.85e-09, -4.88e-09 |
| L_mid | 2.1586e+35 (0.99941) | -1.44e-06, 2.62e-07 | -9.66e-07 | -9.66e-07 | -9.66e-07 | 1.8, 1.8 | 9.66e-07, 9.66e-07 | -2.48e-05, -2.48e-05 | -4.00e-07, -4.00e-07 |
| L_int | 2.1601e+35 (1.00014) | -6.73e-07, 2.76e-07 | -5.90e-07 | -5.90e-07 | -5.90e-07 | 0.3, 0.3 | 1.12e-06, 1.12e-06 | -2.40e-05, -2.37e-05 | -4.43e-07, -4.35e-07 |
- R50 vs C50 (same rst): max rel |L_top| 9.07e-07, |L_bot| 2.90e-08, |L_mid| 1.24e-09, |L_int| 1.11e-06 over the 50 rows; at row 0: 3.21e-12, 1.23e-10, -1.19e-11, 1.05e-12


### m35mms: rst m35mms_m35mms.00004.rst, cycle 22144, t 2.592005e+05; last production row before it t 2.592005e+05 (dt 12.0067); production rows after it: 3356; R50 first row t 2.592125e+05, C50 2.592125e+05
| L | value at the rst (L/L_in) | own step-to-step noise: mean, rms dev (rel per step, 50 rows before) | jump at the 1st step: rel-f | running binary | production after | (jump - mean) / rms: rel-f, running | steps 1-10 max |rel step change|: rel-f, running | row 50 vs the rst: rel-f, running | mean rel change per step rows 11-50: rel-f, running |
|---|---|---|---|---|---|---|---|---|---|
| L_top | 1.0390e+37 (1.01313) | 2.82e-05, 7.69e-05 | 3.47e-05 | 3.47e-05 | 3.47e-05 | 0.1, 0.1 | 1.73e-04, 1.73e-04 | 1.29e-03, 1.30e-03 | 1.74e-05, 1.74e-05 |
| L_bot | 1.0256e+37 (1.00001) | -2.96e-07, 1.68e-08 | -2.59e-07 | -2.56e-07 | -2.56e-07 | 2.2, 2.4 | 2.68e-07, 2.75e-07 | -1.06e-05, -1.06e-05 | -2.02e-07, -2.02e-07 |
| L_mid | 1.0023e+37 (0.97729) | 1.26e-04, 4.63e-06 | 1.32e-04 | 1.32e-04 | 1.32e-04 | 1.3, 1.3 | 1.33e-04, 1.33e-04 | 6.56e-03, 6.56e-03 | 1.30e-04, 1.30e-04 |
| L_int | 1.0669e+37 (1.04035) | 7.05e-05, 3.23e-04 | 3.19e-04 | 3.19e-04 | 3.19e-04 | 0.8, 0.8 | 9.12e-04, 9.12e-04 | -2.21e-03, -2.17e-03 | -4.77e-05, -4.65e-05 |
- R50 vs C50 (same rst): max rel |L_top| 4.08e-05, |L_bot| 1.45e-08, |L_mid| 5.26e-10, |L_int| 1.36e-04 over the 50 rows; at row 0: -1.68e-14, -2.54e-09, 1.14e-12, -7.69e-14


### bsg: rst bsg_bsg3d.00008.rst, cycle 29127, t 4.838566e+06; last production row before it t 4.838566e+06 (dt 169.4271); production rows after it: 530; R50 first row t 4.838736e+06, C50 4.838736e+06
| L | value at the rst (L/L_in) | own step-to-step noise: mean, rms dev (rel per step, 50 rows before) | jump at the 1st step: rel-f | running binary | production after | (jump - mean) / rms: rel-f, running | steps 1-10 max |rel step change|: rel-f, running | row 50 vs the rst: rel-f, running | mean rel change per step rows 11-50: rel-f, running |
|---|---|---|---|---|---|---|---|---|---|
| L_top | 4.7424e+37 (1.01961) | -2.40e-05, 2.01e-05 | -2.45e-05 | -2.45e-05 | -2.45e-05 | -0.0, -0.0 | 2.45e-05, 2.45e-05 | -8.36e-04, -8.18e-04 | -1.64e-05, -1.59e-05 |
| L_bot | 4.6779e+37 (1.00574) | 1.13e-05, 6.27e-06 | -9.92e-07 | -4.03e-07 | -4.03e-07 | -2.0, -1.9 | 4.38e-06, 4.47e-06 | -4.38e-04, -4.37e-04 | -1.03e-05, -1.03e-05 |
| L_mid | 4.7600e+37 (1.02339) | 6.55e-06, 2.97e-06 | 8.97e-06 | 8.97e-06 | 8.97e-06 | 0.8, 0.8 | 9.03e-06, 9.03e-06 | 3.44e-04, 3.44e-04 | 6.36e-06, 6.36e-06 |
| L_int | 4.7251e+37 (1.01589) | 4.32e-05, 4.48e-05 | 7.76e-06 | 7.76e-06 | 7.76e-06 | -0.8, -0.8 | 5.29e-05, 5.29e-05 | -6.79e-04, -6.83e-04 | -2.47e-05, -2.48e-05 |
- R50 vs C50 (same rst): max rel |L_top| 3.03e-05, |L_bot| 7.26e-06, |L_mid| 2.99e-06, |L_int| 4.14e-05 over the 50 rows; at row 0: 9.40e-12, -5.89e-07, -4.97e-13, 3.35e-14


### t9l62: rst t9l62_t9l62.00003.rst, cycle 6332, t 5.184146e+05; last production row before it t 5.184146e+05 (dt 54.3925); production rows after it: 939; R50 first row t 5.184690e+05, C50 5.184690e+05
| L | value at the rst (L/L_in) | own step-to-step noise: mean, rms dev (rel per step, 50 rows before) | jump at the 1st step: rel-f | running binary | production after | (jump - mean) / rms: rel-f, running | steps 1-10 max |rel step change|: rel-f, running | row 50 vs the rst: rel-f, running | mean rel change per step rows 11-50: rel-f, running |
|---|---|---|---|---|---|---|---|---|---|
| L_top | 1.8772e+39 (2.14197) | -1.46e-04, 5.46e-06 | -1.51e-04 | -1.51e-04 | -1.51e-04 | -1.0, -1.0 | 1.54e-04, 1.54e-04 | -6.86e-03, -6.86e-03 | -1.36e-04, -1.36e-04 |
| L_bot | 8.7882e+38 (1.00280) | 3.74e-04, 3.06e-05 | 4.13e-04 | 4.13e-04 | 4.13e-04 | 1.3, 1.3 | 4.13e-04, 4.13e-04 | 1.82e-02, 1.82e-02 | 3.47e-04, 3.47e-04 |
| L_mid | 1.8358e+39 (2.09483) | 4.44e-05, 4.48e-04 | -8.89e-05 | -8.89e-05 | -8.89e-05 | -0.3, -0.3 | 7.03e-04, 7.03e-04 | 1.21e-03, 1.57e-03 | 4.63e-05, 5.55e-05 |
| L_int | 2.9649e+39 (3.38321) | -9.82e-05, 1.16e-06 | -9.75e-05 | -9.75e-05 | -9.75e-05 | 0.6, 0.6 | 9.75e-05, 9.75e-05 | -4.66e-03, -4.66e-03 | -9.26e-05, -9.26e-05 |
- R50 vs C50 (same rst): max rel |L_top| 8.56e-07, |L_bot| 4.29e-08, |L_mid| 3.66e-04, |L_int| 1.74e-08 over the 50 rows; at row 0: 2.22e-16, -1.03e-10, -3.00e-15, 3.23e-13


### t9mlt: rst t9mlt_t9l62.00001.rst, cycle 2683, t 2.592432e+05; last production row before it t 2.592432e+05 (dt 96.6222); production rows after it: 532; R50 first row t 2.593398e+05, C50 2.593398e+05
| L | value at the rst (L/L_in) | own step-to-step noise: mean, rms dev (rel per step, 50 rows before) | jump at the 1st step: rel-f | running binary | production after | (jump - mean) / rms: rel-f, running | steps 1-10 max |rel step change|: rel-f, running | row 50 vs the rst: rel-f, running | mean rel change per step rows 11-50: rel-f, running |
|---|---|---|---|---|---|---|---|---|---|
| L_top | 1.0206e+39 (1.16459) | 4.28e-05, 2.51e-05 | 1.35e-05 | 1.35e-05 | 1.35e-05 | -1.2, -1.2 | 1.44e-05, 1.44e-05 | -9.96e-04, -9.95e-04 | -2.64e-05, -2.64e-05 |
| L_bot | 8.7675e+38 (1.00043) | -1.00e-05, 1.04e-06 | -8.34e-06 | -8.33e-06 | -8.33e-06 | 1.6, 1.6 | 8.55e-06, 8.55e-06 | -2.37e-04, -2.37e-04 | -3.96e-06, -3.96e-06 |
| L_mid | 1.0050e+39 (1.14679) | 1.59e-04, 8.88e-05 | 4.53e-05 | 4.53e-05 | 4.53e-05 | -1.3, -1.3 | 1.24e-04, 1.24e-04 | 3.26e-04, 3.26e-04 | -5.73e-06, -5.74e-06 |
| L_int | 1.0747e+39 (1.22630) | 8.75e-05, 5.51e-04 | 1.10e-03 | 1.10e-03 | 1.10e-03 | 1.8, 1.8 | 1.10e-03, 1.10e-03 | 2.19e-02, 2.19e-02 | 4.48e-04, 4.47e-04 |
- R50 vs C50 (same rst): max rel |L_top| 5.12e-07, |L_bot| 6.30e-08, |L_mid| 2.35e-05, |L_int| 6.09e-05 over the 50 rows; at row 0: -1.11e-16, -1.03e-08, -3.33e-16, -3.04e-13

### Reading of the numbers (no interpretation beyond them)
- **First step:** in all five stars and all four luminosities the rel-f jump is within -2.0 and +2.4 of the star's own step rms (largest: m35mms L_bot +2.2; bsg L_bot -2.0), and equals the running binary's and the production run's own next step to 2 or 3 digits or better. **No jump like the +4.6 % L_bot jump of RSG.** The largest first-step change in L_bot is 4.1e-4 (t9l62, whose own steps are 3.7e-4 at that point, outburst state); elsewhere 1e-6 to 8e-6 (bsg -9.9e-7, its running binary -4.0e-7).
- **Steps 1-10, and the slope of rows 11-50:** equal between rel-f and the running binary to 2-3 digits (table).
- **R50 vs C50 over 50 steps:** max |L_top| 9.1e-7 (m35zams), 4.1e-5 (m35mms), 3.0e-5 (bsg), 8.6e-7 (t9l62), 5.1e-7 (t9mlt); |L_bot| 2.9e-8, 1.5e-8, 7.3e-6, 4.3e-8, 6.3e-8.
- **Radius and T_eff: not in this table.** No dump exists at cycle 0, 1 or 10 for any star; only the end dump (rst + 50 steps) of R50 and C50. The page agent reduces those (directories below).
- t9mlt on restart prints the ramp 'from restart file (computed at the fresh start)', no 'no stored numbers'.

### End dumps (kept, nothing deleted)
- m35zams R50: m35zams.hydro_w.00099.bin m35zams.vet.00099.bin in b1_1124/rc/m35zams/R50/bin
- m35zams C50: m35zams.hydro_w.00099.bin m35zams.vet.00099.bin in b1_1124/rc/m35zams/C50/bin
- m35mms R50: m35mms.hydro_w.00145.bin m35mms.vet.00145.bin in b1_1124/rc/m35mms/R50/bin
- m35mms C50: m35mms.hydro_w.00145.bin m35mms.vet.00145.bin in b1_1124/rc/m35mms/C50/bin
- bsg R50: bsg3d.hydro_w.00113.bin bsg3d.vet.00113.bin in b1_1124/rc/bsg/R50/bin
- bsg C50: bsg3d.hydro_w.00113.bin bsg3d.vet.00113.bin in b1_1124/rc/bsg/C50/bin
- t9l62 R50: t9l62.hydro_w.00049.bin t9l62.vet.00049.bin in b1_1124/rc/t9l62/R50/bin
- t9l62 C50: t9l62.hydro_w.00049.bin t9l62.vet.00049.bin in b1_1124/rc/t9l62/C50/bin
- t9mlt R50: t9l62.hydro_w.00025.bin t9l62.vet.00025.bin in b1_1124/rc/t9mlt/R50/bin
- t9mlt C50: t9l62.hydro_w.00025.bin t9l62.vet.00025.bin in b1_1124/rc/t9mlt/C50/bin

## Branch clean-up (zzzzj6 item 3)
- The t9mlt arm and the continuity table are done; none of my tests still needs a fetch from a caltech-* branch on GitHub. My local worktrees keep their objects. I deleted nothing. Branches: 144; with a local worktree at their tip: 97 (they would not be lost by deleting the branch).
```
caltech-c2-d4ae795a4  wt/c2
caltech-c3-ce559f800  wt/c3
caltech-c4-a000d803c  wt/c4
caltech-c4speed-d979d5268  wt/c4speed
caltech-fofcflag-86f90f362  wt/fofcflag
caltech-idortchk-29046c3f2  wt/wt_cpu
caltech-latwb2-c01bc242f  wt/latwb2
caltech-latwb-83663792f  wt/latwb
caltech-lzleak-0b8ff1caa  wt/lzleak
caltech-rel-a-c880e2703  wt/rela
caltech-rel-b-c8a0fc099  wt/relb
caltech-rel-e-03e874599  wt/rele2
caltech-rel-e-08c48ac10  wt/rele
caltech-rel-f-b09a7b3c1  wt/relf
caltech-runs/beamfix12  wt/bfx12
caltech-runs/beamfix13  wt/bfx13
caltech-runs/beamfix16  wt/bfx16
caltech-runs/beamfix21  wt/bfx21
caltech-runs/beamfix7  wt/bfx7
caltech-runs/beamfixc20  wt/bfxc20
caltech-runs/beamfixc  wt/bfxc
caltech-runs/beamfixf  wt/bfxf
caltech-runs/beamfixg  wt/bfxg
caltech-runs/beamfix  wt/bfx
caltech-runs/conserv-ctl2  wt/conserv-ctl2
caltech-runs/conserv-ctl  wt/conserv-ctl
caltech-runs/conserv-tie  wt/conserv-tie
caltech-runs/convdiag  wt/convdiag
caltech-runs/conv  wt/conv
caltech-runs/eleak  wt/eleak2
caltech-runs/emawb  wt/emawb
caltech-runs/geom-ca  wt/geom-ca
caltech-runs/geom-fm2  wt/geom-fm2
caltech-runs/h200-p4  wt/h200-p4
caltech-runs/h200spd  wt/h2c wt/h2k wt/h2m
caltech-runs/h200  wt/h200 wt/h200mb2
caltech-runs/icbal  wt/icbal wt/icbaln
caltech-runs/kern-fw  wt/kfw
caltech-runs/ksplitdiag  wt/ksplitdiag
caltech-runs/ksplit  wt/ksplit7
caltech-runs/lag20  wt/lag20
caltech-runs/lean2  wt/lean2
caltech-runs/lean  wt/lean
caltech-runs/lzbook  wt/lzbook
caltech-runs/lzcs2  wt/lzcs2
caltech-runs/lzcs-tel  wt/lzcstel
caltech-runs/lzcs  wt/lzcs
caltech-runs/lztel2-study  wt/lztel2s
caltech-runs/lztel2  wt/lztel2
caltech-runs/lztel3-fixstudy  wt/lztel3f
caltech-runs/lztel3-study  wt/lztel3s
caltech-runs/lztel3  wt/lztel3
caltech-runs/lztel  wt/lztel
caltech-runs/pcsym  wt/pcsym
caltech-runs/realf  wt/realf
caltech-runs/realz2  wt/rz2
caltech-runs/realz3  wt/rz3b
caltech-runs/realz  wt/realz2
caltech-runs/spd2b  wt/spd2b
caltech-runs/speedhead  wt/speedhead
caltech-runs/thinf  wt/thinf
caltech-runs/tie-check  wt/tie-check
caltech-runs/tie-prof  wt/tie-prof
caltech-runs/v201b  wt/v201b
caltech-runs/v201d  wt/v201d
caltech-runs/v201e  wt/v201e
caltech-runs/v201  wt/v201
caltech-runs/v202  wt/v202
caltech-runs/v204c  wt/v204c
caltech-runs/v204  wt/v204
caltech-runs/v20b  wt/v20b
caltech-runs/v20f  wt/v20f
caltech-runs/v20g2  wt/v20g2
caltech-runs/v20g3  wt/v20g3
caltech-runs/v20g  wt/v20g
caltech-runs/v20j  wt/v20j
caltech-runs/v20l  wt/v20l
caltech-runs/v20m  wt/v20m
caltech-runs/v20nk  wt/v20nk
caltech-runs/v20n  wt/v20n
caltech-runs/xcap  wt/xcap
caltech-runs/xpos2  wt/xpos2
caltech-runs/xpos  wt/xpos2
caltech-runs/xtloc  wt/xtloc
caltech-sdocell-883d79d97  wt/sdocell
caltech-sdomom-373e2cc1c  wt/sdomom
caltech-sdorder2-09d1e4df4  wt/sdo2
caltech-sdorder3-95813cf85  wt/sdo3
caltech-sdorder5-10bf5513c  wt/sdorder5
caltech-sdorder-a840d7e92  wt/rel38f
caltech-sdosdo-53290b5af  wt/sdosdo
caltech-speed24-e7b2eea39  wt/speed24
caltech-starname-4a26282a6  wt/star
caltech-t9hot-compton-f8e55fee7  wt/wt_cpu_c
caltech-test-fwbx  wt/v204n
caltech-thinfcs-87af1d65c  wt/thinfcs
caltech-zamsfast-3b44a361d  wt/v203
```
- Not fetched by me since: 47 further caltech-* branches have no local worktree of mine: caltech-dflfix-5b6159ad0 caltech-halfday-83b71abed caltech-ia11-bd9553537 caltech-ia121-21e176c67 caltech-ia121-fa59bffd6 caltech-ia122-465aff485 caltech-ia123rc1-da92efe70 caltech-ia123rc1-wd-8197e0c4f caltech-ia123rc2-5b493acfa caltech-ia123rc3-af836f5a3 caltech-ia123rc3-fcapspin-6f2973dc8 caltech-ia123rc4v-262affd6f caltech-ia9-77551e823 caltech-m35fast-13679216a caltech-prod-462210d51 caltech-prod-bf4f917d4 caltech-rel-d-78cf8880b caltech-release-38f3d3c95 caltech-rel-e-f62c8879e caltech-runs/beamreal caltech-runs/c3c caltech-runs/c3h caltech-runs/joint3 caltech-runs/joint4d caltech-runs/rsgdie caltech-runs/stall2 caltech-runs/stall2v3 caltech-runs/v203 caltech-runs/v204n caltech-runs/v2
