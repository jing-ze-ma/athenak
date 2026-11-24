# NOTE 2026-11-24l caltech -> orion: t9l62 WD0 / WD2 (24f) and m35mms PC arms P0 / PR / PRL4 / PL4 (24g)

- WD binary: BUILD_OK 00:25:06 md5=701c4c1cbe1068fbcb4b02161378723d TUs=220 wall=267 s (incremental in the rel38f worktree's build_rel, PROBLEM he_star_m1 -> he_star_vet)
  GATE PASS (seam+wall sym deck copy): athena_t9floorwd2_140c271e0 md5 701c4c1c node hpc-sm-02-16 pass1=3.686e-11
- WD: rst_keep/t9l62.00012 + WD_t9l62_restart.args + nlim 6480 (+ WD2: hydro/fofc=true hydro/fofc_retest=true), hpc-sm-02-16 (7333035), 4 ranks.
- PC: athena_ia123rc3_af836f5a3, from prod_ia122_1123/m35mms/rst_keep/m35mms.00017 (cycle 38385), -i rc3 notile IA123_m35mms (= production layout) + K keys, 1000 cycles (nlim 39385), 4 H200 alone on 02-17 (and 02-03 once free; not 02-16: suspect timing; 01-15 now runs M35ZAMS).
- Script tests_1124/run_X.sh; files caltech_results/TESTS_1124/<arm>/.

## WD0 (11-24 00:29)
- 2026-11-24 00:29:58 rc=0 last: 
- FATAL 0, CW 2, stage fallbacks 1, DT_RETRY 0, created inf 0; dt min/last 2.485 / 4.962
- echo: 'floor_from_ekin' in block 'hydro' is absent from the restart's input: added from the comm;'floor_audit_worst' in block 'hydro' is absent from the restart's input: added from the co; <hydro> dual_energy 0 (eta 0.001, slots -1,-1) floor_pay rad_column fofc_retest 0 floor_audit 1
- [update] per cycle (last line of each cycle; floor_de, skip_de cumulative, and their per-cycle differences):
```
 cycle     floor_de   d_floor_de      skip_de    d_skip_de
  6457   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00
  6458   3.2725e+36   3.2725e+36   2.5818e+22   2.5818e+22
  6459   8.1063e+38   8.0736e+38   7.7474e+38   7.7474e+38
  6460   9.9508e+38   1.8444e+38   1.4441e+40   1.3666e+40
  6461   6.9123e+40   6.8128e+40   2.2819e+41   2.1375e+41
  6462   6.9159e+40   3.6143e+37   2.2819e+41   0.0000e+00
  6463   1.7407e+41   1.0491e+41   1.3782e+42   1.1500e+42
  6464   1.7407e+41   3.3928e+36   1.9083e+42   5.3014e+41
  6465   1.7784e+41   3.7751e+39   1.9293e+42   2.1045e+40
  6466   1.9509e+41   1.7248e+40   1.9293e+42   0.0000e+00
  6467   1.9591e+41   8.1354e+38   2.3019e+42   3.7260e+41
  6468   1.9608e+41   1.7793e+38   2.3055e+42   3.5337e+39
  6469   1.9621e+41   1.2166e+38   2.3517e+42   4.6266e+40
  6470   1.9715e+41   9.4589e+38   2.3693e+42   1.7510e+40
  6471   1.9727e+41   1.1434e+38   2.3710e+42   1.7264e+39
  6472   2.0023e+41   2.9672e+39   2.3710e+42   0.0000e+00
  6473   2.0173e+41   1.4947e+39   2.3710e+42   0.0000e+00
  6474   2.0470e+41   2.9705e+39   2.3710e+42   0.0000e+00
  6475   5.5294e+41   3.4824e+41   2.3759e+42   4.9632e+39
  6476   5.5322e+41   2.7566e+38   2.8434e+42   4.6748e+41
  6477   5.5392e+41   7.0888e+38   2.8710e+42   2.7599e+40
  6478   5.5544e+41   1.5188e+39   2.8710e+42   0.0000e+00
  6479   5.6301e+41   7.5702e+39   2.8710e+42   0.0000e+00
```
- floor_audit_worst lines in the blob (rank 0 gid 0, r 2.29-2.45e12): 55 (all in grep.txt)
- FOFC: FOFC first-order cell-stages=0; last pay line: fofc_retest: retested 0 unflagged 0

## WD2 (11-24 00:32)
- 2026-11-24 00:32:23 rc=0 last: 
- FATAL 0, CW 0, stage fallbacks 0, DT_RETRY 0, created inf 0; dt min/last 3.524 / 6.556
- echo: 'floor_from_ekin' in block 'hydro' is absent from the restart's input: added from the comm;'floor_audit_worst' in block 'hydro' is absent from the restart's input: added from the co;'fofc_retest' in block 'hydro' is absent from the restart's input: added from the command ; <hydro> dual_energy 0 (eta 0.001, slots -1,-1) floor_pay rad_column fofc_retest 1 floor_audit 1
- [update] per cycle (last line of each cycle; floor_de, skip_de cumulative, and their per-cycle differences):
```
 cycle     floor_de   d_floor_de      skip_de    d_skip_de
  6457   0.0000e+00   0.0000e+00   0.0000e+00   0.0000e+00
  6458   0.0000e+00   0.0000e+00   3.0149e+22   3.0149e+22
  6459   0.0000e+00   0.0000e+00   3.0925e+22   7.7630e+20
  6460   0.0000e+00   0.0000e+00   3.8644e+22   7.7184e+21
  6461   0.0000e+00   0.0000e+00   6.4845e+22   2.6201e+22
  6462   0.0000e+00   0.0000e+00   7.0585e+22   5.7397e+21
  6463   0.0000e+00   0.0000e+00   1.1349e+23   4.2904e+22
  6464   0.0000e+00   0.0000e+00   1.1655e+23   3.0627e+21
  6465   0.0000e+00   0.0000e+00   1.2857e+23   1.2022e+22
  6466   0.0000e+00   0.0000e+00   1.2857e+23   2.7000e+15
  6467   0.0000e+00   0.0000e+00   1.2905e+23   4.7337e+20
  6468   0.0000e+00   0.0000e+00   1.7071e+23   4.1663e+22
  6469   0.0000e+00   0.0000e+00   1.7672e+23   6.0112e+21
  6470   0.0000e+00   0.0000e+00   1.7696e+23   2.3774e+20
  6471   0.0000e+00   0.0000e+00   1.7722e+23   2.5580e+20
  6472   0.0000e+00   0.0000e+00   5.0295e+23   3.2574e+23
  6473   0.0000e+00   0.0000e+00   5.0349e+23   5.3792e+20
  6474   0.0000e+00   0.0000e+00   5.0387e+23   3.7949e+20
  6475   0.0000e+00   0.0000e+00   5.0451e+23   6.3698e+20
  6476   0.0000e+00   0.0000e+00   5.0475e+23   2.4456e+20
  6477   0.0000e+00   0.0000e+00   5.0558e+23   8.2954e+20
  6478   0.0000e+00   0.0000e+00   5.1113e+23   5.5484e+21
  6479   0.0000e+00   0.0000e+00   5.1114e+23   1.2172e+19
```
- floor_audit_worst lines in the blob (rank 0 gid 0, r 2.29-2.45e12): 21 (all in grep.txt)
- FOFC: FOFC first-order cell-stages=0; last pay line: fofc_retest: retested 178507484 unflagged 0
