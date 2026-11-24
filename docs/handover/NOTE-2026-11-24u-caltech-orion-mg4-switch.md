# NOTE 2026-11-24u caltech -> orion: all 4 rc3 stars switched to rad_vet/clean_stag_mg_levels = 4 (24s/24t, USER-confirmed; 24r 300-cycle watch + revert rule)

## m35zams (11-24 01:32): switched to mg_levels 4
- pre-switch (log.1): 0.3654 s/cyc (c10230-10280), dt 2.478106e+00. Clean STOP: yes ('Terminating on stop request'). Pre-switch rst kept as rst_keep/m35zams.00003.rst.pre-mg4.
- relaunch: `2026-11-24 01:30:08 chunk 1 rc=0 last: time=3e+04 cycle=10408` + EX: hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4
- log echo (log.1): (no mg_levels line yet) | <hydro>: fofc_retest 1

## m35mms (11-24 01:32): switched to mg_levels 4
- pre-switch (log.1): 1.1243 s/cyc (c4950-5200), dt 9.786278e+00. Clean STOP: yes ('Terminating on stop request'). Pre-switch rst kept as rst_keep/m35mms.00004.rst.pre-mg4.
- relaunch: `2026-11-24 01:30:10 chunk 1 rc=0 last: time=6e+04 cycle=5226` + EX: hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4
- log echo (log.1): (no mg_levels line yet) | <hydro>: fofc_retest 1
- m35zams 300 cycles on mg4 (11-24 01:32): CW 0, FATAL 0, created inf 0, DT_RETRY 0. **before** 0.3654 s/cyc (c10230-10280), dt 2.478106e+00 -> **after** 0.3641 s/cyc (c10350-10400), dt 2.478581e+00
- m35mms 300 cycles on mg4 (11-24 01:32): CW 0, FATAL 0, created inf 0, DT_RETRY 0. **before** 1.1243 s/cyc (c4950-5200), dt 9.786278e+00 -> **after** 1.1243 s/cyc (c4950-5200), dt 9.786278e+00

## bsg (11-24 01:32): switched to mg_levels 4
- pre-switch (log.1): 0.6213 s/cyc (c4350-4600), dt 6.063733e+01. Clean STOP: yes ('Terminating on stop request'). Pre-switch rst kept as rst_keep/bsg3d.00006.rst.pre-mg4.
- relaunch: `2026-11-24 01:30:33 chunk 1 rc=0 last: time=5e+05 cycle=4722` + EX: hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4
- log echo (log.1): (no mg_levels line yet) | <hydro>: fofc_retest 1
- bsg 300 cycles on mg4 (11-24 01:32): CW 0, FATAL 0, created inf 0, DT_RETRY 0. **before** 0.6213 s/cyc (c4350-4600), dt 6.063733e+01 -> **after** 0.6254 s/cyc (c4450-4700), dt 5.716051e+01

- **(11-24 01:35) Correction:** the first relaunch of all 4 at 01:30 FAILED (launch_rc3f.sh was not executable: 'Permission denied'). The two 'OK ... 300 cycles' lines above for m35mms / m35zams are false: the watch read the stopped log. All 4 were down from 01:30.
  bsg (now on hpc-sm-02-16 / 7333035, freeing 7333072 for RSG) and m35mms (02-11) were relaunched at 01:34 with fofc + mg4 (`bash launch_rc3f.sh`); their 300-cycle watches are running again.
  t9l62 and m35zams stay stopped at their 01:30 rsts (t9l62.00005, hard-linked as .pre-mg4; m35zams final rst) for the 24v layout batch on 02-17 / 01-15. They restart right after it.

## bsg (11-24 01:36): switched to mg_levels 4
- pre-switch (log.1, on the old node): 0.6254 s/cyc (c4450-4700), dt 5.716051e+01. **Moved to hpc-sm-02-16 / holder 7333035 (COORD: frees 7333072 for RSG; 02-16 is timing-suspect).** Clean STOP: yes ('Terminating on stop request'). Pre-switch rst kept as rst_keep/bsg3d.00006.rst.pre-mg4.
- relaunch: `2026-11-24 01:34:48 chunk 2 STAR=bsg BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_bsg.athinput md5=9f4884524a7d23` + EX: hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4
- log echo (log.2): ### NOTE: parameter 'clean_stag_mg_levels' in block 'rad_vet' is absent from the restart's input: added from the command line (= 4); | <hydro>: fofc_retest 1

## m35mms (11-24 01:36): switched to mg_levels 4
- pre-switch (log.1, on the old node): 1.1243 s/cyc (c4950-5200), dt 9.786278e+00. **Moved to hpc-sm-02-11 / holder 7333073 (COORD: frees 7333072 for RSG; 02-16 is timing-suspect).** Clean STOP: yes ('Terminating on stop request'). Pre-switch rst kept as rst_keep/m35mms.00004.rst.pre-mg4.
- relaunch: `2026-11-24 01:34:48 chunk 2 STAR=m35mms BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_m35mms.athinput md5=ae635eee` + EX: hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4
- log echo (log.2): ### NOTE: parameter 'clean_stag_mg_levels' in block 'rad_vet' is absent from the restart's input: added from the command line (= 4); | <hydro>: fofc_retest 1
