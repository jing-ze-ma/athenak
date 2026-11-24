# NOTE 2026-11-24u caltech -> orion: m35mms switched to rad_vet/floorcap_lz = true (orion 24zc, USER) with the 24r 300-cycle watch + revert rule

## m35mms (11-24 02:02): switched to floorcap_lz
- pre-switch (log.2): 1.1079 s/cyc (c6300-6550), dt 1.055138e+01. Clean STOP: yes ('Terminating on stop request'). Pre-switch rst kept as rst_keep/m35mms.00006.rst.pre-lz.
- relaunch: `2026-11-24 02:00:12 chunk 3 STAR=m35mms BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc3_af836f5a3 md5=09f50abee4ffa6cd399ce15879626353 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc3_1124/decks/tile/IA123_m35mms.athinput md5=ae635eee` + EX: time/cfl_number=0.3 hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4 rad_vet/floorcap_lz=true
- log echo (log.3): (no mg_levels line yet) | <hydro>: fofc_retest 1
- m35mms 300 cycles on floorcap_lz (11-24 02:06): CW 0, FATAL 0, created inf 0, DT_RETRY 0. **before** 1.1079 s/cyc (c6300-6550), dt 1.055138e+01 -> **after** 1.0046 s/cyc (c6700-6950), dt 5.605165e+00
