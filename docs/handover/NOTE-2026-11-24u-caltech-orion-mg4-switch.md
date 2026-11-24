# NOTE 2026-11-24u caltech -> orion: all 4 rc3 stars switched to rad_vet/clean_stag_mg_levels = 4 (24s/24t, USER-confirmed; 24r 300-cycle watch + revert rule)

## m35zams (11-24 01:32): switched to mg_levels 4
- pre-switch (log.1): 0.3654 s/cyc (c10230-10280), dt 2.478106e+00. Clean STOP: yes ('Terminating on stop request'). Pre-switch rst kept as rst_keep/m35zams.00003.rst.pre-mg4.
- relaunch: `2026-11-24 01:30:08 chunk 1 rc=0 last: time=3e+04 cycle=10408` + EX: hydro/fofc=true hydro/fofc_retest=true rad_vet/clean_stag_mg_levels=4
- log echo (log.1): (no mg_levels line yet) | <hydro>: fofc_retest 1
