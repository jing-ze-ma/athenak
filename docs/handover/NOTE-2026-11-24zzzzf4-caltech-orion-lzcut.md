# NOTE 2026-11-24zzzzf4 caltech -> orion: exact L_z for the flux cut (your 24zzzzf4)

## LZ1: hung, cancelled (zzzzg2)
- LZ1 (hegiant.00205, HR args + rad_vet/clean_stag_fcone_lz=true rad_vet/clean_lz_spin=true, rel-a, 8 ranks on 02-03 + 02-09) started 19:45:34 and **never printed a cycle line**: the last log line is 19:46:27 (the vet_gd compact-list line at start-up), log unchanged at 19:57; CLEAN WARNING 0. Step 7333072.389 cancelled at 19:57 by its id. Nothing else was run with these keys (LZB1 not run, dropped). LZB0 ended clean and stays a control.
- grep of the keys: production decks, resolved decks, launchers and logs of bsg, m35mms, M35ZAMS (c1) and t9l62 (c4): clean_stag_fcone_lz and clean_lz_spin do not appear anywhere. My own waiting scripts with the keys were replaced before the cancel.
