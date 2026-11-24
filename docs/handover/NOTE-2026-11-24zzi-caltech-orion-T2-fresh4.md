# NOTE 2026-11-24zzi caltech -> orion: T2 = caltech-ia123rc4v-262affd6f build + gate, FRESH bsg / m35mms / t9l62 / M35ZAMS (your 24zzh)

## m35zams (11-24 11:12) on hpc-sm-01-15 / 7333036 
- start: `2026-11-24 11:10:14 chunk 0 STAR=m35zams BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/decks/IA123_m35zams.athinput md5=d2090d163`
- 300-cycle watch: FATAL 0, CW 0, created inf 0, DT_RETRY 0; never-read: ### output1/dt = 0;
- s/cyc 0.2738 (c50-320); dt c0/c100/c200/c300: 3.4805 3.4805 3.4805 3.4805 ; L_top/L_in, L_bot/L_in at c300: 1.00028 1.00003; BiCGStab it 7.60, solves 2.16 per stage line
- binder (last DTREP): c310 dt_hyd 3.481 gasonly 3.877 dir 1 r 4.7306e+11 rho 1.834e-04 T 1.911e+14 v -3.038e-01 ceff 1.988e+07 tauT_h 2.6e+04


## m35mms (11-24 11:14) on hpc-sm-02-17 / 7333073 
- start: `2026-11-24 11:10:14 chunk 0 STAR=m35mms BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/decks/IA123_m35mms.athinput md5=1c0709cfdb5`
- 300-cycle watch: FATAL 0, CW 0, created inf 0, DT_RETRY 0; never-read: ### output1/dt = 0;
- s/cyc 0.6467 (c50-300); dt c0/c100/c200/c300: 12.008 12.008 12.008 12.008 ; L_top/L_in, L_bot/L_in at c300: 1.00065 1.00011; BiCGStab it 5.46, solves 2.11 per stage line
- binder (last DTREP): c330 dt_hyd 12.01 gasonly 14.16 dir 3 r 6.7509e+11 rho 1.743e-04 T 2.212e+14 v 6.886e-03 ceff 2.265e+07 tauT_h 8.8e+04


## bsg (11-24 11:14) on hpc-sm-02-11 / 7333073 
- start: `2026-11-24 11:10:14 chunk 0 STAR=bsg BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/decks/IA123_bsg.athinput md5=875f06c8efda9256e`
- 300-cycle watch: FATAL 0, CW 0, created inf 0, DT_RETRY 0; never-read: ### output1/dt = 0;
- s/cyc 0.6867 (c50-300); dt c0/c100/c200/c300: 165.46 168.46 168.72 168.84 ; L_top/L_in, L_bot/L_in at c300: 1.06355 0.99551; BiCGStab it 2.66, solves 3.75 per stage line
- binder (last DTREP): c320 dt_hyd 168.8 gasonly 96.65 dir 3 r 1.4023e+12 rho 5.742e-05 T 1.422e+14 v 3.984e+04 ceff 1.764e+07 tauT_h 1.34e+05


## t9l62 (11-24 11:24) on hpc-sm-02-16 / 7333035 
- start: `2026-11-24 11:10:14 chunk 0 STAR=t9l62 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/decks/IA123_t9l62.athinput md5=22d549e14df6e`
- 300-cycle watch: FATAL 0, CW 0, created inf 0, DT_RETRY 0; never-read: ### output1/dt = 0;
- s/cyc 2.4333 (c50-300); dt c0/c100/c200/c300: 96.852 96.852 96.851 96.849 ; L_top/L_in, L_bot/L_in at c300: 1.00341 0.99779; BiCGStab it 4.59, solves 2.39 per stage line
- binder (last DTREP): c300 dt_hyd 96.85 gasonly 123.3 dir 3 r 1.6909e+12 rho 9.984e-05 T 1.923e+14 v 7.524e-02 ceff 2.279e+07 tauT_h 3.74e+05

- bsg (11-24 11:30): at its first periodic rst (bsg3d.00002.rst) relaunched without DTREP (hydro/dt_report=0): `2026-11-24 11:28:13 chunk 1 STAR=bsg BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingze/p`
- m35mms (11-24 11:36): at its first periodic rst (m35mms.00002.rst) relaunched without DTREP (hydro/dt_report=0): `2026-11-24 11:34:00 chunk 1 STAR=m35mms BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingz`
- t9l62 (11-24 11:36): at its first periodic rst (t9l62.00002.rst) relaunched without DTREP (hydro/dt_report=0): `2026-11-24 11:34:39 chunk 1 STAR=t9l62 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingze`
- m35zams (11-24 11:47): at its first periodic rst (m35zams.00002.rst) relaunched without DTREP (hydro/dt_report=0): `2026-11-24 11:45:41 chunk 1 STAR=m35zams BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jing`

## t9l62 rst cadence 12 h -> 24 h model time (USER via orion 24zzp) (11-24 12:51)
- At the next periodic rst: clean STOP (yes), relaunched from t9l62.00006.rst (all rsts hard-linked into rst_keep/): `2026-11-24 12:37:54 chunk 2 STAR=t9l62 BIN=/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/athena_ia123rc4v_262affd6f md5=7d7c4a71a07d3dc3878e7dea1403df23 DECK=/resnick/groups/carnegie_poc/jingze/prod_rc4v_1124/decks/IA123_t9` EX: hydro/dt_report=0 output5/dt=86400 rad_vet/clean_verbose=2 (24zzt section 0; CLEAN PASS lines in this chunk: 1348). Deck copy prod_rc4v_1124/decks/IA123_t9l62.athinput now has output5 dt = 86400 too.
- s/cyc before 2.8667 s/cyc (c1550-1800); **after (300 cycles) 2.3078 s/cyc (c1850-2150)**. CW 0, FATAL 0, DT_RETRY 0.
