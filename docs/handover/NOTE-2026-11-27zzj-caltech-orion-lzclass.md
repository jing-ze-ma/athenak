# NOTE 2026-11-27zzj caltech -> orion: LZCLASS, the sp L_z ledger classes on three stars (your 27zzj)

## Build (worktree wt_orion1121/lzc at 920f2644d; CUDA 12.9 + g++13, star_vet; 02-11)
- BUILD lzc rc=0 wall=270 s TUs=280
- BUILD_OK lzc md5=59d2a9b8e1c2b4c910009c62ad057cef prov=920f2644de93b9fb68827213617177c8325e3b04
- GATE PASS lzc: athena_lzc_920f2644d md5 59d2a9b8 node hpc-sm-02-11 pass1=3.720e-11 log date 11-26 05:00:26 end: Terminating on cycle limit

## Arms
- LZC templates (920f2644d), fresh, 100 cycles, clean_verbose 1, output1/dcycle 1.
- LZC-bsg, LZN-bsg and LZC-m35: one GPU each, side by side on 02-11 with nothing else there.
- LZC-t9: 4 GPUs on 02-17 (t9l62 does not fit on one H200).
- Sums are over the arm, as rate x dt (the LZLEDGER values are rates: CLEAN LZB growth = sum(lzb x dt) to 4e-5 on bsg, 6e-3 on m35, 1e-6 on t9).
### LZC-bsg: 100 LZLEDGER lines; s/cycle 2.6970 (median of the last elapsed intervals); CW 0 FATAL 0 DT_RETRY 0
- sums over the arm: hflx 6.201e+34  ntel -1.206e+42  src 1.206e+42  rad -5.719e+34  cut 0  fwb 0  flr -1.087e+25  mrg -5.567e+31  gap 0  lzb 5.719e+34  rbt 5.612e+34
- rbt/lzb = 0.9813; rad + lzb = 3.836e+29; END max |.| 1.07e+27, RMS 3.29e+26 (relative to |Lz|: max 2.25e-08)
- CLEAN LZB ledger first 2.555108e+30 last 5.718926e+34 growth 5.718670e+34; sum of lzb 5.718924e+34 (ratio 0.999956)
### LZC-m35: 100 LZLEDGER lines; s/cycle 2.5365 (median of the last elapsed intervals); CW 0 FATAL 0 DT_RETRY 0
- sums over the arm: hflx 0  ntel -5.44e+33  src 5.44e+33  rad 1.076e+35  cut 0  fwb 0  flr 0  mrg 0  gap 0  lzb -1.076e+35  rbt -1.076e+35
- rbt/lzb = 1.0000; rad + lzb = 0; END max |.| 2.22e+20, RMS 7.71e+19 (relative to |Lz|: max 1.28e-14)
- CLEAN LZB ledger first 6.981671e+32 last -1.075904e+35 growth -1.082886e+35; sum of lzb -1.075904e+35 (ratio 1.006489)
### LZC-t9: 100 LZLEDGER lines; s/cycle 3.0774 (median of the last elapsed intervals); CW 0 FATAL 0 DT_RETRY 0
- sums over the arm: hflx 7.654e+23  ntel 7.659e+39  src -7.659e+39  rad 1.448e+38  cut 0  fwb 0  flr 0  mrg 0  gap 0  lzb -1.448e+38  rbt -1.448e+38
- rbt/lzb = 1.0000; rad + lzb = 1.758e+28; END max |.| 3.2e+24, RMS 1.1e+24 (relative to |Lz|: max 2.01e-14)
- CLEAN LZB ledger first 7.809699e+30 last -1.448059e+38 growth -1.448059e+38; sum of lzb -1.448058e+38 (ratio 1.000000)
### LZN-bsg: 0 LZLEDGER lines; s/cycle 2.6258 (median of the last elapsed intervals); CW 0 FATAL 0 DT_RETRY 0
- LZN-bsg vs LZC-bsg: hydro.hst byte-identical yes; old user.hst columns 1-32 equal yes (as printed text: yes); final hydro_w data identical yes

## Reading
- **rbt/lzb = 0.981 (bsg), 1.0000 (m35mms), 1.0000 (t9l62).** On every star the ledger is the radiation's L_z leaving through the physical x1 faces, as you expected.
  - On bsg the remaining 1.9 % (lzb - rbt = 1.1e33 over 100 cycles) is the part that is not boundary torque.
  - rad + lzb is 3.8e29 (bsg), 0 (m35) and 1.8e28 (t9): rad is -lzb to 1e-5 or better.
- ntel and src cancel exactly on all three stars (1.2e42 / 5.4e33 / 7.7e39).
- hflx (the gas x1 torque) is 6.2e34 on bsg, 0 on m35 and 7.7e23 on t9. mrg is -5.6e31 on bsg only.
- END: max 1.1e27 (bsg, 2e-8 of |Lz|), 2.2e20 (m35), 3.2e24 (t9).
- **LZN-bsg vs LZC-bsg:** hydro.hst byte-identical yes; the 32 old user.hst columns identical yes; final hydro_w data identical yes.
- **Cost:** LZC-bsg 2.697 vs LZN-bsg 2.626 s/cycle (+2.7 %; 1 GPU each, side by side). That includes p_ledger itself, which p_ledger_lz turns on.
- Files: docs/handover/caltech_results/LZCLASS/<arm>/ (user.hst with LzL_*, hydro.hst, lz_lines.txt, logs).
