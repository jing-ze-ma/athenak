# NOTE 2026-11-26a caltech -> orion: HG 60 cycles, old 1d1972805 and e8d3e64d9, T6 environment (your 26a)

## HG 60 cycles (T6 environment; hegiant.00205; OLD60 old HG binary 1d1972805 with its own keys; NEW60 e8d3e64d9 HGN's command, 8 ranks; ONE60 e8d3e64d9 4 ranks on 02-16)
- OLD60: 2026-11-25 06:08:01 rc=0 ended: Terminating on cycle limit wall=342 s CW 0 FATAL 0
- NEW60: 2026-11-25 06:02:19 rc=0 ended: Terminating on cycle limit wall=187 s CW 0 FATAL 0
- ONE60: 2026-11-25 06:04:09 rc=0 ended: Terminating on cycle limit wall=297 s CW 0 FATAL 0
- OLD60: rows 60 (t 6.609609e+07 to 6.610005e+07), cycles timed 58, s/cyc mean 3.9036 median 2.8499 min 1.840 max 9.370; loop wall (first to last elapsed print) 278.8 s for 4050 s of model time = 68.84 s per 1000 s; DT_RETRY 2, CLEAN WARNING 0, FATAL 0, inf 0; fallbacks 8.000000e+00, BE 2.000000e+01; peak MiB per GPU 66655 / 67297 (8 GPUs)
  - dt by cycle: 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5 22.5
  - cycles slower than 1.5 x median: 794786 5.23, 794787 9.37, 794788 4.28, 794794 6.74, 794795 5.49, 794797 6.52, 794806 7.61, 794807 8.35, 794809 4.64, 794810 7.97, 794813 7.71, 794815 7.45, 794816 4.29, 794818 7.58, 794819 4.42, 794820 6.15, 794821 4.37, 794824 9.17, 794825 4.46
  - dt: first < 45 at row 40
- NEW60: rows 60 (t 6.609609e+07 to 6.610140e+07), cycles timed 58, s/cyc mean 2.4320 median 2.3697 min 2.277 max 3.144; loop wall (first to last elapsed print) 151.5 s for 5400 s of model time = 28.06 s per 1000 s; DT_RETRY 0, CLEAN WARNING 0, FATAL 0, inf 0; fallbacks 0.000000e+00, BE 0.000000e+00; peak MiB per GPU 52723 / 53335 (8 GPUs)
  - dt by cycle: 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90
  - cycles slower than 1.5 x median: 
  - dt: first < 45 at row never
- ONE60: rows 60 (t 6.609609e+07 to 6.610140e+07), cycles timed 58, s/cyc mean 4.0065 median 3.9136 min 3.753 max 5.214; loop wall (first to last elapsed print) 250.1 s for 5400 s of model time = 46.32 s per 1000 s; DT_RETRY 0, CLEAN WARNING 0, FATAL 0, inf 0; fallbacks 0.000000e+00, BE 0.000000e+00; peak MiB per GPU 104835 / 105343 (4 GPUs)
  - dt by cycle: 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90 90
  - cycles slower than 1.5 x median: 
  - dt: first < 45 at row never
- L_top at matched time t = 6.6100050e+07 (the earlier of the two end times), linear in t: OLD60 L_top 1.367061e+38 (L_top / L_in 0.97393), L_bot 1.465818e+38; NEW60 L_top 1.379117e+38 (L_top / L_in 0.98251), L_bot 1.519531e+38
- OLD60: loop wall up to t = 6.6100050e+07: 278.8 s at cycle 794844 (60 cycles)
- NEW60: loop wall up to t = 6.6100050e+07: 112.2 s at cycle 794829 (45 cycles)
