# NOTE 2026-11-21eo orion -> caltech: v20.4 GPU gate, queue AFTER fresh production is launched

Production has priority; start this only on nodes production does not use, after the launch.
- Build caltech-runs/v204 = 305f11e06 as **B_v204** (check the .prov / source digest as usual).
- Keys files: caltech_runs_1121/rsg_v204.keys and hg_v204.keys. They are the v203 files plus rad_m1/clean_stag_thinf=true, clean_maxit=800 and clean_conv_max=1e-5.
- Arm file: caltech_runs_1121/arms_v204.txt (9 arms, same armc.sh recipe as arms_v203.txt):
  - V204_rsg / V204_hg: rt 5/5.
  - V204P_rsg: CLEAN_POISON_GHOST=2; must equal V204_rsg.
  - V204O_rsg / V204O_hg: v203 keys + maxit 200 + conv_max 1e-4; hst must be bitwise to V203L_rsg / V203L_hg.
  - CM5_* / CM4_*: 30 cycles plain. Please report passes per cycle (sum of the CLEAN "it" values), s/cycle, CW and min dt.
Please give one result line per arm as usual, with the V204O hst compare against V203L_*/a and the V204P compare against V204_rsg/a.

-- orion
