# TASK 11-13 viper->Caltech (USER "give something for Caltech to run, maybe dhj mhd"): two dhj MHD runs

This supersedes UPDATE2 (cf8116fda): **WASP-189b 30 G now runs on Caltech, not on Raven.**
Its Raven queue task is withdrawn and the Raven run dir has CANCEL, so nothing runs it twice.

Bundle on viper: /viper/ptmp2/jinma/caltech_xfer_1113/ (the user copies it). md5s are in md5_bundle_1113.txt.
The decks and run.cfg carry /raven paths; repoint the pkg/ck/ic/start/remap paths to your copy (same md5).
Binary: our dhj build of prod-1118 1d1972805, or dhjfast c3ee3a77 (the Raven runs use it; same dhj physics).
Run the 30 G recipe as for toi1518b_mhd30_hyd_1112, with the KEYS from run.cfg (incl. mhd/valfven_max_cgs=5e6).

## 1. WASP-189b 30 G (priority 1): w189g30/
- The 30 G boot is already done from our hydro rst dhj.00116 (rot 109.000, md5 433731ba).
- Restart from w189g30/dhj.00118.rst with run.cfg KEYS; TLIM 28478085.8 s (= rot 109 + 12).
- bbot_gauss 35.167587 (30 G polar at R_p 1.619 R_J); the boot measured max |B| 35.14 G, div B 3.8e-15.
- At TLIM run fe30.sh on the last 4 bins (fe3d pc6 + the ESPRESSO/MAROON-X w189 analysis). Observed: 1.05 +- 1.59 km/s; hydro gives 9.02.
- (w189g30/start/remap.dat + run.athinput can redo the boot if you prefer: time/nlim=30 from the deck, then -r.)

## 2. TOI-1518b 20 G (priority 2): toi20g/
- Same as toi1518b_mhd30_hyd_1112 except bbot_gauss = 24.377582 (20 G at R_p).
- Fresh boot from toi20g/start/remap.dat (hydro rot ~190; start_time 31397735.019 s): boot -i run.athinput time/nlim=30, then -r.
- TLIM 33370362.2 s (as the 30 G run); Fe I via fe30.sh on the last 4 bins. Also report KE and ME vs time: the 30 G run was still decelerating.
- Why: Fe I terminator shift vs field: 1 G -7.74, 10 G -6.78, 30 G -3.00 km/s; observed -3.82 +- 0.50. 20 G brackets where the field matches.

Reply by NOTE (Caltech->viper) with smoke status, s/cycle, and Fe I results when done.
