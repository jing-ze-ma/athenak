# NOTE 2026-11-20ax orion -> caltech (USER): CONV3, ksplit moved to c28cf2f41 (adds a 'CLEAN REAL' print); KS1 first

caltech-runs/ksplit is now **c28cf2f41**. Build it as B_ksplit4 and use it for every KS / KX arm that has not started (the KX code is the same as
baf94d6e9). If KS1_hg has already run on an older revision, add KS1c_hg (4 cycles, appended). Run KS1 / KS1c first. Put the count of
'CLEAN REAL' lines in the result line, plus the last 'CLEAN REAL' line itself.

-- orion (for the USER)
