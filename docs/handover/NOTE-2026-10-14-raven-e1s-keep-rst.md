# NOTE Raven -> Caltech (10-14): E1S -- please PROTECT restart files now (remap source for a 2x-lateral wedge run, USER)

USER chose to start a 2x-lateral narrow-wedge E1 (pi/6 x pi/6, 128x128, 1 block on 1 H200) by REMAPPING the running e1s state.
Raven is building the remap tool + start path; a TASK follows.
Now, please:
1. Copy (do not move) every rst still present in e1s_1014/prod/rst to e1s_1014/keep_rst/ and from now on copy each NEW rst there too
   (the link script's prune keeps only 3). Do not stop or restart e1s.
2. Reply with the list of rst times (and bin dump times / cadence and whether bins are float32 or float64) that exist.
Preferred source: the rst just BEFORE the MLT scaffold reaches w = 0 (t ~2.0-2.4e7 s) if it still exists; else the newest.
The fresh-start TASK e1w (not pushed) is replaced by this remap plan.
