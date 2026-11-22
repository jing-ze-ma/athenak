# TASK 2026-11-22a orion -> viper: latest cell-centred (clean-cc) state, for a momentum comparison

USER (11-22): Orion is measuring whether the CELL-CENTRED clean scheme conserves radiation momentum better than the staggered
one (stag2). Please post, as a NOTE on this branch:
1. The newest clean-cc commit you consider runnable (branch + sha, pushed?), and its recommended production keys
   (the cell-layout equivalent of docs/clean_stag2/{hg,rsg}_latest.txt).
2. How it does thin transport / realizability (limiters, clips, cuts, AP-HLL?), and whether radiation momentum telescopes.
3. Its known state on RSG/HG RESTARTS — USER recalls it never ran cleanly there: what failed (CW, stalls, dt collapse,
   realizability), on which commit, and any notes/run dirs.
No runs needed on viper; text + pointers are enough. Orion's agent CELLEVAL builds and runs it on Orion/Raven.

-- orion (coordinator)
