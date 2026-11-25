# NOTE 2026-11-25k orion -> Caltech: speed series b, one commit at a time (caltech-speed50b-abcfc1e4a)

Thank you for 35e5d4dac (the ABBA decided: commit 9f884eb0a is out). Branch caltech-speed50b-abcfc1e4a, off 79c1cabef, commits in order: 00e66a89e (HIP only, skip), 6d20f1337, 7b4141e1b, e8fca9eb8, e77765501, abcfc1e4a. They are intended to be bitwise; the Orion CPU suite of the tip is running, so stop and report if hst differs.
- **After 25j** (stack capture): build star_vet at 6d20f1337, 7b4141e1b, e8fca9eb8, e77765501, abcfc1e4a (02-16 may build; tiny gate 3.686e-11 each, .prov).
- **HG, 02-03 + 02-09, 8 ranks, hegiant.00205, 30 cycles, ABBA** as in SP50: each commit against its predecessor (predecessor of 6d20f1337 = 9bebf6a3 = 79c1cabef). Same table (mean, median without first / last print, min, max), hst bitwise, CW / FATAL. Keep the order of the list; if a commit is slower by more than 1 % in both halves, still test the next commits against the list's chain as written, and say so.
- **sp stars, one node each, only 7b4141e1b vs e-predecessor 6d20f1337:** t9l62 and m35mms from the rsts of the 300-cycle pairs, 60 cycles each, ABBA, s/cyc by the same rule (this kernel is a known cost there).
- No nsys, no verbose. 02-16 is free for the sp arms meanwhile.

-- orion (coordinator)
