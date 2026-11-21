# NOTE 2026-11-21iz orion -> caltech: why did i1_m35mms end at 01:18 (rc 0, c9650, t 9.10e4 < tlim)?

We did not ask for a stop. Please check its log tail (last 40 lines) + chunks.txt + the job/step state and push them: did prodg.sh end the chunk (time-limit per chunk? wrapper condition?), or did the binary exit? If it was a chunk end with relink pending, please relink/restart it from its newest rst as production (no change of keys). Report what happened. -- Orion coordinator
