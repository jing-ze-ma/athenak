# NOTE 2026-11-21ii orion -> caltech: the two rebuilt-IC runs (I1_bsg, I1_m35) are ALSO production (USER)

USER: I1_bsg and I1_m35 (rebuilt ICs, decks docs/handover/caltech_runs_1121/icrel/{bsg,m35mms}.athinput with clean_init_dtfac 1.0e6 + flim + thinf_held + floorcap) are production too, alongside the four J4 runs. Same production rules as NOTE ih: paper tlim (BSG 57 d = 4.9248e6 s, m35mms 10 d = 8.64e5 s), one EXCLUSIVE node each, prodg.sh chunks before the holder ends, KEEPRST 3, chunks.txt with md5s, the production output cadences of the matching _fresh decks (hst, bins, rst, m1_vet), progress every ~20 min. If already running, continue from their state with the paper tlim and production outputs (or restart from t = 0 if that is simpler while still early) — say which. The USER will confirm in your session. That makes 6 production nodes; the C2F test arms use the remaining 2.

-- Orion coordinator
