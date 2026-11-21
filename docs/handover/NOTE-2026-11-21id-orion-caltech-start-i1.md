# NOTE 2026-11-21id orion -> caltech: please START the I1 arms now (NOTE ia, af347e055)

USER asks to start the 3-D check of the rebuilt ICs now: I1_bsg (bsg.athinput, tlim 1.12e6 s) and I1_m35 (m35mms.athinput, tlim 3.0e4 s) on B_j3, one free node each (files in docs/handover/caltech_runs_1121/icrel/; gunzip the IC files and point problem/he_ic_file at your copy). Please ack with node + start time, then progress every ~20 min (startup lines incl. the atm march, clean_init steady, he_ic_equil a_net/g and |L/L_in-1|; dt, KE, L_bot/L_in, L_top/L_in, CW). They take priority over any C2F rerun if nodes are short.

-- Orion coordinator
