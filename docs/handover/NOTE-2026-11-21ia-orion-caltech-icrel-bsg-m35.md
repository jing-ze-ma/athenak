# NOTE 2026-11-21ia orion -> caltech: 3-D check of the REBUILT ICs (viper builder) — BSG + m35, on B_j3

Files in docs/handover/caltech_runs_1121/icrel/: bsg.athinput, m35mms.athinput (start decks: viper ic_star column, ideal gas, viper floors/sponge/vceil/outflow,
he_gm_column, sph1d re_atm hse + he_ic_balance + he_ic_equil, clean_init dtfac 1.0 n 40, flim + thinf_held + floorcap, seed n2 1e-3),
and the gzipped IC files (gunzip; point problem/he_ic_file at your copy — the decks hold an Orion path). Tables = your fresh_prod_1120 ones.
Run on B_j3 (a5257b7ee), production layout, one node each, free test nodes:
  I1_bsg : bsg.athinput, tlim 1.12e6 s (13 d)
  I1_m35 : m35mms.athinput, tlim 3.0e4 s
Report every ~20 min: the startup lines (atm march, clean_init steady, he_ic_equil a_net/g, |L/L_in-1|), dt, KE, L_bot/L_in, L_top/L_in,
CW/FATAL, FLIM/floorcap lines, never-read list. Results to caltech_results/I1/. T19/T9 follow when their setups pass.

-- Orion coordinator
