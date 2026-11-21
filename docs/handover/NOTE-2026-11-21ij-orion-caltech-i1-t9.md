# NOTE 2026-11-21ij orion -> caltech: 3-D check of the rebuilt T9L6.2 IC (one test node)

Files: docs/handover/caltech_runs_1121/icrel/t9l62.athinput + ic_t9l62_thin4m_dfl20.txt.gz (gunzip; point problem/he_ic_file at your copy). Rebuilt with viper's ic_star
(Jiang's init_outburst with the 7-digit-r noise rows removed + mass fix), radiative base, R_ph 88.0 Rsun; deck keys: clean_init_dtfac 1e6,
clean_stag_flim, clean_stag_thinf_held, floorcap 1e-14/3e7/25, he_gm 7.4032825532e27, dfloor 1e-20.
Run I1_t9 on B_j3 (joint3 a5257b7ee), production layout, one free TEST node (not a production node), tlim 2.0e4 s. Progress every ~20 min:
startup lines (atm march, clean_init steady, he_ic_equil), dt, KE, L_bot/L_in, L_top/L_in, CW, FLIM/floorcap. Results to caltech_results/I1/.
If it is healthy the USER may make it production like the other rebuilt runs (that needs the USER in your session).

-- Orion coordinator
