# NOTE 2026-11-21iq orion -> caltech: rebuilt T9/T19 become PRODUCTION on the 2 test nodes (USER); C2X moves to viper

USER: the rebuilt Yan-Fei stars stay on Caltech and can be production. So the 2 test nodes become production nodes (one run per node,
no sharing), and C2X (NOTE il) is dropped here (it runs on viper instead):
  P_t9r  : rebuilt T9L6.2 (I1_t9 deck, B_j3). Continue from I1_t9's final rst (t = 2e4) if it was written, else start fresh; tlim = paper
           745 h = 2.682e6 s; production outputs/rst cadence as t9l62 in 1121b.
  P_t19r : rebuilt T19L6.4 (I1_t19 deck from NOTE in: icrel/t19l64.athinput + ic_t19l64_thin4m_dfl19.txt.gz, he_ic_rad_atm_relax 0.25), B_j3,
           fresh, tlim = paper 350 h = 1.26e6 s, production outputs/rst cadence as t19l64 in 1121b.
Prepare both run dirs + launch commands now; LAUNCH ONLY after the USER confirms in your session (production rule). Then add both to the PROD
monitor lines. Holders end 11-21 17:46 PST like the others.

-- Orion coordinator
