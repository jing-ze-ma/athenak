# NOTE 2026-11-21ih orion -> caltech: the J4 runs ARE the new fresh PRODUCTION (USER) — paper tlim, production rules

USER: the four J4 runs (NOTE ig: fresh from t = 0 on B_j3 joint3 a5257b7ee, *_wngl_smoke decks -> use the *_wngl_fresh production decks
with the same added keys: rad_m1/clean_init_dtfac = 1.0e6, rad_m1/clean_stag_flim = true, rad_m1/clean_stag_thinf_held = true,
rad_m1/floorcap_rho = 1e-14, rad_m1/floorcap_v = 3e7, rad_m1/floorcap_report = 25) are the new production. The USER will confirm in your
session as required. Please run them as production:
  - paper tlim: BSG 57 d (4.9248e6 s), m35mms 10 d (8.64e5 s), T9L6.2 745 h (2.682e6 s), T19L6.4 350 h (1.26e6 s);
  - production rules: one EXCLUSIVE node per star (no other steps on it), prodg.sh with chunk wall ends before the holder ends, KEEPRST 3,
    the deck output cadences (hst/bins/rst/m1_vet as in the _fresh decks), chunks.txt with binary md5 + deck md5;
  - if J4 was already started from the smoke decks with tlim per NOTE ig, either restart now from t = 0 with the _fresh decks (preferred,
    still early) or continue from its rst with the _fresh deck's outputs and the paper tlim — say which.
  - successor holders: the current holders end 2026-11-21 17:46 PST; plan the chunk chain (the USER may need to OTP).
Progress lines every ~20 min as now. Results/monitor as for the 11-20 production.

-- Orion coordinator
