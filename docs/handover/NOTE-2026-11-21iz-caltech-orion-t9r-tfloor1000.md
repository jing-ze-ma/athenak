# NOTE 2026-11-21iz Caltech -> orion: t9r production restarted with T floor 1000 K (USER) + why

USER saw grid-scale artifacts in t9r (rebuilt IC) that t9l62 (original IC) does not have at the same t = 16.88 d: one-cell hot blobs in T_gas outside the photosphere (zoom 0.25-1.5 R_ph), fine speckle + sharp radial striations in log rho - <log rho>, and regular angular bands in the outer domain (> ~5.5 R_ph). Also L_bot/L_in rings (std 0.19-0.24 since day 9, vs 0.07-0.11 in t9l62).
Deck differences t9r vs t9l62 (beyond the IC): tfloor 3100 K vs 1000 K (5.2 % of cells at the floor), dfloor 1e-20 vs 3.6e-19, he_bc_outer outflow + radial sponge (r0 5.07e13) vs copy, vceil 1e8 vs 3e8, he_opac_logd_min -21, he_gm_column.
USER decision: change the T floor to 1000 K and restart from the current rst. Done ~15:40 PST: t9r from t9l62.00034 (written 15:31 under cfl 0.6 + kp keys), B_fwbs, kp keys, cfl 0.6, + hydro/tfloor = 123359972487.0939 (= 1000 K, same code-unit factor as t9l62). rst hard-linked in t9r/rst_backup_f1k/. Running (c 26800 at 15:42). dfloor / sponge / vceil unchanged (ask the USER if you want those aligned too).

-- Caltech
