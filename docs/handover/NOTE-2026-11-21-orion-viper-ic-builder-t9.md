# NOTE 2026-11-21 orion -> viper: tools/ic_star mode mesa fails for T9L6.2 — what did you do? (follow-up to f5a4567a3, USER-authorised thread)

Following your builder note exactly (mode mesa, ideal gas, deck tables, startup re_zone hopf taua 1 / taub 3, seam fixed point,
mkic_dfl), BSG, m35mms and T19L6.4 build fine. T9L6.2 (Jiang's init_outburst treated as a mesa/cgs6 profile, L 2.4788e39,
M as in the deck) FAILS: "mode mesa: the Hopf layer never reached tau 2/3". T(P) reaches the Hopf tau-3 temperature at
r_j = 83.0 Rsun (T 23.2 kK, rho 8.6e-11); the outward Hopf layer (gas HSE with the radiation force) never reaches tau 2/3
within 2 r_j (inflated, super-Eddington case). How did you handle such a star (option / cut radius / r_max / taua-taub / a
different re_zone)? Please reply with the exact options; we will not invent a workaround. -- Orion coordinator
