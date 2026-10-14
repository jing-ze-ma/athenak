# NOTE Raven -> Caltech (10-14): He giant Mdot_top = 0 is expected (no bug)

- Mdot_top is computed in HeStarHist (src/pgen/he_star_m1.cpp:3303, FIN4 5c504e49). It is the outer-face IDN flux x area, split by sign:
  Mdot_top (hst col 17) = fm > 0, Min_top (col 20) = fm <= 0. There is no threshold.
- The He giant top is x1max = 200 Rsun with he_bc_outer = outflow (ghost v1 = max(v1_edge, 0)), and the sponge starts at 180 Rsun.
  R_ph ~ 75 Rsun, so nothing outflowing has reached the edge yet. Only floor gas falls in: Min_top ~ -2..-6e8 g/s.
- Raven smoke/bg user.hst (9 files) also give Mdot_top = 0 in every row. (Header columns are 1-based; awk $16 is M_tot.)
- Action: no change. Monitor Min_top next to Mdot_top. Mdot_top switches on once ejecta reach ~180-200 Rsun with v_r > 0.
  An earlier mass-loss signal can be taken from the dumps afterwards (rho v_r r^2 at ~150 Rsun); no mid-production code change.
