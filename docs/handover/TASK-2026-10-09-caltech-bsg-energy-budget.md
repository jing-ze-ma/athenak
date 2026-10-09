# TASK viper -> Caltech: energy-budget check in every BSG link NOTE (user 10-09)

Question from the user: does the new radiation scheme give the correct output luminosity for BSG? Page status at 6.0 d
(unseeded run): flux at each column's own tau_R = 1 = 1.021 L (5.5 d), but L_top (Marshak top face, 80 Rsun) mean over
0-6 d = 0.93 L (range 0.91-1.19, start-up transient). Settle it as the seeded run relaxes. In EVERY link NOTE of the
seeded run (and once now for the unseeded prod_noseed at its last dump), report:
1. L_top / L at the latest hst rows (mean over the last 0.5 d and over the last 2 d), L = 5.581e38 erg/s.
2. L_in check: 4 pi r_in^2 F at the inner face (from hst or the first-cell face flux) / L.
3. L(r) / L from the newest dump, solid-angle mean, at r = 25, 30, 35, 40, 45, 50, tau_R=1 (own column), 55, 60, 70,
   79 Rsun: radiative part = lab-frame F_r minus (4/3) E v_r (comoving-equivalent), plus the gas enthalpy + kinetic
   advective part (5/2 P_gas v_r + 1/2 rho v^2 v_r), and their sum (the page's Fig.-5c definition). Above the
   photosphere use the FACE luminosity if the run writes m1_face (cell-centred F in thin cells is not a transport flux);
   say which you used.
4. Budget over the link: L_in - <L_top> - dE_tot/dt, with E_tot = internal + kinetic + gravitational + radiation
   energy of the domain (hst columns), plus the boundary mass/enthalpy flux at both walls if non-zero. Give it as a
   fraction of L.
Pass criteria (in steady state, t > 30 d): |L_top/L - 1| <= 0.03 averaged over 2 d (convective variability aside),
budget residual <= 1 % of L, L(r) flat within the convective noise between 25 Rsun and the photosphere and continuing
to the top without a drop > 2 %. Flag at once if after ~10 d L_top stays < 0.97 L while L at tau=1 is ~1: that would
point at the thin atmosphere (53-80 Rsun) / the half-range top faces, which viper will then investigate.
