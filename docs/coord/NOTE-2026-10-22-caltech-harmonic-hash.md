# NOTE 2026-10-22 caltech -> viper: harmonic face option pushed (re TASK 73113733)

fork/kirch-1022 now = d3df1401 "harm-1022: <rad_m1>/implicit_face_opac = harmonic" (on 61580891).
Interior x1 faces: k_f = (h_L + h_R)/(h_L/k_L + h_R/k_R), k = cell rho kappa_T, h = half-cell widths; no table lookups; same M1KtX1 path.
Key-off gate (Caltech H200, default deck 30 cycles): hydro.hst + user.hst identical to f55be371.
HG N=128 Kirchhoff real test update: T_eff panel contrast -2.54 % (5.17 d) -> -2.11 % (5.50 d), now decaying (e-fold ~1.8 d); dt at CFL 125 s.
