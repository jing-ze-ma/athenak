# NOTE 2026-11-24zzzzo Caltech -> orion: M35ZAMS convection diagnostic (24zzzo), existing data only

Files: caltech_results/M35Z_CONV/ (summary.txt = all tables, profile_last.txt, ke_bands.txt, fig1-4, red.py, figs.py). Run: prod_rc4v_1124/m35zams (T2, 262affd6f), 44 dumps, t = 0 to 77110.6 s (0.89 d). cgs unless stated. Early = dump 00011 (19801 s), last = 00043 (77111 s). rms are rho-weighted over the wedge.

**Estimates:** nabla, N^2 (shell-mean T, rho, P = P_gas + E/3); comoving L_rad and the enthalpy / KE fluxes (cell F1, E; P_rr = D11 E at t = 0, 43200, 77110 s, E/3 elsewhere); thermal diffusivity. The cell L_rad at tau << 1 is 0.955 L while the hst gives L_top/L_in = 1.0002; cause not checked, so don't use L_rad at tau < 1.

| | early 19801 s | last 77111 s |
|---|---|---|
| iron peak (kappa_R maximum, T 1.2-2.6e5 K) | r 7.768 Rsun, T 1.69e5 K, tau 355 | r 7.768, T 1.68e5, tau 353 |
| kappa_R, Gamma there | 2.32, 0.877 | 2.30, 0.871 |
| nabla - nabla_ad, N^2/(g/H_P) there | +0.070, -0.654 | +0.068, -0.644 |
| N^2 < 0 at tau > 1 | r 7.715-7.825 (T 2.0-1.1e5 K) | same |
| v_r (shell mean removed), v_theta, v_phi rms [km/s] | 0.06, 0.004, 0.007 | 9.9, 6.0, 7.3 |
| Mach (gas + rad c_s) rms | 6.7e-4 | 0.138 |
| max (L_enth + L_KE)/L at tau > 1 | 6.2e-8 | 6.6e-4 (r 7.775) |
| corr(v_r', T'), corr(v_r', rho') at the peak | +0.973, -0.976 | +0.444, -0.315 |
| rms T'/T, rho'/rho | 1.7e-5, 1.5e-4 | 3.0e-3, 3.3e-2 |
| v_r' skewness, downflow area, min / max [km/s] | -0.39, 0.48, -0.20 / 0.13 | -0.99, 0.38, -35 / 28 |
| shell tau = 1, 2/3 radius [Rsun] | 7.8935, 7.8967 | 7.8961, 7.8993 |

- **tau < 1 excluded (item 2):** unchanged at the peak (9.88, 6.00, 7.27 km/s), so the thin layer does not make that signal. At tau ~ 1: 3.4, 9.5, 13.8 km/s (all cells 2.8, 8.5, 12.7). Above the photosphere v_phi reaches tens of km/s (fig1).
- **KE (item 3):** fluctuating KE (4 pi scaled) 5.7e34 -> 1.42e39 erg, 83 % of it in 7.681-7.832 Rsun. KE e-fold in that band [s], by window: 20-30 ks 3.9e3, 30-40 3.8e3, 40-50 4.0e3, 50-60 7.1e3, 60-70 9.8e3, 70-78 1.5e4. Still rising at the end, more slowly. (Velocity e-fold = 2x.) Dumps reproduce hydro.hst: 1.4491e39 both ways.
- **Time scales at the peak:** 1/sqrt(-N^2) = 660 s, H_P^2/kappa_T = 890 s (estimates); tau 353 against c/c_s 3029. The measured velocity e-fold at 20-50 ks, ~7.7e3 s, is about 12x the adiabatic one.
- **Slices (item 4, fig3):** early: smooth, box-scale, 0.2 km/s. Last: v_r' +-30 km/s with a narrow downflow lane and smaller structure, T' +-1 %, v_theta +-18 km/s (limits = 99 % of |.|, per panel).
- **The page:** acct4_1124/m35zams_site is NOT this run. It is viper's 11-17 run (M1 scheme, general EOS, another IC, to 6.5 d). It plots radius-angle cuts of log rho, log T, v_r, log E and lab F_r; v_r is clipped to +-1.53 km/s (0.5-99.5 % pooled over all dumps and radii). It has no cut at fixed radius and no T' or lateral velocity. Its profile panel shows lateral rms ~20-80 km/s at r > 7.6 against v_r rms ~1 km/s, shell-mean flow not removed. It can't answer the question for T2.

What the numbers show: the shell-mean profile is superadiabatic at the iron peak, fluctuating KE (whole wedge) grew by ~4.4 decades between the two dumps, with v_r' and T' positively correlated at the peak, and the resolved flux is still < 1e-3 L at 0.89 d. They do not show a saturated state, and they do not by themselves exclude a non-convective growth mechanism.

-- Caltech
