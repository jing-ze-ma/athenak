# S5 derivation (2026-10-09): Plaskett progenitor at P = 3.81 d, physical ambient, hot optically thick stream

User decisions 2026-10-09 (design doc section 11): (1) physical ambient (floor gas at radiative
equilibrium, absorption on), (2) the hot optically thick one-zone stream with its radiation injected,
(3) the initial period 3.81 d of Wade et al. 2026 Sect. 5. Code units as before: Rsun, km/s,
Rsun/(km/s) = 6.957e5 s, density unit rho_u = 5.11e-8 g/cc (unchanged, so arad and the other
<rad_m1> constants stay the same).

Scripts (this directory) and their outputs (run directory `accretor_1008/rhd/s5/setup/`):

| script | output | what |
|---|---|---|
| plaskett_setup_1009.py (`OUTDIR=. NX1=640 PDAY=3.81`) | setup_1009.out, faces_nx640.txt | binary, L1, fan, photosphere, grid, r_meas, wb depth |
| s5/setup/traj_launch.py | (stdout) | central trajectory vs the launch speed; curvature at the r_out crossing |
| stream_physics_1009.py EOS ROSS PLANCK | stream_physics_1009.out | L1 state (adiabatic Q), one-zone core at r_out, window state |
| make_ic_accretor_column.py ... a_sep=33.2596 | column_z0.02_p381.txt/.out | S0 column (Z 0.02 tables, unchanged) |
| band_luminosity.py 33.2596 | band_1009.out | von Zeipel band factor |

## 1. Inputs and approximations

| input | value | source / approximation |
|---|---|---|
| P_orb | 3.81 d | Wade+2026 Sect. 5 (initial period; their Discussion says "a 3.7-d orbit") |
| M_gainer, M_donor | 16, 18 Msun | unchanged (Wade+2026 Sect. 6.4 "M ~ 16 Msun, R ~ 9 Rsun"; donor 18.2 less winds) |
| R_acc, T_eff,gainer | 9.00129 Rsun (a grid face), 28.3 kK | unchanged (HR figure read) |
| donor T_eff | 26.2 kK | unchanged (HR figure read) |
| donor R | 12.943 Rsun = R_L | APPROXIMATION: the donor fills its Roche lobe (Eggleton 1983 at q 1.125); at fixed T_eff its L rises to log 4.852 (env13 4.84) |
| Mdot | 1e-4 Msun/yr | Wade+2026 Sect. 6.4 |
| L1 rate factor | 0.649 | Ryu+2025 Eq. 4, ADIABATIC case (optically thick overflow; the L1 gas comes from tau ~1.6e3) |
| spin | 1 (corotation) | unchanged, our choice (Wade+ models are non-rotating) |
| tables | X 0.70 Z 0.02 GS98 ext2 (md5 ad516c7e / 99d668f6) | unchanged; copies in s5/setup |

## 2. Geometry (setup_1009.out)

| quantity | 3.81 d (S5) | 3.69 d (env13) |
|---|---|---|
| a | 33.2596 Rsun | 32.5575 |
| Omega | 1.908715e-5 /s = 13.278932 code | 1.9708e-5 |
| P_orb (code) | 0.473169 | 0.458266 |
| a Omega | 441.65 km/s | 446.39 |
| d_L1 | 16.22664 Rsun | 15.8841 |
| r_out = 0.85 d_L1 (x1max) | 13.7926 | 13.5015 |
| L1 curvature A, B, C | -16.988, 6.994, 7.994 (q only) | same |
| photosphere r_ph (phi 0 / 90 / 180) | 9.4452 / 9.0013 / 9.2841 | 9.4847 / 9.0013 / 9.3040 |
| H_p,min (c_ph 19.41, env13 criterion) | 0.01020 Rsun -> dr <= 2.551e-3 | 0.01021 |
| envelope top psi = -15 c_ph^2 (phi 0 / 90) | 9.6412 / 9.1575 | 9.6855 / 9.1577 |
| grid (nx1 640, plateau stretch, c_k = 0) | fine zone 8.756..9.598 (340 cells, dr 2.474e-3), face 218 = R_acc, ratio <= 1.095, dr(r_in) = dr(r_out) = 2.82e-2 | 8.76..9.638, face 215 |
| stretch keys | p_amp -2.1938589633103867, p_xa 0.13152342165045811, p_xb 0.77203350400482496, p_w 0.01875 | |
| r_meas | 9.572888 (face 449 >= max r_ph + 10 H_p) | 9.615675 |
| env_wb_depth | 0.396 (1.5 d_pen 0.264) | 0.41 |
| tlim (10 orbits), bin dt P/20, rst dt P/10 | 4.73169, 0.023658, 0.047317 | 4.58266, 0.022913, 0.045827 |

## 3. S0 column and the bottom flux

Same EOS dump and Rosseland table, a_sep = 33.2596 (only the phi = 90 Roche cut changes):
L = 4 pi R_acc^2 sigma T_eff^4 = 4.6822e4 Lsun (unchanged), rho_ph 6.0277e-10 g/cc (3.69 d: 6.0196e-10),
Gamma_edd(ph) 0.350, psi(r_in) 1.43565e5 (km/s)^2, HSE residual 4.1e-6; Schwarzschild-unstable zones
8.961-8.993 and 8.585-8.791 (as before; radiative column kept).

Band factor (band_luminosity.py at psi = 0): **0.963286** (3.69 d: 0.960504; the lobe is less distorted
at the wider orbit). F_col(r_in) = 7.4249e13 erg/cm^2/s = 1453030.13 code, so
**implicit_flux_x1min = 0.963286 x 1453030.13 = 1399683.58 code** (L 45103 Lsun at r_in, phi 90).
3-D: factor 1.

## 4. The stream state at r_out (stream_physics_1009.out)

1. **L1** (donor's grey atmosphere, EOS + Rosseland table, g = G M_d/R_L^2 = 2.95e3; Kolb & Ritter form
   Mdot = Q 2 pi/(Omega^2 sqrt(BC)) int c_T dP_gas with Q = 0.649): the overflow reaches tau 1646:
   T_L1 155.3 kK, rho_L1 5.25e-8 g/cc, P_rad/P_gas 1.35, c_T = sqrt(P_gas/rho) 45.45 km/s,
   sigma_y,z 1.21 Rsun. APPROXIMATION (predecessor's model, kept): the integrand uses the gas pressure
   and the isothermal c_T although P_rad/P_gas ~1.35 at L1; a K&R integral over P_tot would give a
   shallower L1 (lower T_L1) for the same Mdot. With g x0.3 near L1: 253 kK (bracket).
2. **Launch** (new, for consistency with 1): the ballistic stream leaves L1 at the L1 overflow speed
   c_T = 45.45 km/s (problem/t_don = 155115 K at mu 0.62; t_don only sets the launch speed in the pgen),
   not at the donor-T_eff 18.68 km/s. Central trajectory (traj_launch.py): at r_out phi 3.410 deg,
   v_r -126.46, v_phi 60.68, |v| 140.27 km/s (25.6 deg from radial), flight 0.0348 code (2.42e4 s);
   impact on the photosphere at phi 15.07 deg, |v| 410.7 km/s. (18.68-km/s launch: 3.798 deg,
   -121.64/54.94, 0.0522 code, impact 14.75 deg.)
3. **One-zone core** L1 -> r_out (adiabatic work on gas + trapped radiation, radiative diffusion loss
   a(T^4 - T_eq^4)/t_diff with t_diff = 3 kappa_R rho s^2/c, T_eq = 19.2 kK the thin equilibrium T
   with both stars; pressure-supported widths s = c_iso/(Omega sqrt(Phi_nn)), c_iso^2 = P_tot/rho,
   Phi_nn 3.6132, Phi_zz 9.222 Omega^2 at the r_out crossing; Mdot conserved):

| quantity | S5 window value |
|---|---|
| T (uniform across the window) | **64767 K** (env_t_stream) |
| rho peak | 1.2305e-8 g/cc = **0.240807 code** (rho_stream) |
| c_iso (P_tot) | 34.96 km/s |
| sigma_perp, sigma_z | 1.3850, 0.8669 Rsun |
| arc width sigma_perp/cos 25.6 deg | **1.5362 Rsun** (stream_width; sigma_phi 6.38 deg; nsig 3 window +-19.1 deg about 3.41 deg) |
| v at r_out | (-126.46, 60.68) km/s (from the pgen's own ballistic integration) |
| P_rad/P_gas (peak) | 0.419 (2 sigma: 3.1, 3 sigma: 38) |
| kappa_R, kappa_P (peak) | 2.72, 349 cm^2/g |
| tau across one sigma_perp / to the surface from 3 sigma | 3.2e3 / 11 |
| t_diff / t_flight | 1.8 |
| injected E = a T^4 | 1.3313e5 erg/cc = **260.52 code** |
| injected F = (4/3) E v (lab) | 2.49e12 erg/cm^2/s = 48725 code; F_r 2.24e12 = 0.145 of the stellar flux at r_out; F/(cE) 6.2e-4 |

   Bounds: pure adiabat 75.0 kK, 0.139 code, arc width 2.02, P_rad/P_gas 1.13.

4. **Profile across the window: APPROXIMATION.** One zone means one T: the window has T uniform and
   the Gaussian density of an isothermal pressure-supported stream. A real cooling stream has T falling
   toward its surface (diffusion), so the uniform T overstates the wing radiation: the injected
   radiation energy flux over the 6-sigma window is 0.96 of the gas enthalpy flux (0.35 of the kinetic),
   about 2.4x what a Gaussian-weighted T^4 would carry. The wings (> 2 sigma) carry 4.6 % of the mass.

## 5. Injection (problem/stream_rad = true)

Window ghosts (r > r_out, |phi - phi_s| < 3 sigma): gas rho_s exp(-dphi^2/2 sigma^2), v_s, T_s (EOS,
gas only, as before); M1 ghosts E = a T_s^4, F_r = (4/3) E v_r, F_phi = (4/3) E v_phi (optically thick
advection, comoving flux 0). The implicit x1 solve does not read x1 ghosts (face-flux BC), so the same
state enters the solve as a per-column incident bath of the outer Marshak face (module key-free API
RadiationM1::SetX1maxBathColumns): comoving face flux c q (E_ie - E_bath), inflow enthalpy flux
(1 + chi) v E_bath = (4/3) v E_bath in the thick limit (implicit_bc_advect). Elsewhere the bath is 0
(vacuum/Marshak as before). The advective radiation energy through r_out is booked in Ebnd.
vet_col: its formal solution still enters from a vacuum top; with tau ~1e2 per radial cell in the window
core, that affects only a skin of the top cell (closure error confined to tau < 1).

## 6. Ambient (problem/env_amb_mode = radeq)

- Density: the floor rho_amb = hydro/dfloor = 1e-8 code = 5.1e-16 g/cc everywhere outside the column's
  atmosphere (IC: the column ends where its density falls to rho_amb). For reference the gainer's wind
  (Mdot ~1e-7 Msun/yr, 2000 km/s) has ~5e-15 g/cc at 10 Rsun: the floor is 10x below it.
- T: the radiation temperature. IC gas T = T_col (clamped at the column top, 23.8 kK) = the radiation
  IC's T; afterwards the M1 coupling (absorption ON, Planck table incl. its low-rho extension) holds
  T_gas = T_rad (coupling time ~0.1/kappa_P s vs dt ~1 s). No T relaxation (env_t_relax = 0).
- Dynamics: not hydrostatic (H_p ~0.01 Rsun at 2e4 K); it falls, the existing floor sponge (rho < 10
  rho_amb: velocity damping 1e-4, |v| <= 100 km/s) keeps it slow; KE removed is booked in Erel.
- Budget: the pgen applies the floor itself (after the sources; added gas at the cell's v and T): mass
  in history column 8 "Mfl" (cumulative; replaces Jstr, the wall-stress diagnostic that has no meaning
  deep in the envelope), energy in Erel. Mass budget dMdom = Min - Mout - Mwal + Mfl.
- Size: P_floor/P_ph ~1e-6; ram of floor gas falling at 100 km/s ~1e-4 code vs P_ph ~7; KE flux ~1e-6 of F.
