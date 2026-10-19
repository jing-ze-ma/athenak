# Plaskett accretor setup, 10-07: derived inputs

Script: `plaskett_setup_1007.py` (same directory). Full output: `out.txt`. Faces of the chosen grid: `faces_nx500.txt`.
`fan_cache.npz` caches the ballistic fan; delete it to recompute (about 70 s).
`NX1=640 python3 plaskett_setup_1007.py` prints the mesh keys for another cell count.
Constants are the pgen's. Code units: Rsun, km/s, Rsun/(km/s). Code density is in units of the stream peak at r_out.

## Final key values

| key | value | formula / source |
|---|---|---|
| problem/env_cs_ph | **19.4132** | sqrt(k 28300/(0.62 m_H)) |
| problem/env_cs_stream | **18.6791** | sqrt(k 26200/(0.62 m_H)) |
| problem/t_don | **26200.0** (mu_don 0.62) | the pgen derives c_s,don = 18.679 from it |
| problem/a_sep, period | 32.5575, 3.69 | Kepler; unchanged |
| problem/r_acc, env_r_spin | **9.00129260997162** | unchanged; face 153 of the new grid (exact) |
| problem/stream_width | **0.735** | sigma_perp/cos(alpha) at r_out (sec. 2); sigma_phi = 0.05444 rad = 3.12 deg |
| problem/env_rho_ph | **0.055** (range 0.042-0.055; +-2x from P_ph) | rho_ph/rho_peak(r_out) = 2.819e-9/5.11e-8 |
| problem/env_r_top | **9.68728512526772** (NEW, needed) | first face above the phi = 0 envelope top 9.6855 |
| env_wb_depth (new key, not in the source) | **0.41** (conservative 0.55) | 1.5 d_pen; today's equivalent is hydro/wb_rmax = 8.59 (8.45) |
| problem/env_amb_rho | **1e-7** (recommended) | P_amb/P_ph 4.3e-4 |
| hydro/dfloor = problem/rho_amb | **1e-8** | keeps rho_amb/floor = 10 at R_acc |
| mesh (nx1 500) | see sec. 6 | plateau stretch, fine zone 8.80-9.638 at dr 2.477e-3 |
| time/tlim | **4.58266** | 10 x 0.458266 |
| output2/dt (bin) | **0.022913** | P/20 |
| output3/dt (rst) | e.g. 0.045827 | P/10 (the 0.0485 in the input is from RY Per) |
| r_meas (analysis) | **9.615674** (face 401) | max r_ph + 10 H_p(phi 0) = 9.6150 |

## 1. Sound speeds, binary, L1

- **Sound speeds:** c_ph = 19.413 km/s and c_s,don = 18.679 km/s. Both are isothermal, with mu = 0.62.
- **Binary:** a = 32.5575 Rsun and Omega = 13.710767 (1.9708e-5 /s), so P_orb = 0.458266 code = 3.69 d. a Omega = 446.39 km/s.
- **L1:** d_L1 = 15.8841 and r_out = 0.85 d_L1 = 13.501495, written as x1max 13.5015.
- **Curvature at L1:** Phi = (1/2) Omega^2 (A x^2 + B y^2 + C z^2). These are numerical second derivatives of the code's potential:
  - A = -16.988 (the Ryu fit gives -16.79)
  - B = 6.994 = -(A+3)/2
  - C = 7.994 = B + 1
- **Widths at L1** from rho_L exp(-Omega^2 (B y^2 + C z^2)/(2c^2)):
  - sigma_y = c/(Omega sqrt B) = 0.5151, times 0.932 (Ryu) = **0.4801**.
  - sigma_z = c/(Omega sqrt C) = **0.4819**, uncorrected. This is my choice: Ryu+25 found that vertical HSE holds and quote no z factor.
  - The peak factor 0.935 is used only for the rho_L report.

## 2. Ballistic fan (3-D restricted three-body problem, Coriolis included)

**Setup.** 121 x 31 launch points cover y in +-4.5 sigma_y and z in 0..4.5 sigma_z, mirrored to z < 0. Each is weighted by the Gaussian flux rho_L c dy dz. The launch is the pgen's: from x_L1 - 1e-4 at speed c toward the accretor. RK4 with h = 2e-5/Omega, integrated to the sphere r = r_out.

**Launch check.** 127 of 3751 trajectories fall back to the donor, all at y < -4 sigma_y. The flux fraction that reaches r_out is 0.99999.

**Central trajectory at r_out** (what the pgen uses): phi 3.799 deg, v_r -122.93, v_phi 55.47, |v| 134.87 km/s. The flight time is 0.0507 code.

**Flux-weighted at r_out:** <phi> = 4.100 deg, a shift of +0.30 deg (0.07 Rsun, about 0.1 sigma), with v_r -122.58 and v_phi 54.61. The central orbit is therefore an adequate stream centre.

**The fan's widths are pressureless caustic artefacts:**
- **Fan numbers:** a Gaussian fit to the flux-weighted phi distribution gives sigma_phi = 0.457 deg, so sigma r_out = 0.108 Rsun. The fitted FWHM is 1.08 deg against 0.61 deg from the histogram. The weighted std is 0.968 deg (0.228 Rsun).
- **In-plane folds:** phi(y) folds at y = -3.5 sigma_y and +1.8 sigma_y. The flux piles up on the +1.8 sigma fold at phi 3.2 deg, and the Gaussian fit locks onto that caustic.
- **Vertical:** z(r_out)/z(L1) = -0.39. The trajectories have crossed the midplane, and the fan sigma_z is 0.176.
- **Cause:** the restoring curvature (sqrt B = 2.64 Omega, sqrt C = 2.83 Omega) acting over the 0.70/Omega flight is a phase above pi/2, so pressureless orbits focus.
- **Why the real stream does not:** an isothermal stream is pressure-supported. Ryu+25 report that it expands laterally at 0.2-0.3 c and stays in vertical HSE; it does not focus.

**Adopted: transverse hydrostatic widths at the r_out crossing** (the same definition as at L1).
- **In plane:** Phi_nn = 4.135 Omega^2 across the stream, so sigma_perp = c/sqrt(Phi_nn) = **0.670**. That is +0.19 over L1, consistent with 0.25 c of expansion over the 0.05 flight.
- **Vertical:** Phi_zz = 9.218 Omega^2, so sigma_z(r_out) = **0.449**.
- **Arc width:** the stream is 24.3 deg from radial, so the width along the r_out arc is 0.670/cos 24.3 deg = **0.735** = stream_width.
- **Coverage:** FWHM 1.73 Rsun, 18 phi cells per sigma at nx3 = 2048, and a nsig-3 window of +-9.4 deg.

**One velocity is acceptable.**
- Launches with |y| <= sigma_y: |v| = 134.2-136.6 km/s (+-1%), v_r -131..-115, v_phi 40..70. The direction spread is about +-6 deg.
- All flux (2.5-97.5%): v_r -154..-105 and v_phi -13..81. The tails are the folded wings.
- v_z rms is 17 km/s, below c.

## 3. Density scale

**Photosphere** (viper's formula, uncertain by +-2x): g = G M_acc/R_acc^2 = 5.418e3. P_ph = (2/3) g/kappa_e = 1.062e4, so **rho_ph = P_ph/c_ph^2 = 2.819e-9 g/cc**.

**Peak density at r_out:** Mdot = 6.302e21 g/s = rho_peak |v| 2 pi sigma_perp sigma_z. This is the same as rho_peak |v_r| 2 pi sigma_arc sigma_z.

| widths at r_out (perp x z) | rho_peak(r_out) g/cc | env_rho_ph |
|---|---|---|
| **pressure-supported 0.670 x 0.449 (adopted)** | 5.11e-8 | **0.055** |
| L1 widths kept, 0.480 x 0.482 | 6.64e-8 | 0.042 |
| L1 analytic, 0.515 x 0.482 | 6.19e-8 | 0.046 |
| ballistic caustic, 0.098 x 0.176 (unphysical) | 8.9e-7 | 0.003 |

- **Against earlier estimates:** viper's 0.3 is 5.5x too high, because it used A = (c/Omega)^2 with no 2 pi, 650 km/s and c_ph 20. The simple estimate of 0.051 agrees, because its widths (0.54 x 0.541) times 135 km/s nearly coincide with these.
- **Ryu check:** Mdot/Mdot_an = 0.721 - 0.149 tanh^2(0.522 log q) = 0.7209 at q = 1.125. Setting 2 pi rho_L c sigma_y sigma_z = Mdot gives rho_L = 4.47e-7 analytic, or **6.20e-7 g/cc** with the 0.721 factor (simulation peak 0.935x = 5.80e-7).
- **Plausibility against the donor photosphere** (same formula, R = 12.8, Teff 26.2 kK):
  - rho_ph,don = 1.69e-9 g/cc, so rho_L/rho_ph,don = 366.
  - Isothermal Bernoulli then puts the donor photosphere 6.4 c^2 above Phi_L1. That is about 0.107 Rsun, or 6.4 H_p,don (H_p,don = 0.017 Rsun), a relative overfill of 8e-3. This is plausible for thermal-timescale Case A.
  - The L1 gas is optically thick (kappa rho_L sigma_y about 7e3), so "isothermal at Teff" is only approximate.

## 4. Photosphere (code RochePot, spin 1)

Equipotential through (R_acc, 90 deg): Phi_s = RochePot(9.00129, pi/2).

| phi | 0 | 15 | 30 | 45 | 60 | 75 | 90 | 105 | 120 | 135 | 150 | 165 | 180 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r_ph | 9.4847 | 9.4238 | 9.2862 | 9.1466 | 9.0476 | 9.0005 | 9.0013 | 9.0398 | 9.1026 | 9.1748 | 9.2410 | 9.2874 | 9.3040 |
| H_p | .01302 | .01258 | .01167 | .01089 | .01041 | .01021 | .01024 | .01044 | .01073 | .01107 | .01139 | .01161 | .01169 |

- **Photosphere:** the maximum r_ph is **9.4847**, at phi = 0. H_p = c_ph^2/g_r ranges from 0.01021 to 0.01302.
- **Envelope top** (psi = -15 c_ph^2): 9.6855 at phi 0, 9.1577 at phi 90, 9.4832 at phi 180.
- **INCONSISTENCY (affects plaskett_env12 too):**
  - The pgen's default r_top searches only above (R_acc, 90 deg), giving 9.167 with c_ph 19.4 and about 9.16 in env12.
  - EnvSetup initialises every cell with r > r_top as ambient (`amb || x1v > rtop`), and Phi_wb is flat in r above r_top.
  - The photosphere lies outside r_top for |phi| < 42.5 deg and |phi| > 133.5 deg, so the caps facing L1 and L2 (up to 0.32 Rsun thick at phi 0) are missing at t = 0.
  - Fix: problem/env_r_top >= 9.6855. Use the face 9.68728512526772.
  - Not checked: whether WB on the hot ambient between 9.17 and 9.69 at phi 90 behaves. Its H_amb = 2.4 Rsun is much larger than dr, unlike the cold-ambient failure.

## 5. Penetration depth

**Impact of the central trajectory on the photosphere:** phi 13.95 deg, r 9.4316, |v| 395.5 km/s (rotating frame). The normal speed is v_n = 339.0 km/s, 31 deg from the inward normal.

**Envelope pressure:** P(psi) = env_rho_ph c_ph^2 (1 + psi/(4 c_ph^2))^4, with env_rho_ph 0.055 and P_ph = 20.78.

**Ram pressure:** rho_imp v_n^2.

| rho_imp (code) | ram/P_ph | psi_pen | d_pen at phi_imp | same psi at phi 90 |
|---|---|---|---|---|
| 2-D continuity, rho ∝ 1/abs(v): 134.87/395.46 = 0.341 | 1.89e3 | 22.4 c_ph^2 | **0.272** | 0.223 |
| conservative 1.0 | 5.53e3 | 30.5 c_ph^2 | **0.367** | 0.302 |

**Recommendation:** env_wb_depth = 1.5 d_pen = **0.41 Rsun** (conservative: 0.55). Applied on the equipotential through R_acc - 0.41 at phi 90, i.e. r = 8.593, that switch-off lies even deeper on the impact side, so it is conservative there.

- **Key not implemented:** env_wb_depth does not exist in the source yet. Today's hydro/wb_rmax (a sphere) = 8.59, or 8.45 conservative; the input has 8.52.
- **Resolution at the switch-off** (plain gravity there): the envelope H_p is 0.107 Rsun at 8.593 and 0.138 at 8.451.

| switch-off | H_p/dr, nx1 500 | nx1 640 | nx1 768 |
|---|---|---|---|
| 8.593 | 5.1 | 8.0 | 8.4 |
| 8.451 | 4.0 | 6.4 | 8.7 |

## 6. Radial grid

- **Requirement:** H_p,min = 0.01021 at phi 75-90, so dr <= H_p/4 = 2.554e-3. The fine zone must span 8.80 to max r_ph + 0.15 = 9.635.
- **Cell count:** that zone alone needs about 330 cells, leaving about 160 of 500 for 6.3-8.8 and 9.64-13.5.
- **Why the 10-07 method fails here:** grid_design.py's c1..c8 polynomial fit to a log target left 15% wiggles inside a zone this large; dr reached 2.95e-3 inside it.
- **Design used: a PURE plateau** (c_k = 0). du/dxi is exactly flat inside, where dr_in = 0.97 H_p/4 = 2.477e-3, and flat outside, with tanh edges of 12 cells.
- **Pinning R_acc to a face:** the plateau is then shifted rigidly by a fraction of a cell so face 153 is exactly 9.00129260997162. The pgen fatals unless a face lies within 1e-6 r_acc of problem/r_acc.

| nx1 | zone dr<=H_p/4 | max neighbour ratio | max dr/r | dr(r_in) | dr(R_acc) | dr(9.5) | dr(r_out) |
|---|---|---|---|---|---|---|---|
| **500** | 8.8002-9.6383 (338 cells) | 1.123 | 1.23e-2 | 7.74e-2 | 2.477e-3 | 2.477e-3 | 7.78e-2 |
| 640 | 8.7605-9.6383 | 1.095 | 4.5e-3 | 2.84e-2 | 2.477e-3 | 2.477e-3 | 2.84e-2 |
| 768 | 8.7804-9.6383 | 1.078 | 2.7e-3 | 1.72e-2 | 2.477e-3 | 2.477e-3 | 1.72e-2 |

**nx1 = 500 works, but is coarse outside the zone.** Compared with r dphi of 0.019-0.041:
- dr is 0.077 at r_in and 0.078 at r_out.
- The neighbour ratio is 1.12, against <= 1.081 for the old grid.

nx1 = 640 is cleaner: ratio 1.095 and dr/r <= 4.5e-3. 768 matches the old ratio. Set meshblock nx1 = mesh nx1.

**Keys for nx1 = 500.** mesh.cpp always reads c1..c4, with default 0. c5..c8 and the plateau are read only when present. faces = r0 + (r1 - r0) u(i/nx1).
```
nx1       = 500
x1min     = 6.3
x1max     = 13.5015
use_grid_stretch_r_poly = true
f_stretch_r_c1 = 0.0
f_stretch_r_c2 = 0.0
f_stretch_r_c3 = 0.0
f_stretch_r_c4 = 0.0
f_stretch_r_p_amp = -5.2336587774582233
f_stretch_r_p_xa = 0.061545115391222878
f_stretch_r_p_xb = 0.90340829556861646
f_stretch_r_p_w = 0.024
```
Do NOT keep the old f_stretch_r_c5..c8 lines: they would be read.

dr in Rsun at selected radii:

| r | dr | dr/r |
|---|---|---|
| R_acc | 2.477e-3 | 2.75e-4 |
| 9.5 | 2.477e-3 | 2.61e-4 |
| 9.65 | 2.63e-3 | 2.73e-4 |
| 10.0 | 3.64e-2 | 3.64e-3 |
| r_out | 7.78e-2 | 5.80e-3 |

Max dr/r is 1.23e-2, at r_in.

## 7. Ambient (code units, P_ph = 20.78, stream ram at r_out = 1.82e4, inflow 17.2 per orbit in the wedge)

| env_amb_rho / dfloor | P_amb/P_ph | cold layer top | P_amb(r_out)/ram | ambient mass / inflow per orbit |
|---|---|---|---|---|
| 1e-6 / 1e-7 (current) | **4.3e-3** | 5.4 H_p | 1.9e-6 | 1.3e-6 |
| 2e-7 / 2e-8 | 8.7e-4 | 7.1 H_p | 3.9e-7 | 2.5e-7 |
| **1e-7 / 1e-8** | 4.3e-4 | 7.7 H_p | 1.9e-7 | 1.3e-7 |
| 5e-9 / 5e-10 | 2.2e-5 | 10.7 H_p | 1e-8 | 6e-9 |

For reference:
- RY Per h5 (PASS) had P_amb/P_ph = 1.9e-5.
- RY Per a3 (FAILED: outflow at phi ~ 0) had 1.9e-3.
- Plaskett env12 had 7.5e-4.

**Decision:** keep the old input's ratio of a 1e-6 ambient to the 1e-7 floor and lower both. At the new density scale the current 1e-6 gives 4.3e-3, which is above the failed a3 level. Recommended: **env_amb_rho 1e-7 with dfloor = rho_amb 1e-8.** Mass and ram criteria pass in every case.

**Side effect:** a lower ambient raises the stream/ambient density contrast from 1e6 to 1e7. That contrast was implicated in the q10 interface failure.

## 8. r_meas

- **Using H_p at phi 0, where r_ph is largest:** max r_ph + 10 H_p = 9.4847 + 0.1302 = 9.6150, giving face 401 = **9.615674**. r_out - r_meas = 3.886 Rsun.
- **Using H_p at phi 90:** 9.5872, giving face 390 = 9.588346.

## 9. Run length

- tlim = 10 P = **4.58266**.
- bin dt = P/20 = **0.022913**.
- The input's 8.525 (= 18.6 Plaskett orbits), bin dt 0.04263 (= P/10.75) and rst dt 0.0485 are RY Per numbers, P = 0.8519.

## Other inconsistencies noticed

- The plaskett_env12 header and some comments are RY Per text: "c_ph 15.5", "564 x 4 x 2048", "rho_ph 20".
- env12 sets no problem/stream_width, so the default c_s,don/Omega applies. With t_don 33000 that is 20.963/13.71 = 1.529 Rsun, versus 0.735 here.
- The pgen's r_top search steps are log-grid, dlg = ln(x1max/x1min)/nx1. They do not match the stretched grid's faces. This is harmless: r_top need not be a face.
- I did not check the effect of the pgen's 1e-4 launch offset and of the RK4 h = 2e-6 integration. The central trajectory here (RK4 with h = 2e-5, plus a solve_ivp check to the photosphere) gives phi 3.799 deg, v_r -122.93, v_phi 55.47. That matches the earlier run.log values (-123.05, 56.08, which used c_s,don = 20.96) up to the c_s change.
