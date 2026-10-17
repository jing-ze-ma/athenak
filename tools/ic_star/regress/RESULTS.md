# tools/ic_star regression and validation results (10-15/16, Raven login node)

Python /u/jinma/miniconda3/envs/magritte (numpy 1.23.5, scipy 1.10.1). Run outputs are in
/raven/ptmp/jinma/ic_star_runs/<name>/ (ic_, deck_, diag_ .txt/.png, prov_ .json).

## 1. Mode global, AG Car A deep70 (the E1 production IC)

`specs/agcar_A_deep70.yaml` (mltpp agcar_legacy, pp_alpha 0.970755, KW, alpha 1.5;
legacy_deep two-stage) vs `/raven/ptmp/jinma/agcmodel_1012/ic/ic_agcar_A_ge_deep70.txt`:
- 8971 = 8971 rows, r column bitwise identical.
- Max relative difference per column:

  | column | max rel diff |
  |---|---|
  | rho | 3.65e-8 |
  | eint | 3.75e-8 |
  | F_r | 3.2e-9 |
  | E | 2.6e-8 |
  | T | 6.5e-9 |
  | F_MLT/F | 1.0e-7 absolute (1 unit in the last printed digit) |

  The 2869 deep rows carry the 3.6e-8; the 6102 rows from 70 Rsun upward stay within 5.4e-9.
- Origin is platform round-off, not the method:
  - the ORIGINAL Viper script make_ic_agcar_ge.py, rerun on Raven in the same Python env as
    the tool, gives a file whose 6102 rows are **text-identical** to the tool's (0 difference);
  - both differ from the 10-08 Viper file by <= 5.4e-9 (scipy/numpy versions);
  - the deep stage restarts from the printed first row, which grows that to 3.6e-8.
- Rules on the E1 values:
  - t_Fe = 2.0149e6 s, so ramp 4.0297e6 / 2.0149e7 (E1 4.032e6 / 2.016e7, -0.06 %);
  - implicit_flux_x1min = 8.142577063662e12 (E1 8.142577965749e12 = L/(4 pi x1min^2), 1.1e-7);
  - floor 1e-16, dmax 2e-16, dfloor 1e-20, sponge 1e-3, tfloor 3000, logt_min 3.4, pfloor 1e-12,
    vceil 1e8: all equal to E1.
- Physicality flags on E1 itself:
  - m_env(83 Rsun..R_ph) = 0.197 M (10.85 Msun of 55): the point-mass gravity flag;
  - a density inversion below R_ph (inflated envelope);
  - wall buffer 3.44 H_p (< 4).

## 2. hg145 ic2 reproduction (scheme mltpp_inward, make_ic_agcar_mltpp.py --mode inward)

`specs/hegiant_14p5M_hg145compat.yaml` vs `/raven/ptmp/jinma/hegiant145_1015/ic2/ic_hg145_mltpp.txt`:
- **BIT FOR BIT**: 14354/14354 rows text-identical.
- This needs the plain extrapolating spline (no clamping), as the original used.
- Table check on that column, the reason for the new hard check:
  - 10017 rows lie outside the EOS dump dd45fb60 (log rho -21..-1 with margin), r 7.0-51.9
    Rsun, log rho up to +0.68;
  - 4529 rows lie outside the opacity grid (log rho > 0), r 7.0-17.3 Rsun.
- The hg145 production candidate ic3/ic_hg145_mesaP.txt (log rho -16.0..-3.79, log T
  3.60..6.15) is inside all tables. The superseded ic2 column was the only one outside.

## 3. Tables

| check | result |
|---|---|
| ext2 X 0.36 Z 0.02 vs E1 rc_1009 rosseland_ext2 b28e97c5 | every value identical |
| ext2 X 0.36 Z 0.02 vs E1 rc_1009 planck_ext2 dc482452 | every value identical |
| E1 planck_ext2orion 493f9d01 | differs in 4462/91357 values (the Orion low-rho stitch, not part of the recipe) |
| ext2 X 0.7 Z 0.008 vs bsg_newsetup_1015 ext2nolo, Rosseland and Planck | identical |
| ext2 14p5M (saha_mixture = hg145 mixture) vs hg145 tables/ext2 Planck | identical |
| same, Rosseland, with the TOPS-surface mixture for Saha | up to 0.40 dex at log rho < -14 only (fixed by `saha_mixture`) |
| EOS dump recipe | byte-identical to AG Car b5e87701, BSG 193fb40b, hg145 dd45fb60 |

## 4. Mode mesa

**14p5M** (`specs/hegiant_14p5M.yaml`, x1min 8, bottom 6 Rsun, ext2 at the hg145 composition):
- R_ph 1.0381 R_MESA; Teff 0.981;
- m(x1min) -6.8e-6 rel; m_env 1.042 x MESA;
- F_rad > F nowhere;
- CZ: ours 6.0-113.9 Rsun vs MESA 18.6-109.4 Rsun. Below 18.6 Rsun our CZ has F_MLT/F only
  0.05-0.16 (our TOPS kappa is above MESA's), so F_MLT/F(x1min) = 0.056;
- T/T_MESA 1.000-1.044 and rho/rho_MESA 0.999-1.073 up to 87 Rsun.
- vs hg145 ic3 (same method, the worker's own builder):
  - rho, T and F_r agree to 0.1-0.7 % at 6.5-100 Rsun;
  - F_MLT/F at the wall: ours 0.056 (F = L_ph/4 pi r^2), theirs 0.061 (F = L_MESA(r)/4 pi r^2).
  - The difference comes from the anchor (they start at 8 Rsun, we at 6 Rsun), the L
    definition and the integrator. Above 105 Rsun the two atmospheres differ by design: theirs
    is the RE atmosphere from MESA's surface, ours is Hopf + iso. In the run, sph1d/re_atm
    replaces both.

**BSG** (`specs/bsg.yaml`, cut at 48.32 Rsun):
- the profile top is reached at 47.78 Rsun, then the outward MLT extension is shot:
  eps = +0.3698 (a +37 % density jump at the profile top), and R_ph = 53.0297 (target 53.0298);
- m(x1min) +4e-6; m_env 0.985 x MESA;
- F_rad > F on 125 rows at 38.7-39.7 Rsun, max excess 0.7 %;
- CZs:
  - ours (our kappa on MESA's T(P)): 26.2-38.3 (F_MLT/F <= 0.13), 40.1-45.0 (<= 0.05),
    46.7-47.6 Rsun;
  - MESA: radiative everywhere;
- T/T_MESA 1.000 -> 0.909 at 48.3 Rsun: a radius offset growing outward (our EOS rho is
  0.98-1.0 x MESA's at equal P, T);
- the r_in rule gives x1min 27.5 Rsun (4.5 H_p below the Fe peak at 41.1 Rsun; MESA has no
  CZ there);
- floor rule 1e-17 (rho_ph 4.16e-12);
- the iso cap above R_ph adds only tau 0.017. In the run, sph1d/re_atm rebuilds this layer.

**15M** (`specs/hegiant_15M.yaml`, rule: MESA Fe-CZ base 11.38 Rsun, x1min 6.26 Rsun):
- R_ph 1.0146 R_MESA; m_env 1.007;
- F_rad > F nowhere;
- our CZ reaches the column bottom (F_MLT/F 0.08 there, interim TOPS);
- T/T_MESA <= 1.044 and rho <= 1.12 below 0.95 R_ph.

**M2pt754** (reference): MESA is not in thermal equilibrium (L(3 Rsun) 418 L_ph). Under
L_ph, F_rad > F on 8558 rows (2.85-15.8 Rsun, up to +127 %), so the scaffold is clipped to 0
there. R_ph 1.055; m_env 1.064.

## 5. Mode global validation on 14p5M (MESA r15140 gated MLT++, alpha 2, Cox, alt H_p, ext2)

`specs/hegiant_14p5M_global.yaml` (validation only; x1min 8, bottom 6 Rsun; L, M, R_ph = MESA):
- MESA gate fixed point gradT_excess_alpha* = 0.528. The rule(alpha) jumps between 0 and 1 near
  0.53: beta_min 0.31-0.42 at 19 Rsun, lambda 0.70-0.87.
- Main column (dlnr 2e-4, rtol 1e-12):
  - CZ base (Schwarzschild) 19.35 Rsun, MESA 18.6-19.4: it coincides;
  - wall buffer 4.32 H_p, radiative;
  - **m_env(8 Rsun..R_ph) = 0.095 Msun = 2.8 x MESA's 0.034**;
  - **at 25 Rsun T 3.33e5 K = 1.31 x MESA, rho 8.70e-7 = 3.86 x MESA**;
  - T/T_MESA 1.17-1.36 and rho/rho_MESA 1.5-4.3 over 0.05-0.8 R_ph;
  - a density inversion at 0.975-0.995 R_ph;
  - classified INFLATED (an inversion) but with m_env/M 0.035, an intermediate state between
    the two branches.
- Branch scan (coarse, MLT++ 'on' = the same alpha* 0.528):

  | alpha_MLT | MLT++ | T(8 Rsun) | rho(8 Rsun) | m_env/M | branch |
  |---|---|---|---|---|---|
  | 1.0 | on | 6.12e5 | 8.2e-6 | 3.2e-3 | inflated |
  | 1.0 | off | 4.78e5 | 2.6e-6 | 1.8e-3 | inflated |
  | 1.5 | on | 3.84e6 | 4.3 | 1.2e4 | DENSE ADIABAT |
  | 1.5 | off | 5.66e5 | 5.8e-6 | 2.3e-3 | inflated |
  | 2.0 | on | 3.84e6 | 4.4 | 1.25e4 | DENSE ADIABAT |
  | 2.0 | off | 6.21e5 | 8.7e-6 | 2.9e-3 | inflated |

  MESA itself has T(8) 9.0e5 K, rho 4.0e-5 and m_env/M 0.0124.
- Verdict: **global mode does NOT reproduce this MESA envelope.** Without MLT++ every alpha
  gives the inflated branch (m_env 0.1-0.2 x MESA). With MESA-gated MLT++:
  - alpha 1.5-2 at coarse tolerance gives the dense adiabat;
  - alpha 2 at full tolerance gives an intermediate envelope 2.8x MESA's mass.
  So in this super-Eddington envelope (Gamma_L up to 5.6) the inward integration sits on a
  separatrix, and the branch depends even on the ODE tolerance.
- Mode mesa is the IC for stars with MESA models. Global mode stays for stars without one
  (AG Car), always with its flags and the scan.

## 6. Photosphere seam removal (10-16)

Binary for the seam passes: /raven/ptmp/jinma/icstd_1016/bin/athena_he_cpu_b16b66db (CPU build of
fin51-r-fix b16b66db, md5 13782668). Gates: /raven/ptmp/jinma/icstd_1016/runs/<name>, production
A100 binary b16b66db (md5 042157c1), 1 A100, 16x16 lateral, seed 0, 3 d, `regress/gate_atm.py`.

**Regressions after the change** (taper and seam off in the legacy specs):
- AG Car A deep70: unchanged (max rel 3.06e-6 = F_MLT/F last digit; rho 3.65e-8);
- hg145 ic2: still bit for bit (14354/14354).

**AG Car A, if regenerated with the defaults** (`specs/agcar_A_deep70_seam.yaml`, seam on the E1 deck):
- the taper removes up to F_MLT/F 0.32 at 0.978-0.995 R_ph;
- the seam fixed point converges in 7 passes, oscillating with ratio -0.48 (+7.1, -2.9, +1.4,
  -0.66, +0.32, -0.15, +0.07 cells) to **R_ph 2.71877e13 = 390.8 Rsun (+0.64 % vs 388.3)**;
- the interior is identical below 370 Rsun; rho at 388-400 Rsun is 1.3-16x the old isothermal
  cap (the RE + HSE atmosphere).

**BSG** (`specs/bsg.yaml`, cut 48.32 Rsun, two-sided fit, taper 30 -> 3, seam on the g480 deck):
- Two-sided fit, free parameter R_ph:
  - P is continuous at the fit point (dlnP 3e-6);
  - the T-continuity root lies just beyond the profile top, so the fit takes min |dlnT| over the
    +-4 % scan;
  - result: **R_ph 3.5786e12 (0.970 x Ma's 53.03 Rsun)**, fit point r* 42.31 Rsun;
  - **continuity residuals dlnT -5.3e-3, dln rho +3.8e-2** (the previous eps-shoot jump was
    +37 %);
  - no root exists with |dR_ph| <= 1 %: the T(P) match needs R_ph ~0.970-0.975.
- Seam: pass 0 gives 3.5786e12 -> 3.5877e12 (+3.9 cells); pass 1 gives 3.5877e12 -> 3.5877e12
  (0 cells). **R_ph = 3.5877e12 = 51.57 Rsun.**
  - In the gate runs (full in-code chain) the startup sph1d gives 3.58803e12 on g480 (+0.14
    cell) and 3.589e12 on the P grid (+0.12 cell), with max |E/E_col - 1| 2e-5 (g480) and 3e-4.
  - The worker's old tool column had a 1.5 Rsun mismatch (53.03 vs 51.5).
- Gate (a), band 0.99-1.1 R_ph, max |<v_r>| in km/s / max |rho/rho_0 - 1|:

  | grid | 0.45 d | 0.75 d | 1.05 d | 1.65 d | 2.25 d | 2.85 d | verdict |
  |---|---|---|---|---|---|---|---|
  | g480 | 0.61/0.18 | 1.13/0.27 | 4.38/0.26 | 14.0/1.0 | 10.2/2.1 | 15.2/129 | FAIL (growth) |
  | P (256) | 0.38/0.41 | 1.28/0.51 | 2.23/0.57 | 7.3/0.94 | 8.7/2.6 | 17.8/201 | FAIL (growth) |

  - The interior 0.5-0.9 R_ph stays at <= 0.4 km/s.
  - The 0.9-0.99 R_ph band grows to 2-3.7 km/s.
  - The growth matches the worker's d0/d1 (4-8 km/s at 1-1.25 d, 40-50 by 2-3 d) at about half
    the amplitude.
  - **The seam and the density jump are removed, and the oscillation survives:** not a setup
    defect of the IC, most likely a physical (or closure) pulsation of the near-Eddington outer
    envelope.

**He giant 14p5M** (`specs/hegiant_14p5M.yaml`, mode mesa):
- The code's own startup (sph1d + re_atm + HSE) does NOT converge in the seam passes for any
  column with a radiative photosphere. Configurations tried:
  1. default taper 30 -> 3;
  2. no taper (the column's F_MLT/F is ~0 by tau 3);
  3. scaffold extended from tau 30 to 2/3;
  4. MESA T(P) to MESA's surface with the taper at tau 1 -> 0.1.

  Every one ran 400 outer HSE iterations with Gamma > 0.9 clamped in 7.8e3-1.3e4 fine cells,
  and R_ph moved +1..+4 % per pass. That is about +31-35 cells in one pass; with the default taper
  R_ph grew another 3.3 % in pass 1, so the iteration diverges.
- hg145 ic3 'converged' only because its tau jumps from 0.76 to > 3000 at 112.6 Rsun: sph1d then
  solved only the thin layer above. That jump is the dense shell seen in the hg145 gate.
- The seam is a code-side limit for this super-Eddington photosphere: the He/H recombination
  opacity cliff, kappa 0.01 -> 5 between 7.5 and 20 kK, Gamma_L up to 6.
- Gate column: specs/hegiant_14p5M.yaml with `top: surface`, taper tau 1 -> 0.1 (extend), seam off.
- In the gate startup the sph1d HSE converged, but to a collapsed atmosphere:
  - R_ph 7.964e12 -> 7.8635e12 (-1.26 %, ~11 cells);
  - rho above x 0.006;
  - photosphere at the r_b of the RE zone.
- Gate (a) on the deck grid 704, band 0.99-1.1 R_ph:
  - max |<v_r>| 1.42 km/s (0.4 d), DECAYING to 0.7-1.2 km/s by 2.85 d (no growth);
  - max |drho| 2.3: the collapsed thin atmosphere refills over the first ~1 d, then stays;
  - interior quiet (<= 0.02 km/s);
  - FAIL on the 10 % rho criterion, no oscillation.
  - Compare the worker gate: +16 km/s and rho x548 by 3.2 d.
- Gate (a), 2x radial (1408, same startup R_ph 7.8636e12), band 0.99-1.1 R_ph:
  - max |<v_r>| 5.9 km/s at 1.0 d (1.03-1.1 R_ph), decaying to 1.8-2.0 km/s by 2.85 d (no growth);
  - max |drho| 14 (the collapsed thin atmosphere refills);
  - FAIL on the thresholds.
- Verdict for 14p5M: no growing oscillation on either grid. The failure is the startup collapse
  of the thin atmosphere (rho x0.006 by the in-code sph1d HSE), which refills. That is a setup
  defect of the in-code chain for this star, which the IC tool alone cannot remove.

## 7. BSG physics vs closure, cheap diagnostics (10-16; regress/bsg_puls.py, gate_atm.py)

**(a) 10-d gate** on the P grid (runs/bsg_P10, same IC/deck as bsg_P, tlim 8.64e5). Band max
|<v_r>| [km/s] / max |rho/rho0 - 1|:

| t [d] | 0.9-0.99 R_ph | 0.99-1.01 | 1.01-1.03 | 1.03-1.1 |
|---|---|---|---|---|
| 0.85 | 0.52/0.02 | 0.37/0.04 | 1.8/0.21 | 0.39/0.53 |
| 1.85 | 2.6/0.24 | 2.2/0.61 | 4.4/2.1 | 6.7/3.6 |
| 2.85 | 2.7/0.10 | 12.6/0.93 | 9.6/94 | 17.8/201 |
| 3.35 | 0.93/0.01 | 39.2/0.95 | 36.5/8.2 | 24.9/39 |
| 4.85 | 0.37/0.05 | 3.4/0.22 | 53.4/0.93 | 7.5/1.1 |
| 6.85 | 0.20/0.04 | 3.3/0.22 | 7.2/7.0 | 13.8/7.7 |
| 8.35 | 0.68/0.13 | 4.5/0.35 | 22.8/1.1 | 22.1/6.8 |
| 9.85 | 0.20/0.07 | 0.24/0.12 | 10.2/0.89 | 0.20/1.0 |

- **It saturates:**
  - growth until about 3.3-4.8 d (peak 54 km/s at 1.01-1.03 R_ph);
  - after that an irregular, shock-dominated 5-50 km/s state in 1.01-1.1 R_ph to 10 d;
  - 0.99-1.01 R_ph falls back to 0.2-6 km/s;
  - 0.9-0.99 R_ph to 0.2-0.9 km/s;
  - the interior (0.5-0.9 R_ph) stays at <= 0.3 km/s throughout.
- Periods (FFT of the shell-mean v_r, 3-10 d): 1.46 d at 51.1-52.8 Rsun (R_ph 51.57),
  1.04 d at 53.6 Rsun, 2.43 d at 50.0 Rsun. In 0-3 d it is ~0.8-1.0 d at 52.8 Rsun (zero
  crossings 0.85-1.0 d); the window is too short to resolve the rest.

**(b) Mode estimates** for the IC column (closed wall at 20 Rsun, outer bound 52.06 Rsun, Gamma_1
of the code EOS + radiation, G m(r)):
- linear adiabatic radial modes: n = 0 3.47 d, n = 1 2.13 d, n = 2 1.72 d, n = 3 1.39 d;
- sound crossing 20 Rsun..top: 3.06 d;
- **acoustic cut-off of the photosphere**: 2 pi/w_c with w_c = c_s/(2 H_P) is 1.35 d at R_ph
  (1.34-1.58 d over tau 0.1-3) with g. With g_eff = g(1 - Gamma_rad) it is 2.7 d (2.6-3.0 d);
- thermal time U/L of the layer above tau 30: 630 s (0.007 d); above tau 100: 0.036 d; above
  tau 1000: 0.57 d.
- So the observed 1.46 d at and just above R_ph sits at the photospheric acoustic cut-off
  (with g), and the 2.4 d at 50 Rsun is near the n = 1 mode and the g_eff cut-off. The early
  ~0.85-1 d is near n = 3 and the cut-off.
- The layer where the amplitude lives has a thermal time << period: strongly non-adiabatic,
  isothermal waves. It is an atmospheric acoustic resonance (waves at the cut-off steepening in
  the falling density), not a deep kappa-mechanism mode.

**(c) Gamma along the column:**
- Gamma_rad = Gamma_L, since the layer is radiative with F_MLT = 0 above 44 Rsun;
- 0.53-0.67 over 43.8-50.5 Rsun, max 0.673 at tau 35 (50.0 Rsun), falling to 0.50 at R_ph;
- the rms <v_r> peaks at 51.8-53.2 Rsun, i.e. at and above R_ph where Gamma 0.39-0.50, NOT where
  Gamma -> 1 (no layer of this IC reaches Gamma > 0.68);
- there is no He II / H convection zone in this column (F_MLT/F = 0 above 44 Rsun; MESA is
  radiative too).

**Closure test (described only; a code path is needed for a clean one).**
- Ingredients: the same P-grid gate with the thin-atmosphere closure changed: vet_col
  (production) vs an M1/Eddington-factor closure in the optically thin layer, keeping the IC chain
  (ic_rad needs vet_gd/vet_col: the worker's e1 run FATALed with vet_gd off).
- A clean test needs a key that switches only the closure above tau ~1 after the startup, e.g.
  `rad_m1/closure_thin = m1` applied for tau < 1 from cycle 1.
- Interpretation:
  - if the oscillation persists at the same 1.4-1.5 d period with an M1 thin closure, it is the
    hydrodynamic atmospheric resonance (physical, set by the cut-off frequency of this
    near-Eddington photosphere);
  - if it disappears, the VET thin-layer coupling drives it.
- Cheaper partial test with existing keys: halve and double the IC atmosphere scale height
  (he_ic_rad_taua/taub or the floor) to see whether the period tracks the cut-off (period ~ H/c_s).

### 7c. Cut-off period test: the atmosphere scale height at R_ph changed x0.7 / x1.4

**Knob.** There is no key that changes only the atmosphere, so the test used three consistent
diagnostic stars (specs/bsgG_{base,H07,H14}.yaml, mode global):
- L and composition of the BSG, R_ph fixed at the BSG fixed point (51.6 Rsun);
- M = 20 / 28.57 / 14.29 Msun, so g at R_ph = 205 / 294 / 146 cm/s^2 and the gas scale height
  at fixed T_ph goes x1 / x0.7 / x1.4;
- each passed through the seam fixed point: 2 passes each, HSE converged, 0 clamped cells.

Gates: P grid, 1 A100, seed 0, 5 d (runs/cut_bsgG_*). The base of sect. 7a (mode mesa, runs/bsg_P10)
is listed for comparison.

**Predicted cut-off** 4 pi H_P/c_s at R_ph, H_P = P/(rho g) with the total P, from the IC column.
Because the radiation pressure fraction changes too, the real H ratios are 0.58 / 1.77, not
0.7 / 1.4.

| star | g_ph | H_P [Rsun] | c_s [km/s] | cut-off (g) [d] | cut-off (g_eff) [d] | observed period [d], 0.5-5 / 2-5 d | max / rms \|<v_r>\| [km/s] (0.99-1.1 R_ph) |
|---|---|---|---|---|---|---|---|
| H07 (M 28.6) | 294 | 0.201 | 23.5 | 0.87 | 1.48 | 0.87 / 0.77 (zero crossings 0.80) | 6.1 / 1.3-2.0 |
| base (M 20, global) | 205 | 0.347 | 25.9 | 1.36 | 2.73 | 1.01 / 0.97 (0.98-1.00) | 11.5 / 4.2-4.3 |
| H14 (M 14.3) | 146 | 0.616 | 29.1 | 2.14 | 5.36 | 1.21 / 1.43 (1.17-1.27) | 28.7 / 8.9-10.7 |
| base (mode mesa, bsg_P10) | 194 | | | 1.35 | 2.72 | 1.39 / 1.27 (1.27-1.68) | 54 / 17-21 |

- Periods are FFTs of the linearly detrended shell-mean v_r (Hann window, zero-padded), median
  over the 5 radii of largest rms in 0.99-1.1 R_ph; the peak sits at 1.02-1.03 R_ph, and at
  1.05 R_ph for H14.

**Result.**
- **The period tracks the scale height monotonically:** H07 < base < H14, ratio 0.8 : 1 : 1.3-1.4.
  It does so sub-linearly: the predicted cut-off ratio is 0.64 : 1 : 1.57 (period ~ H^0.5-0.6).
- The observed period is 0.9-1.0x the g cut-off for H07 and the mesa base, 0.7x for the global
  base, 0.6-0.67x for H14. The g_eff cut-off is too long in every case.
- **The saturated amplitude grows with the scale height:** max |v| 6 -> 11.5 -> 29 km/s, rms
  ~1.6 -> 4.3 -> ~10 km/s. In H14 the oscillation sits higher (1.05 R_ph).
- The amplitude also depends on the interior: global base 11.5 km/s vs mesa base 54 km/s at
  similar g; the deeper structure differs.

**Interpretation.** The atmospheric oscillation is set by the stratification of the photosphere
and thin atmosphere, consistent with an acoustic resonance near the cut-off frequency that is
excited at the photosphere and steepens in the falling density.
- This test does NOT separate physics from closure: a thin-layer VET-coupling instability would
  also scale with H.
- The described closure-switch test (sect. 7, a key for an M1/Eddington closure at tau < 1) is
  still the discriminator.
- Caveat: g also changes Gamma (0.47 / 0.67 / 0.94 x the base's at the Fe peak), and with it the
  interior structure (m_env/M 0.08 for H07).

## 8. Exact TOPS tables and the startup gamma cap (10-16, USER: BSG R_ph 51.6 accepted, 14p5M option A)

### 8a. Tables rebuilt by the recipe; the cache key now includes the md5 of the TOPS raw files

- **14p5M**: TOPS of the exact mixture (hg145 tables/tops_hg145.dat, LANL OK on try 35 10-16
  10:25, rho 1e-18..1). Cache opac_416a61196c64, Saha mixture = the 14p5M mixture.
  Change vs the interim table (old-surface-mixture TOPS, opac_ef4dcdfaea5b), in log10 kappa:

  | region | Rosseland | Planck |
  |---|---|---|
  | log T >= 4.2, log rho >= -14 | max 0.126, mean -0.011 | max 1.12, mean -0.036 |
  | log T >= 4.2, log rho < -14 (extension) | max 0.005 | max 1.44, mean -0.10 |
  | log T < 4.2 (Ferguson, blend) | max 0.40 (mostly the old surface-mixture Saha-es extension) | max 1.16 |

- **BSG**: the TOPS low-density file (Viper retry OK on try 37, 10:18) added. Cache
  opac_5b8e2fdfefea. Change vs ext2nolo:
  - Rosseland: max 0.006 dex, only at log rho < -14;
  - Planck: max 0.43 dex at log T >= 4.2 / log rho < -14 and 0.52 dex at log T < 4.2 (the
    extension), mean +0.003.
  - The BSG IC was not regenerated: the Rosseland change is below log rho -14 only.

### 8b. Option A: problem/he_ic_rad_gamma_cap (branch ic-gammacap-1016, commit 843ce3e4, from b16b66db)

Implementation: src/pgen/he_star_m1.cpp, about 60 lines in the sph1d block.
- The key is read only when named.
- Inside the HSE outer loop, after kt and the march gravity are known, the column's F_r is capped
  at cap c g/kt per fine cell above r_b, and the excess goes into the frozen scaffold:
  F_MLT/F = 1 - F_r/F_tot, with F_tot = F_r0/(1 - F_MLT0/F) of the column, so F_r + F_MLT = L.
- The cap is re-applied from the column's own F_r each outer iteration (no compounding).
- The dump writes the capped scaffold.
- It needs he_ic_rad_hse = true and mlt_flux_frozen = true.
- Binaries: /raven/ptmp/jinma/icstd_1016/bin/athena_he_{cpu,a100}_843ce3e4.

**Key off is bitwise** (b16b66db vs 843ce3e4, 16x16 lateral, seed 0; every .bin and .hst output
and the sph1d R_ph line compared):

| deck | CPU (3 cycles) | A100 (10 cycles) |
|---|---|---|
| E1 (e1_prod_rol2u) | 9 files, 0 differ | 16 files, 0 differ |
| BSG (bsg_seam_g480) | 13 files, 0 differ | 20 files, 0 differ |
| He giant (hg145_seam) | 11 files, 0 differ | 18 files, 0 differ |

**AG Car B, on vs off** (documented only; shakeB deck + ic_agcar_B_ge.txt, 0 cycles, CPU,
cap 0.9):

| | off | on |
|---|---|---|
| R_ph | 7.04718e12 -> 7.12966e12 | 7.04718e12 -> 7.12992e12 (+2.6e8 cm) |
| capped cells | - | 542 fine cells, max F_MLT/F added 0.061 |
| Gamma > 0.9 clamped cells | 538 | 60 (equality at round-off) |
| HSE | converged (7.0e-7) | converged (9.2e-7) |

**14p5M** (mode mesa, exact TOPS tables, taper 30 -> 3; CPU 0-cycle passes on the hg145 deck):

| cap | R_ph 7.96114e12 -> | HSE last max \|dln rho\| | F_r capped (last outer) |
|---|---|---|---|
| 0.9 | 8.10596e12 (+17 cells) | 2.0, not converged | 6193 cells, max +0.88 F_MLT/F |
| 0.7 | 8.00198e12 | 1.85, not converged | 1975 cells |
| 0.5 | 8.11212e12 | 1.59, not converged | 6475 cells |
| 0.3 | 7.96504e12 (+0.5 cell) | 3.8, not converged | 0 cells |

- Seam iteration at cap 0.3: pass 0 gives 7.961 -> 7.965 (not converged); pass 1 gives
  7.965 -> 7.814 (HSE converged, but the collapsed branch, -22 cells); pass 2 gives
  7.814 -> 8.024 (not converged, +26 cells). **It flips between a collapsed and an inflated
  branch; there is no fixed point.**
- The same column and key in the A100 gate startup gives 7.961 -> 8.094 (inflated), while the CPU
  pass gives 7.965. The branch is picked by round-off (device vs host opacity evaluation): the
  startup problem is ill-conditioned for this photosphere.
- With cap 0.3 no cell is capped in the final outer iteration, yet the loop still oscillates. The
  Gamma limit is not what drives it. The driver is the Picard alternation (E solve at fixed rho
  <-> HSE at fixed T) on the steep kappa(T) of the He/H recombination cliff.
- **Option A therefore does NOT give a converged, reproducible startup for 14p5M.** The next step
  is option B of README sect. 2c: a damped or coupled E-HSE update, or a continuation that detects
  the fold. It needs the USER.

### 8c. 14p5M gate (a) with the cap (cap 0.3, exact TOPS tables, column = the tool's mode-mesa column after the taper)

Setup: A100 843ce3e4, 16x16 lateral, seed 0, 3 d. Deck hg145_cap2 = the hg145 deck + tool block +
he_ic_rad_gamma_cap 0.3 + the new tables + implicit_flux_x1min 3.6060e13 (the new tables change
F_r(8 Rsun) by +4.6 %).

Startup: R_ph 7.961e12 -> 8.094e12 (704) and 8.084e12 (1408): the inflated branch, HSE not
converged.

Band 0.99-1.1 R_ph, max |<v_r>| [km/s] / max |rho/rho0 - 1|:

| grid | 0.45 d | 1.05 d | 1.65 d | 2.25 d | 2.85 d | verdict |
|---|---|---|---|---|---|---|
| 704 | 2.8/0.56 | 3.8/1.4 | 7.0/5.6 | 12.0/25 | 19.3/182 | FAIL (growth) |
| 1408 | 5.5/1.3 | 6.9/5.6 | 11.1/24 | 20.3/161 | 27.3/3456 | FAIL (growth) |

- The interior below 0.9 R_ph stays <= 0.3 km/s.
- This is worse than the uncapped run of sect. 6, whose startup collapsed and then refilled
  without growth.
- **The cap neither converges the startup nor passes the gate. No 14p5M production IC.**

## 9. MESA-anchored 14p5M startup (10-17, branch ic-mesascaf-1017; README sect. 2d)

Work dir: /raven/ptmp/jinma/hg145_mscaf_1017. Scripts are in scripts/: d1/d2/d3 diagnostics,
mkcol, runpass/seamloop, mkkapm, cmpdump, diffchk, discflux, prof, mkgate. CPU passes use
/raven/ptmp/jinma/icstd_1016/bin/athena_he_cpu_b16b66db; gates use the A100 b16b66db build
(/raven/ptmp/jinma/overlap_1015/bin, md5 042157c1).

### 9a. MESA's own outer envelope (opacity-free; scripts/d1.py)
- kappa_Edd = 4 pi c G M/L = 0.969 cm^2/g. MESA's kappa: 0.003 at the photosphere, 0.12 at
  tau 10, 0.77 at tau 27, a peak of 4.9 at 18.7 kK (tau 1600), 3.4 at 42 kK (tau 8000).
- So Gamma_L = 1-5 for tau 30-1.5e4, and MESA's F_MLT/F there is 0.88-0.99.
- Gamma_rad of MESA's actual gradient, 4 (P_rad/P) nabla, is at most **1.009** (108.4 Rsun, tau
  2900): a density inversion, 1.3e-8 -> 2.4e-8 g/cc outward over 108.6-109.3 Rsun.
- MESA is radiative (F_MLT/F < 0.05) only at tau < 3.4. Gamma_L < 0.12 at tau < 10.
- Our kappa over MESA's along MESA's track: 1.0-1.1 at 35-60 kK; 1.4-1.5 at 24-33 kK; 1.1-1.3
  at the 18 kK peak; **2.0-2.5 at 7.6-11 kK** (Ferguson low-T table, GS98-scaled He-rich).

### 9b. Does a static solution exist?
- **Fully radiative: no.** In a plane-parallel grey RE + HSE march (scripts/d2.py, Teff 7504 K),
  Gamma = 1.2-3.4 at tau 20-1300 with our kappa. The gas pressure must then fall inward
  (rho 1.4e-8 -> 1.2e-9). Deeper, where Gamma_L is 1.5-3.5 down to 25 Rsun, P_gas -> 0.
  The photosphere itself, tau < 10 with Gamma < 0.1, is a normal static radiative layer.
- **With MESA's convective flux: yes**, and this is MESA's model. In our microphysics the
  scaffold must be MESA's gradient re-expressed in OUR kappa:
  - F_MLT/F = 1 - nabla_M/nabla_rad(ours), i.e. the mode-mesa column;
  - its F_r is the diffusion flux of the column's own T (checked to 1e-3, scripts/diffchk.py);
  - its Gamma_rad is max 1.0098 (12 rows > 0.9, all at tau > 3500).
- MESA's OWN F_MLT/F at the same P_tot with our kappa gives **Gamma_rad 1.49**: 29 rows > 1 at
  26-32 kK, because our kappa is 1.4-1.5x there. That split has no static solution (the tool's
  diag now prints both).
- So the right IC is static (not a wind), with the convective scaffold carrying the flux down to
  tau ~3.

### 9c. Why the code's startup failed: the RE re-solve on the cliff, for either opacity
CPU 0-cycle passes on the hg145 deck (passes/):

| column / keys | result |
|---|---|
| raw column (no taper), default taua/taub 30/3000 (p1L_raw) | 400 outer, NOT converged, Gamma clamp in 9486 cells, R_ph +25 cells (inflated) |
| raw, taua/taub 3/30 (p2L_t3_30) | 400 outer, NOT converged (max dln rho 0.79) |
| MESA-like kappa (diag table: ours x kappa_M/kappa_ours(T) on MESA's track), default 30/3000 (p3L_kapM) | converged (7e-7), but to the COLLAPSED branch: tau(r_b) 3000 -> 0.2, rho x0.01 at R_ph, R_ph -14 cells |
| raw, taua/taub 1/3 (p4L_t1_3) | converged, 9 outer, 6.4e-7, no clamp, R_ph -0.58 cell |

- After ONE outer iteration (hse_maxit 1, p3L_kapM_o1) the VET solution at the column's own
  rho is already 25 % cooler than T_col at tau_col 300-700. The inner VET loop converged in 25
  passes.
- Cause: with a frozen F_r(r), T^4 ~ int kappa(T) rho F_r dr. With dln kappa/dln T ~13 in the
  cliff, the Picard iteration has a cool (collapsed) and a hot (inflated) attracting branch, and
  the column is the unstable one between them.
- The column IS a solution (diffusion-consistent F_r), but the RE problem marched from the top
  is so ill-conditioned through the cliff that a 7 % T offset at tau 3 picks a branch.
- This is independent of the opacity table: MESA-like kappa collapses as well.
- The 2b taper (F_MLT -> 0 by tau 3) adds an inconsistent F_r on top (RESULTS 6). The 8b cap
  addresses Gamma, not this feedback.

### 9d. The method and its fixed point (tool, no code change)
Spec `specs/hegiant_14p5M.yaml`:
- `taper: {enable: false}`, `startup: {re_zone: hopf}`;
- he_ic_rad_taua 1, he_ic_rad_taub 3;
- seam on the hg145 deck with he_nfine 262144. The old seam passes ran the deck's 2^20 while
  the gate decks had 262144.

Run (sbatch small, 381 s): /raven/ptmp/jinma/ic_star_runs/hegiant_14p5M_anchor, IC
`ic_hegiant_14p5M.txt` md5 2863f18b.
- seam pass 0: R_ph 7.96114e12 -> 7.95671e12 (-0.579 cell), HSE 7.6e-7;
- seam pass 1: 7.95671e12 -> 7.95704e12 (**+0.043 cell**), **HSE 5.84e-7**, 0 Gamma-clamped
  cells. **FIXED POINT** R_ph = 7.95704e12 = 114.38 Rsun (1.041 R_MESA).
- L consistency 4 pi r^2 (F_r + F_MLT)/L - 1 = 4.3e-11. Before the fix, the dump interpolation
  gave 0.9 % at the F_MLT jump at the join; each pass now sets F_r = F (1 - F_MLT/F).
- implicit_flux_x1min 3.606021e13; F_MLT/F(x1min) 0.012; m_env 1.042 x MESA;
  T/T_MESA <= 1.04 up to 87 Rsun.
- **Reproducible:**
  - a second build on another node gives the same R_ph lines (7.95671 / 7.95704);
  - the A100 b16b66db run of the final pass deck (runs/a100pass1g) gives R_ph 7.95671e12 ->
    7.95704e12, identical, with the dump within 3.3e-10 (print precision) of the CPU dump in
    every column;
  - the gate startups (A100, atm_force module, 704/1408 meshes) give 7.95721e12 / 7.95717e12
    (+0.02 cell), with HSE 5.84e-7 / 5.87e-7.
- MESA-like-kappa diagnostic build (hegiant_14p5M_kapMdiag, same method): fixed point in 2
  passes (-0.012 cell). With MESA's kappa the column's photosphere is self-consistent:
  tau_R(r_j) 2.67 from the mesh top, vs **5.45 with our kappa**. Our 2.4x photospheric kappa puts
  2.4 extra tau into the column's isothermal cap.
- Default off is bitwise for the tool:
  - AG Car A deep70 IC data byte-identical to the 10-16 file;
  - BSG (mode mesa, taper + seam on) built with 3c3c3107 and with this branch: IC data and deck
    identical (the seam passes included);
  - only the new diag Gamma_rad line is added.
  No he_star_m1 change, so E1/BSG/HG decks run the same binary.

### 9e. Gate (a'): FAIL (regress/gate_aprime.py; runs/hg704, hg1408; 9 d, seed 0, 16x16)
Deck /raven/ptmp/jinma/hg145_mscaf_1017/decks/hg145_anchor.athinput (md5 915c1d80 = hg145e1 + tool block):
- [1] resolution 2.25 cells/H_rho at 1.027 R_ph on 704 (pair 4.74): FAIL;
- [2] rms 8.47 vs 8.51 km/s, PASS; period 4.45 vs 4.47 d, PASS; profile L2 17 %, FAIL;
- [3] no growth, PASS; but rms 8.5 km/s > c_s,iso 4.0, FAIL;
- [4] FAIL on both runs:
  - 0.5-0.9 R_ph 4.0/4.3 km/s;
  - 0.9-0.99 R_ph 8.1/7.7 km/s;
  - R_ph +5.7/+4.6 %;
  - L_top/L_in 0.46-1.08;
  - dM/M 1e-8 and Mdot 1.6e-9 Msun/yr pass;
- [5] PASS: NC 0, 0 retries, rc 0.

What happens (scripts/prof.py):
- the <rad_m1> ic_rad relaxation already drops L_top to **0.595 L_IC (704)**, 0.769 (1408);
- 4x4 lateral resolution runs (2 d): 0.595 (704), 0.774 (2816), 0.895 (5632);
- MESA-like kappa: 0.906 (704).
- The deficit sits at the join (113.7-114.1 Rsun, the tau 3-5 layer). There the frozen
  scaffold falls from F_MLT/F 0.98 to 0 within 0.3 Rsun (3 cells at 704), and the column's E at
  the join is about 0.6 of what RE would put under 5.4 tau of our kappa.
- Within 0.25 d the cliff top (113.2-113.7 Rsun, T 16-19 kK, F_r 1-5 % F) heats to 28-31 kK by
  1-2 d. Then the 0.9-0.99 R_ph band contracts at 2-5 km/s and the atmosphere above swells
  (rho x1e4-1e5 at 1.01-1.1 R_ph). Bounded, with a ~4.5 d period.
- The same happens with MESA-like kappa (kapM704: 0.9-0.99 R_ph 1.2 km/s at 1 d, 5 km/s at
  2 d; atmosphere x1e4) and at 2816 (0.9-0.99 R_ph 2.2-2.9 km/s by 0.45-0.75 d, vs 2.3-4.4 on
  704).
- So neither the startup nor the photospheric opacity nor the radial resolution (up to 4x) is
  the remaining blocker.
- What remains is the run itself: the frozen, r-fixed scaffold hands 0.5-0.98 L to radiation
  across the thin cliff top, which has a tiny heat capacity, kappa ~T^13, and MESA's nabla up to
  2.9 (superadiabatic) with a Gamma_rad ~ 1 density inversion. The static MESA envelope is not
  held by frozen-scaffold radiation hydro.
- **No 14p5M production IC.**

## 10. 14p5M convective-flux hand-over (10-18, USER option 1; code branch mltp-1017 @ d2239c5e off fin51-r-fix bc5d24c1)

Arms (704, 9 d, seed 0, 16x16, the sect. 9 IC/deck). he_base_heat is removed from the deck for the closures that
need mlt_closure != frozen (it acts only once the ramp starts, at 14.8 d, i.e. after the run).
- **(a) scaffold frozen in pressure**: new opt-in `problem/mlt_closure = pressure` (he_star_m1, mltp-1017).
  - F_sub = fm_IC(ln P_face) F_tot,IC(face), with P_face the shell-mean P_gas + E/3;
  - relaxed over mlt_relax_time (1500 s), every cycle;
  - the adaptive closure's restart state is reused.
- **(b) live/relaxing closure**: the existing `mlt_closure = adaptive`, the shell-mean deficit
  max(0, L - L_rad - L_conv)/A relaxed over 1500 s, masked by fm_IC/0.01.
- **(c) the hand-over at the join** is what (b) adjusts itself. No separate taper was built: a taper over a few H_p
  would leave MESA's T(P) inconsistent with F_r, i.e. the ill-conditioned cliff of sect. 9c.

Single-run (a') items. Columns: startup L_top/L_IC; band 0.99-1.1 saturated rms / max |v| [km/s];
0.5-0.9 R_ph max; 0.9-0.99 R_ph saturated max; R_ph drift; L_top/L_in after 2 d:

| arm | L_top(0) | band rms / max | 0.5-0.9 | 0.9-0.99 | R_ph | L_top/L_in |
|---|---|---|---|---|---|---|
| frozen-r, our kappa (hg704) | 0.60 | 8.5 / 24 | 4.0 | 8.1 | +5.7 % | 0.46-1.08 |
| adaptive, our kappa (adaptg704) | 0.60 (1.00 by 0.35 d) | 6.8 / 21 | 3.0 | 5.0 | +2.3 % | 0.36-1.42 |
| pressure, our kappa (pres704) | 0.60 | 8.0 / 20 | 1.7 | 5.4 | +3.2 % | 0.76-1.14 |
| frozen-r, MESA-like kappa (kapM704) | 0.91 | 7.2 / 24 | 3.5 | 5.1 | +2.8 % | 0.57-1.28 |
| pressure, MESA-like kappa (presKM704) | 0.91 | 4.3 / 19 (grows) | 0.6 | 4.0 | +1.1 % | 0.61-1.17 |
| **adaptive, MESA-like kappa (adaptkM704)** | 0.91 | **1.5 / 7.5** | **0.39** | **0.42** | **+0.06 %** | **0.97-1.04** |

- With our kappa, adaptive keeps the star quiet for ~3 d (0.9-0.99 R_ph <= 0.8 km/s). After that the
  photosphere layer restructures: the E deficit at the join from our 2-2.5x low-T kappa.
- The pressure-frozen scaffold alone does not fix it.
- **Adaptive + MESA-like kappa removes the failure mode:** no cliff-top heating, no envelope contraction, star items 3/4
  pass on 704.
- Gate (a') pair adaptkM704/adaptkM1408: FAIL items 1, 2, 3, 4.
  - 704 has 2.46 cells/H; the pair has 5.07.
  - The atmosphere is NOT resolution-converged: band rms 1.5 (704) vs 4.65 km/s (1408), periods 4.5 vs 3.1 d.
  - The 1408 run has rms 4.65 > c_s,iso 4.0 and L_top/L_in down to 0.91.
  - 1408 star: 0.5-0.9 R_ph 0.60, 0.9-0.99 R_ph 0.79 km/s, R_ph 0.43 %.
  - A deck that passes item 1 is 1408-class, with a pair of ~2816 (a restart chain, ~5 h A100).
- The MESA-like table is a diagnostic only (ours x kappa_M/kappa_ours(T) along MESA's track).
- Default off (mlt_closure frozen) is bitwise, b16b66db vs mltp-1017:
  - CPU 3 cycles: E1 9, BSG 13, HG 11 files, 0 differ;
  - A100 10 cycles: 16 / 20 / 18 files, 0 differ;
  - the sph1d R_ph lines are identical in all six.
- **Not production-ready.** Needed:
  1. a physical He-rich low-T opacity (7.5-11 kK) consistent with MESA's kap set-up;
  2. adaptive + that table, screened with 1-2 d arms (the failures show by 0.25-1 d);
  3. a final (a') pair 1408/2816.
