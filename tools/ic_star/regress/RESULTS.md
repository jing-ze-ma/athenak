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
