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

See `/raven/ptmp/jinma/ic_star_runs/hegiant_14p5M_global/diag_hegiant_14p5M_global.txt`
(filled in below when the run ends).
