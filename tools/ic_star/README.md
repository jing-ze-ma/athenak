# tools/ic_star: the standard initial condition for he_star_m1 envelope runs

ONE pipeline and ONE in-code IC chain for every stellar-envelope run (BSG, He giants, AG Car /
LBVs, future stars). There are two input modes:

| mode | when | structure | microphysics |
|---|---|---|---|
| **mesa** (default when a MESA profile exists) | BSG, He giants | MESA's T(P), i.e. MESA's actual nabla, re-integrated in HSE outward from the column bottom | OUR EOS, kappa_R and gravity (G m(r)) |
| **global** (no MESA profile) | AG Car / LBVs | the AG Car surface-inward MLT/MLT++ column from L, M, Teff/R_ph | OUR EOS and kappa_R, point mass |

```
star spec (YAML/JSON) --make_star_ic.py--> tables (cached by composition hash; hard range check)
                                         -> IC column  ic_<name>.txt   (+ provenance header)
                                         -> deck block deck_<name>.athinput
                                         -> diagnostics diag_<name>.txt / .png, prov_<name>.json
```

```bash
cd tools/ic_star
export OPENBLAS_NUM_THREADS=1        # login nodes: OpenBLAS thread limit
PY=/u/jinma/miniconda3/envs/magritte/bin/python3      # numpy 1.23, scipy 1.10, PyYAML, mpl
$PY make_star_ic.py specs/bsg.yaml --out /raven/ptmp/jinma/ic_star_runs/bsg
bash regress/run_regress.sh          # AG Car A deep70 (global) + hg145 ic2 reproduction
```

A run takes 10-30 min on a loaded login node (LSODA with a Python right-hand side). In mode
global, the MLT++ gate and the branch scan add several coarse columns; skip the scan with
`--no-scan`.

## 1. Mode `mesa` (icstd/mesacol.py)

1. MESA profile: format `mesa` (profile*.data: logP, opacity, luminosity, mass), or `cgs6`
   (r T rho m kappa F_r; P is rebuilt from MESA's own HSE). `r_max_rsun` drops non-envelope
   rows at the top (L. Ma's BSG file: its last row at 52.49 Rsun follows a 4.2 Rsun gap, at
   T 63 kK and rho 2.6e-11).
2. Structure: T = T_MESA(P), PCHIP in ln P, i.e. MESA's nabla(P) = dlnT/dlnP.
3. HSE re-integrated OUTWARD from the column bottom r_bot with OUR EOS: rho solves
   P - aT^4/3 = P_gas,code(rho, T). Gravity is the code's: G m(r) = G m_MESA(r_bot) +
   G int 4 pi r^2 rho dr. The deck gets **he_gm = G m(x1min) of this column and
   he_gm_column = true**.
4. Stop where T(P) reaches the Hopf tau = 3 temperature (3F/4sigma (3 + q))^1/4. If the MESA
   profile ends first (BSG), continue OUTWARD with local MLT (spec mlt) and OUR kappa from
   rho_top (1 + eps). eps is SHOT so that R_ph = the target (spec R_ph or MESA's). That is the
   separatrix eigenvalue of an outward integration (Ma+2026 footnote 5); the density jump eps is
   reported.
5. Hopf layer outward to tau 2/3: OUR R_ph and rho_ph. Above: the isothermal gravity-only
   atmosphere, column floor and floor rows, as in mode global.
6. F_rad from OUR kappa: F_rad/F = nabla/nabla_rad(ours, g = G m/r^2, F = L/4 pi r^2).
   F_MLT/F = 1 - F_rad/F, **clipped to [0, 1]**. Where F_rad > F (our kappa below MESA's), the
   column has F_r = F and F_MLT = 0, so there is no negative scaffold. The rows, radii and maximum
   excess are reported: those layers will lose heat in the run.
7. Checks: m(x1min) and m_env vs MESA, R_ph and Teff vs MESA, our convective zones (our kappa
   on MESA's T(P)) vs MESA's own (F_rad,M/F from MESA's opacity, or F_r for cgs6), T/rho/m ratios,
   L_MESA(x1min)/L_ph flat.

## 2. Mode `global` (icstd/column.py; the AG Car method)

Same equations as make_ic_agcar_ge.py (AG Car A/B) and make_ic_agcar_deep.py (the E1 deep
wall):
- R_ph (tau_R 2/3) from L and Teff or given; grey Hopf T(tau);
- isothermal gravity-only atmosphere above R_ph, with rho_ph SHOT so that its tau(R_ph) = 2/3;
- Hopf layer with gas HSE and the radiation force of F, inward to tau 3;
- then inward dP_tot/dr = -rho GM/r^2 and dlnT/dlnP = nabla (local MLT + MLT++, = nabla_rad
  where radiative), with F_r = F nabla/nabla_rad, F_MLT/F = 1 - nabla/nabla_rad, L constant,
  GM a point mass;
- rows on a log grid dln r 2e-4 with x1min as an exact row; atmosphere rows every 2e-4 R_ph;
  floor rows every 50 x that.

Validation on 14p5M (regress/RESULTS.md sect. 5): even with MESA's gated MLT++, global mode
does not reproduce a super-Eddington MESA envelope. The inward integration sits on a separatrix
between an inflated branch and a dense adiabat (2.8x MESA's m_env at alpha 2; the branch flips
with alpha, MLT++ and even the ODE tolerance). Use mode mesa whenever a MESA model exists.

**Mandatory physicality checks (they flag, never accept silently):**
- m_env(x1min..R_ph) << M (FLAG at > 5 % M);
- a smooth base (second difference of ln rho in the bottom rows vs above);
- branch classification:
  - DENSE ADIABAT: convective to the base, or m_env > 0.3 M;
  - INFLATED: density inversion below R_ph;
  - normal;
- **branch scan**: alpha 1.0/1.5/2.0 x MLT++ on/off with coarse columns. Each row gives T and rho
  at x1min, F_MLT/F there, m_env/M and the branch.

### 2a. MLT and MLT++: MESA r15140, from the source

The source is MESAHub/mesa tag `0.0.15140` (= r15140). Copies are in
`/raven/ptmp/jinma/ic_star_cache/raw/mesa_src`, with md5s:
mlt_info.f90 95da7347, mlt_get_results.f90 3cfae1d9, controls.defaults ada07997,
star_utils.f90 e999399f.

- **On/off**: `okay_to_reduce_gradT_excess = .false.` (star/defaults/controls.defaults:2292).
  "This is the treatment of convection, referred to as MLT++ in Section 7.2 of Paxton et al.
  (2013)": **MLT++ is OFF by default in MESA.**
- **Global gate** `set_gradT_excess_alpha` (star/private/mlt_info.f90:1330-1489). It runs once
  per step and gives ONE alpha for the model:
  ```fortran
  beta = 1d0 ! beta = min over k of Pgas(k)/P(k)            (:1382-1390)
  beta = beta*(1d0 + s% xa(1,nz))                            (:1392)
  lambda = 0d0 ! lambda = max over k of Lrad(k)/Ledd(k)       (:1397-1405, k = 2..k_beta)
     tmp = get_Lrad_div_Ledd(s,k)
  lambda = min(1d0,lambda)
  if (lambda1 < 0) then; alpha = 1
  else if (lambda >= lambda1) then
     if (beta <= beta1) then; alpha = 1
     else if (beta < beta1 + dbeta) then; alpha = (beta1 + dbeta - beta)/dbeta
     else; alpha = 0
  else if (lambda >= lambda2) then
     beta_limit = beta2 + (lambda - lambda2)*(beta1 - beta2)/(lambda1 - lambda2)
     ... same with beta_limit ...
  else if (lambda > lambda2 - dlambda) then
     if (beta <= beta2) then; alpha = 1
     else if (beta < beta2 + dbeta) then; alpha = (lambda - (lambda2 - dlambda))/dlambda
     else; alpha = 0
  else; alpha = 0                                            (:1430-1460)
  if (s% gradT_excess_alpha < 1d-4) s% gradT_excess_alpha = 0d0      (:1480)
  if (s% gradT_excess_alpha > 0.9999d0) s% gradT_excess_alpha = 1d0  (:1481)
  ```
  - L_rad/L_Edd = -(area^2 crad (dT^4/dm)/3)/(4 pi G m) (star_utils.f90:1816-1835) is the
    radiative luminosity of the actual T gradient, = kappa F_rad/(c g).
  - The time smoothing (:1462-1477, gradT_excess_age_fraction 0.9) does not apply to a
    static IC.
- **Local application** `adjust_gradT_excess` (mlt_info.f90:1077-1157):
  ```fortran
  if (gradT_excess_alpha <= 0.0 .or. s% gradT_sub_grada(k) <= s% gradT_excess_f1) return
  if (s% lnT(k)/ln10 > s% gradT_excess_max_logT) return
  log_tau = log10(s% tau(k)); if (log_tau < s% gradT_excess_max_log_tau_full_off) return
  alfa = s% gradT_excess_f2 ! for full boost, use this fraction of gradT
  if (gradT_excess_alpha < 1) then; alfa = alfa + (1d0 - alfa)*(1d0 - gradT_excess_alpha)
  beta = 1.d0 - alfa
  s% gradT(k) = alfa*s% gradT(k) + beta*s% grada_face(k)
  ```
- **Defaults** (controls.defaults):
  - gradT_excess_f1 1d-4 and f2 1d-3 (:2304-2305);
  - lambda1 1.0, beta1 0.35 (:2342-2343);
  - lambda2 0.5, beta2 0.25 (:2355-2356);
  - dlambda 0.1, dbeta 0.1 (:2368-2369);
  - max_logT 8 (:2399);
  - min_log_tau_full_on and max_log_tau_full_off -99 (:2412-2413, inactive);
  - max_center_h1 1, min_center_he4 0 (:2379, 2389).
- **Implemented** (`mlt.mltpp: mesa_gated`): alpha* = rule(column(alpha*)), solved by bisection
  on [0, 1] with coarse columns, as the AG Car A fixed point of 10-08.
  - beta and lambda are taken over OUR column (envelope only), not the whole star.
  - (1 + xa(1, centre)) is a knob (`xa1_center`, default 0): xa(1) is the net's first isotope
    at the centre, h1 in most nets, so 0 for H-free cores.
- **MLT** (MESA defaults):
  - MLT_option 'Cox' (controls.defaults:2191; "'Cox' option assumes optically thick
    material", :2186), i.e. Cox & Giuli 14.78-14.82 (mlt_get_results.f90:1088-1385): A =
    K_conv/K_rad, a0 = 9/4, Gamma from the cubic, Zeta = Gamma^3/B^3,
    nabla = (1 - Zeta) nabla_rad + Zeta nabla_ad.
  - It has no radiative-loss / optically-thin correction. ML1/ML2/Mihalas/Henyey carry the
    omega = Lambda rho kappa factors, :1105-1170.
  - Our Kippenhahn-Weigert cubic (the AG Car builders, `mlt_option: KW`) gives the SAME nabla
    and v_MLT to 1e-10 at four AG Car points: it is the same model in other variables.
  - The real differences are two:
    - mixing_length_alpha = 2 in MESA (:2018) vs 1.5 in the AG Car builders;
    - alt_scale_height_flag = .true. (:2148): MESA's mixing length uses
      H_p = min(P/(rho g), sqrt(P/G)/rho) (mlt_get_results.f90:308-318), the AG Car builders
      P/(rho g).
  - Both are spec knobs (`alpha`, `alt_scale_height`) with MESA's defaults.
- **Legacy** (`mltpp: agcar_legacy`): a FIXED pp_alpha applied wherever nabla - nabla_ad > f1.
  pp_form 'mesa' is MESA's alfa blend (AG Car A, pp_alpha 0.970755); 'f2' is
  nabla_ad + f2 (nabla - nabla_ad), the make_ic_agcar_mltpp.py form. The AG Car A deep70
  regression uses it.

## 3. Step by step for a new star

1. Write `specs/<name>.yaml` (copy `bsg.yaml` with a MESA profile or `agcar_A_deep70.yaml`
   without). Give L, M (or GM), Teff or R_ph, X Y Z (+ GS98 or a TOPS mixture), and the profile.
2. Run `make_star_ic.py`. Tables are built once per composition. The run stops (FATAL) if any
   column row leaves the EOS dump or the opacity grid. A recipe-made EOS dump is widened
   automatically, once, to cover the column with 0.1 dex T / 0.5 dex rho margin.
3. Read `diag_<name>.txt`:
   - coverage (the opacity source of every row);
   - the wall (H_p buffer >= 4-5 below the Fe-CZ base, radiative buffer, L flat);
   - L consistency;
   - background (floor-layer extent; sponge ACTIVE: column floor < dmax);
   - regime flags;
   - mode checks (mesa: m_env/R_ph/T_ph/CZs vs MESA, F_rad > F rows; global: physicality
     flags and the branch scan).
   Do NOT patch the method per star: one method for all.
4. Paste `deck_<name>.athinput` into a copy of the current production deck (E1:
   /raven/ptmp/jinma/overlap_1015/decks/e1_prod_rol2u.athinput). Add the grid fit, tlim (long-tlim
   rule), vet_gd_rotate_dt (~20 dt) and the seed amplitude/k band.
5. Run the background checklist on the full-size deck before any production submit.

## 4. Spec reference (`icstd/spec.py` DEFAULTS)

| block | keys | default / meaning |
|---|---|---|
| | name, mode | mode: `mesa` if mesa_profile else `global` |
| star | L/logL/`mesa`; M/GM/`mesa`; Teff/R_ph_rsun/R_ph_cm/`mesa` | M x 1.3271244e26, or `mesa` = G_code M_MESA |
| composition | X Y Z mixture a_metal tops_mixture | GS98 via planck_tools/mixspec.py; a_metal 16 (E1) |
| mesa_profile | path, format (mesa/cgs6), r_max_rsun | |
| r_in | rule (hp_below_cz / explicit), n_hp, r_rsun/r_cm, bottom_rsun, rlo_fac, cz_fmlt, explore_tmax, l_flat_tol | 4.5 H_p below the Fe-CZ base (mesa: MESA's CZ; the Fe peak if radiative); bottom 0.95 x1min |
| r_top | rtop_fac, col_top_fac | x1max 3 R_ph; column top 1.15 x1max |
| mlt | alpha, mlt_option, alt_scale_height, mltpp, use_mltpp, f1, f2, lambda1/2, beta1/2, dlambda, dbeta, xa1_center, max_logT, pp_alpha, pp_form | MESA r15140 defaults (sect. 2a): 2.0, Cox, true, mesa_gated, **false**, 1e-4, 1e-3, 1.0/0.5, 0.35/0.25, 0.1, 0.1, 0, 8 |
| grid | dlnr, drn_atm_rph, fl_step, tauj, rtol, rho_lo/hi, scheme | 2e-4, 2e-4, 50, 3, 1e-12; scheme agcar_ge (mltpp_inward = reproduce make_ic_agcar_mltpp.py only) |
| opacity | tops_hi/lo, rosseland/planck (explicit), ic_rosseland, ld_ext | ext2 recipe; ld_ext -21 |
| eos | dump (explicit) or logd_min/max, logt_min/max, dlog | -21..-1 (mesa: >= MESA base + 0.5 dex), tfloor rule..6.8 (mesa: >= MESA base T + 0.1 dex), 0.02 |
| background | rho_floor_rel, dmax_rel, dfloor_rel, sponge_rate_rel, sponge_r0_rel, tfloor_rel, wb_rmax_rel, pfloor_rel, vceil_rel, balance_rmax_rel (+ explicit rho_floor/tfloor/sponge_rate) | sect. 6 |
| ramp | start_turnovers, time_turnovers | 2, 10 |
| | allow_outside_tables, report_radii_rsun | reproduction/validation only; extra MESA comparison radii |
| config | cache, eos_binary | /raven/ptmp/jinma/ic_star_cache, bubbles_1012 athena_cpu_m4 |

## 5. Tables: one recipe, any composition (`icstd/tables.py`)

- **Opacity = the E1 ext2 recipe** (Viper lbv_1008/agcar/tables_ext/scripts/build_ext.py fergR):
  - TOPS/ATOMIC (high T) + Ferguson+2005 GS98 Rosseland and Planck (low T), linear in X and
    in log Z between the bracketing tables;
  - blended over log T 4.0-4.2, with the Saha-es + power-law low-rho extension;
  - grid log T 2.6-8.0 / log rho -21..0.
  - Cache `<cache>/tables/opac_<hash>/` (hash of X, Z, mixture string, recipe version) with
    mixture.txt, provenance.json and build_fields.npz (TOPS floors and the Ferguson box, used
    by the coverage map).
  - TOPS raw: `opacity.tops_hi/lo` copied into `<cache>/raw/tops/<hash>_{hi,lo}.dat`, or
    fetched (planck_tools/fetch_tops.py: hi 0.0005-1 keV x 1e-14..1 (71), lo 1e-21..1e-14
    (36)). **LANL was down on 10-15.** Without 'lo', the extension starts at log rho -14.
  - Ferguson raw: `<cache>/raw/ferguson05/{ross,planck}` (from Viper bench/m1_opac).
  - Checks:
    - X 0.36 Z 0.02 rebuilt = E1 rosseland_ext2 b28e97c5 and planck_ext2 dc482452, every
      value (E1's planck_ext2orion 493f9d01 differs in 4462/91357 values, the Orion low-rho
      stitch, which is not in the recipe);
    - X 0.7 Z 0.008 = bsg_newsetup_1015 ext2nolo, both tables, every value.
  - He-rich, non-GS98 mixtures: TOPS uses the given mixture (`tops_mixture`). Ferguson is
    GS98-scaled at (X, Z), which is an approximation, but only below log T 4.2.
- **EOS = the code's own dump** (1-cell shock_tube deck with eos_xh/yhe/a_metal and the range
  keys, `config.eos_binary`). Cache `<cache>/tables/eos_<hash>/`. Checks, byte-identical:
  - AG Car b5e87701;
  - BSG 193fb40b;
  - hg145 dd45fb60.
- **Hard range check**: FATAL if any column row lies outside the EOS dump (margin 0.2 dex rho /
  0.05 dex T) or the opacity grid. The tool never relies on clamping or extrapolation.

## 6. Rules in units of the star (hot-star calibrated; they reproduce E1)

| deck key | rule | AG Car A (E1) |
|---|---|---|
| x1min | n_hp H_p below the Fe-CZ base | 83 Rsun = 3.4 H_p (explicit); default 4.5 |
| x1max | rtop_fac R_ph | 3.0 |
| column floor | decade(4e-6 rho_ph) | 1e-16 (rho_ph 2.56e-11) |
| he_sponge_dmax | 2 x column floor | 2e-16 |
| dfloor | 1e-4 x column floor | 1e-20 |
| he_sponge_rate | decade(80 g/c_s,iso at R_ph) | 1e-3 /s (g/c_s 1.26e-5) |
| he_sponge_r0 | 0.9 x1max | |
| tfloor_kelvin | 100-K floor of 0.76 T_rad(x1max), T_rad = (E_col/a)^1/4 | 3000 K (T_rad 3956 K) |
| eos_logt_min | floor(10 log10 tfloor)/10 | 3.4 |
| pfloor | decade(6e-8 P_gas(floor, tfloor)) | 1e-12 |
| vceil | decade(4.3 v_esc(R_ph)) | 1e8 |
| wb_rmax | 1.1 R_ph | |
| he_ic_balance_rmax | 0.8808 R_ph | 342 Rsun |
| mlt_ramp_start / _time | 2 / 10 x t_Fe = H_p/v_MLT at the Fe opacity peak | 4.0297e6 / 2.0149e7 (E1 4.032e6 / 2.016e7) |
| implicit_flux_x1min | the column's F_r(x1min), log-linear in r (= he_star_m1's fine grid; the code fatals if the two differ by > 1e-3, he_star_m1.cpp ~1552) | 8.1426e12 |
| he_gm (+ he_gm_column) | global: GM point mass; mesa: G m(x1min) of the column + he_gm_column true | |

**Caveats for big, cool stars (RGs, RSGs).** These rules are tuned on hot, compact-ish
envelopes and are probably NOT good for them:
- **Floor and dmax** follow rho_ph. An RSG has an extended atmosphere, possibly dynamically
  supported, so tie the floor to the atmosphere at a chosen tau or radius instead.
- **Sponge rate** follows g/c_s at R_ph. With g ~ 1 cm/s^2, an absolute 1e-3 /s would freeze an
  RSG atmosphere.
- **tfloor**: 3000 K is near an RSG Teff. Below Teff ~5000 K, H2 and molecules matter (the EOS
  dumps have eos_h2 = false), Ferguson includes grains below ~2000 K, and the grey Hopf/RE
  atmosphere is crude (molecular, non-grey).
- **x1max**: extended atmospheres reach several R_ph, so choose rtop_fac from the atmosphere.

The diagnostics FLAG these regimes (Teff < 5000 K, the Ferguson/grain range, H_ph/R_ph > 0.03,
super-Eddington layers, x1max > 5 R_ph). These are flags only.

## 7. Keys every production deck must name (the deck block)

- **<mesh>**: x1min, x1max (+ the grid fit).
- **<hydro>**:
  - eos general/table, eos_h2 false, eos_ionization true, eos_metal_ionization false,
    eos_radiation false;
  - eos_xh, eos_yhe, eos_a_metal, eos_logd_min -21, eos_logd_max, eos_logt_min, eos_logt_max,
    eos_dlog;
  - wb_rmax, dfloor, tfloor_kelvin, pfloor, dfloor_keep_velocity, vceil.
- **<rad_m1>**: implicit_bc_x1min flux, **implicit_flux_x1min = F_r(x1min) of the column**,
  implicit_bc_x1max marshak, implicit_flux_x1max 0.
- **<problem>**:
  - he_gm (+ he_gm_column in mode mesa), he_ic_file, he_opac_table, he_planck_table,
    he_opac_logd_min, he_nfine 262144, he_tau_int 1, he_ic_cols 5;
  - **the IC chain**: he_ic_rad sph1d, he_ic_rad_gas re_atm, he_ic_rad_hse true,
    he_ic_rad_atm_force module, he_ic_rad_atm_post false, he_ic_eint_from_t true,
    he_ic_balance true, he_ic_balance_rmax, he_ic_equil true (+ dtfac 1, relax_n 40,
    relax_tol 1e-4, maxit 1, tol 1e-9, last_march true);
  - **the scaffold**: mlt_flux_frozen true, mlt_closure frozen, mlt_ramp_start, mlt_ramp_time;
  - **walls and top**: he_wall_zero_flux, he_wall_closed, he_bc_inner wall,
    he_bc_outer outflow, he_sponge_rate, he_sponge_r0, he_sponge_mode radial, he_sponge_dmax,
    he_seed_rlo/rhi, he_grid_dump.

## 8. Regression and validation

`regress/RESULTS.md` holds the numbers of 10-15/16 and `regress/run_regress.sh` the commands:
- AG Car A deep70 (global, agcar_legacy) vs /raven/ptmp/jinma/agcmodel_1012/ic/ic_agcar_A_ge_deep70.txt;
- hg145 ic2 (global, scheme mltpp_inward, allow_outside_tables) vs its ic2 file;
- mode mesa on 14p5M and the BSG: self-consistency against MESA (m_env, R_ph, T_ph, CZs);
- global-mode validation on 14p5M (gated MLT++) against its MESA model.

## Files

- `make_star_ic.py`: the driver.
- `icstd/`:
  - `micro.py`: EOS, opacity, MLT (KW and Cox), the MESA gate rule;
  - `column.py`: mode global, the gate solve;
  - `mesacol.py`: mode mesa;
  - `tables.py`, `saha_es.py`: tables;
  - `analyze.py`: rules and diagnostics;
  - `deck.py`: the deck block;
  - `diag.py`: the PNG;
  - `mesa.py`: profile readers.
- `specs/`:
  - `agcar_A_deep70`: the regression;
  - `bsg`;
  - `hegiant_14p5M`;
  - `hegiant_14p5M_global`: validation;
  - `hegiant_14p5M_hg145compat`;
  - `hegiant_15M`;
  - `hegiant_M2pt754`: reference.
