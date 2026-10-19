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

## 2b. Photosphere seam removal (10-16; spec `taper`, `seam`, `mesa_profile.top_fit`)

The hg145 and BSG gates failed in the atmosphere for the same reason, a seam at the
photosphere:
- the code's startup (he_ic_rad = sph1d + re_atm + he_ic_rad_hse) put tau 2/3 at a different
  radius than the column (hg145: 112.6 vs 113.9 Rsun; BSG: 51.5 vs 53.03 Rsun);
- the frozen F_MLT (~0.5 F for the He giant) ended abruptly there.

The standard now removes it:
1. **Taper** (`taper: {enable, tau_hi 30, tau_lo 3}`, on by default):
   - F_MLT/F x smoothstep(ln tau): 1 at tau >= 30, 0 at tau <= 3 (tau from the mesh top with the
     run kappa_R);
   - F_r = F (1 - F_MLT/F), so F_r + F_MLT = L holds exactly.
   - Why 30 -> 3: sph1d sets T to radiative equilibrium with weight 1 for tau < he_ic_rad_taua =
     30 (he_star_m1.cpp ~821). The taper therefore lives entirely in the code's full-RE zone, and
     radiation carries the whole flux by the Hopf join at tau 3, as it does in the RE atmosphere.
2. **Code-consistent atmosphere** (`seam: {deck, binary, tol_cells 0.1, maxit 8}`): the column
   is passed through the code's OWN startup solve:
   - the binary is the CPU build of the production commit, /raven/ptmp/jinma/icstd_1016/bin/athena_he_cpu_b16b66db;
   - it runs 0 cycles on the run deck with lateral 4x4, he_ic_equil and ic_rad off (they come
     after the dump), and he_ic_rad_dump;
   - the dumped column (sph1d E, RE T, HSE rho; F_r and F_MLT kept), with the taper re-applied,
     becomes the new column;
   - the loop iterates to the fixed point R_ph(column) = R_ph(sph1d) within 0.1 of the radial
     cell at R_ph. The code prints both: 'R_ph (tau_R 2/3, run kappa) rph0 -> rph1'.
   - The deck must have the spec's x1min and x1max, because the fine grid and tau are measured
     from the mesh top.
3. **Two-sided fit** for a truncated MESA profile (`mesa_profile.top_fit: two_sided`, the
   default; `shoot` = the old eps density jump):
   - outward: MESA T(P) + our HSE up to the profile top r_e;
   - inward: the mode-global surface column (Hopf + iso atmosphere, local MLT of the spec =
     MESA defaults, point mass G m(r_e)) from R_ph into [0.9 r_e, r_e];
   - **free parameter R_ph** (within 1 %): for each R_ph the fit point r* is where P_in = P_out,
     and R_ph is chosen so that T is continuous there too (rho then follows from the EOS). The
     continuity residuals dlnP, dlnT and dln rho are reported. If no root exists within 1 %, the
     smallest |dlnT| is used (reported).
4. **Results** (regress/RESULTS.md sect. 6):
   - **BSG:**
     - seam fixed point in 2 passes, R_ph 51.57 Rsun;
     - the run's startup reproduces it to 0.12-0.14 cell;
     - two-sided fit dlnT -5e-3, dln rho +3.8e-2 at R_ph 0.970 x 53.03;
     - gate (a) still fails: the oscillation survives without the seam, so it is likely
       physical/closure.
   - **AG Car A** if regenerated: +0.64 % R_ph.
   - **He giant 14p5M:** the code's own startup does not converge (a super-Eddington RE
     photosphere at the He/H opacity cliff), so `seam: {enable: false}` and
     `mesa_profile.top: surface` with the taper at tau 1 -> 0.1 (`extend: true`). This is a
     code-side limit, not an IC defect. **Superseded 10-17 by sect. 2d** (taper off, startup
     re_zone hopf, seam on).
5. Options:
   - `taper.extend` continues the scaffold outward at its tau_hi value before tapering;
   - `mesa_profile.top: surface` follows MESA's T(P) to MESA's surface zone (R_ph there, no
     Hopf layer), for stars whose MESA convective flux reaches the photosphere.
6. **Legacy**: `taper: {enable: false}, seam: {enable: false}` reproduces the old columns (AG Car A
   deep70 regression). `specs/agcar_A_deep70_seam.yaml` shows how the E1 IC would change if
   regenerated with the defaults.

### 2c. The code's startup solve for a super-Eddington photosphere (14p5M): diagnosis

**Status 10-17: solved without a code change by `startup.re_zone: hopf` (sect. 2d,
RESULTS sect. 9).** Options A-C below are kept for the record; A was tried (RESULTS 8b/8c,
branch ic-gammacap-1016) and does not converge.


**What fails.** he_star_m1.cpp ~813-1125 runs an outer Picard loop of two steps:
- (i) the E solve at fixed rho (HsSph1dPass VET passes, T = (E/a)^1/4 in the RE zone tau <
  taua);
- (ii) the HSE march at fixed T, with g_eff = g - kt F_r/c. Where Gamma = kt F_r/(c g) >
  he_ic_rad_hse_gmax (0.9), g_eff is set to (1 - gmax) g.

For 14p5M the photosphere lies on the He/H recombination opacity cliff (kappa 0.01 -> 5 between
7.5 and 20 kK, Gamma_L up to 6):
- a slightly hotter RE T raises kappa 100x, so Gamma > 0.9;
- the clamp gives an artificial 0.1 g scale height, so rho is pushed out, tau grows and RE T
  rises again. The loop has no contraction when dln kappa/dln T is large: 400 outer iterations,
  max |dln rho| 0.9-1.9, R_ph +1..+4 % per pass (seam passes), or the opposite collapse
  (rho x 0.006 above r_b) in the gate startup.
- hg145 ic3 "converged" only because a tau jump (0.76 -> > 3000) confined the solve to the thin
  top layer.

**Options** (each opt-in through a new problem key, default off = bitwise):
- **A. Scaffold takes the super-Eddington excess** (`he_ic_rad_hse_scaffold`, suggested):
  - in step (i)/(ii), cap the radiative flux of each fine cell at F_r <= gmax c g/kt(T, rho) and
    move the excess into the frozen scaffold (hM), so that F_r + F_MLT = L is kept;
  - the column's F_MLT and the frozen deposit then follow the code's own RE solution. This is
    what MESA does there: convection, F_MLT/F 0.76-0.95 at tau 30-100;
  - the change is local to the sph1d block (cap + bookkeeping, also written by he_ic_rad_dump);
  - about 150-250 lines; 1-2 days including the CPU seam passes (2-20 min each) and the gates.
- **B. Robust nonlinear solve:**
  - replace the Picard alternation by a damped coupled update (damp T as well as rho), or by a
    continuation in the flux (F_r from 0.5 F to F over the outer iterations);
  - report a fold, i.e. when no static solution exists, instead of clamping to 0.1 g;
  - this also answers whether a static radiative photosphere exists at all for this L/M;
  - 3-5 days; higher risk for the existing productions.
- **C. No code change:** a tool pre-check (Gamma_rad of the RE T > ~0.9 near the photosphere)
  that flags the star and keeps the scaffold through the photosphere (mesa_profile.top =
  surface, taper above tau 1). The seam then remains, as now.

**Risks.**
- A changes the IC energy bookkeeping wherever Gamma > gmax. AG Car B (Gamma 1.3 at tau 0.5-3)
  would change if opted in; E1 and the BSG (Gamma <= 0.8 and 0.67) would not.
- If no static solution exists, the right IC is not static: it is an outflow, and neither A nor
  B makes it static.
- B can change convergence for every star.

**Gates.**
1. Option off: E1 and BSG startup bitwise (sph1d dump and R_ph line identical).
2. 14p5M seam fixed point: HSE dln rho < 1e-6 and R_ph pass to pass < 0.1 cell.
3. Gate (a) 14p5M on 704 and 1408, band 0.99-1.1 R_ph: |<v_r>| < 1 km/s, drho < 10 %, no growth.
4. AG Car B startup with the option on vs off (documented change).

### 2d. MESA-anchored startup: `startup.re_zone: hopf` (10-17; RESULTS sect. 9)

**Why the startup failed (10-17 diagnosis).**
- The startup re-solves radiative equilibrium (RE) with the column's FROZEN F_r(r) for tau_R < 30
  and blends it to the column's T up to tau_R 3000. For 14p5M that zone is the He/H opacity cliff
  (kappa 0.01 -> 5 at 10-20 kK, dln kappa/dln T up to ~13).
- There the RE T at fixed rho follows from tau(r) = int kappa(T) rho dr. The Picard loop
  (VET passes with kappa lagged, then HSE at fixed T) has two attracting branches:
  - a cool, collapsed one: lower T, kappa 10x lower, tau lower;
  - a hot, inflated one.
  The column, i.e. MESA's structure, sits between them and is unstable under this iteration.
- This holds with OUR kappa and also with MESA's kappa (diagnostic table, RESULTS 9):
  - with the column's own scaffold and no taper: 400 outer iterations, inflated, +25 cells;
  - with MESA's kappa: converged, but on the collapsed branch (-14 cells, rho x 0.01 at R_ph).
  The taper of 2b (F_MLT -> 0 by tau 3, where MESA still has F_MLT/F 0.2-0.9) made it worse.
- MESA itself carries 76-99 % of F by convection at tau 30-3000. The static structure exists only
  with that convective flux. The column's scaffold F_MLT/F = 1 - nabla_MESA/nabla_rad(ours)
  re-expresses it in our kappa, so its Gamma_rad = kappa F_r/(c g) equals MESA's (max 1.01).
  MESA's own F_MLT/F with OUR kappa would give Gamma_rad 1.49 (no static solution with that
  split; our kappa is 1.4-1.5x MESA's at 24-33 kK).

**The option** (tool only; the code keys he_ic_rad_taua/taub already exist; default 'deck' =
unchanged, bitwise):
- `startup: {re_zone: hopf}`: he_ic_rad_taua = max(2/3, tauj/3) = 1 and he_ic_rad_taub = tauj = 3
  (tau_R of the run kappa from the mesh top). They go into the seam passes and the deck block.
- The startup then re-solves RE only in the column's photosphere layer. Below tau 3 it keeps the
  column's T (MESA's T(P)) and the scaffold, and its HSE march starts at r_b(tau 3). The
  Gamma > 0.9 layer (tau ~3500-4000) is not marched.
- `taper: {enable: false}` is allowed with the seam: the scaffold stays the column's own, and each
  pass sets F_r = F (1 - F_MLT/F) exactly on the dump rows.
- The seam passes use he_nfine 262144 (= the deck block; the hg145 seam deck had 2^20).
- The diagnostics print Gamma_rad of the scaffold vs that of MESA's own split with our kappa.

**When to use it:** mode mesa, when MESA's convective flux reaches tau ~ tauj and the default
startup does not converge (the 2c symptoms). E1/BSG/AG Car keep the default.

**14p5M** (`specs/hegiant_14p5M.yaml`; the old cap spec is `hegiant_14p5M_cap.yaml`):
- seam fixed point in 2 passes: R_ph 7.96114e12 -> 7.95671e12 (-0.58 cell) -> 7.95704e12
  (+0.043 cell) = 114.38 Rsun;
- HSE 5.8e-7, no Gamma clamp;
- CPU = A100 startup (R_ph identical, dump 3e-10); gate (a') FAILS: the run restructures the cliff
  top within 0.25 d, also with MESA-like kappa and at 4x radial resolution (RESULTS sect. 9e).

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
6. Run the atmosphere gate (a') (sect. 8a): the gate run on the deck grid and on its
   2x-atmosphere grid, then `regress/gate_aprime.py DECK_RUN PAIR_RUN`.

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
| startup | re_zone | 'deck' (the deck's he_ic_rad_taua/taub, E1 30/3000); 'hopf': RE only at tau_R < tauj (sect. 2d) |
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

## 8a. Atmosphere gate (a') (ADOPTED by the USER 2026-10-17; `regress/gate_aprime.py`)

Gate (a') replaces the absolute atmosphere thresholds of gate (a) (`regress/gate_atm.py`:
0.99-1.1 R_ph |<v_r>| < 1 km/s, |drho| < 10 %), which no resolution of the BSG new setup can pass
(its photosphere carries a converged, bounded, self-excited acoustic resonance; bsgres-1017,
/raven/ptmp/jinma/bsgres_1017/NOTE.md). Runs: seed 0, 16x16 lateral, >= 8 d, 0.05-d hydro_w +
m1 dumps, on the deck grid AND on its 2x-atmosphere grid. All items must pass:

| item | criterion |
|---|---|
| 1 resolution | >= 4.4 cells per atmospheric (density) scale height over 0.95-1.1 R_ph on the deck grid: min of H_rho/dr, H_rho = -dr/dln<rho> of the t = 0 startup atmosphere, floor-plateau cells (rho < 1.5 x min) excluded |
| 2 convergence | the SAME floor rule on both runs (column floor = min rho of he_ic_file, startup plateau within 1 %, dfloor and he_sponge_dmax identical), and deck vs 2x-atmosphere grid, saturated t >= 4 d, band 0.99-1.1 R_ph: band mean rms <v_r> within 15 %, time-rms profile L2 diff <= 10 % (L2 over 0.99-1.1 R_ph, pair profile interpolated to the deck radii, \|\|p_d - p_p\|\|/\|\|p_p\|\|), dominant period (FFT, median of the 5 top-rms radii) within 10 % |
| 3 bounded | no growth (last-third max <= 1.5 x middle-third max) and saturated band rms <= c_s,iso(R_ph) |
| 4 star unaffected | 0.5-0.9 R_ph max \|<v_r>\| < 1 km/s over the whole run; 0.9-0.99 R_ph < 3 km/s on the SATURATED state t >= 4 d (the growth-phase max is printed for information); R_ph(tau_R 2/3) within +-2 % of t = 0 at every dump; L_top/L_in within +-5 % for t >= 2 d; \|dM/M\| < 1e-6; Mdot_top < 1e-8 Msun/yr |
| 5 robustness | summary NON-CONVERGED = 0, no FATAL/DIVERGED/NaN, rc 0, dt_min >= 0.25 dt_median (hst dt WITHOUT the startup row and the final step truncated to tlim; dt_retry retries reported) |

The three readings in items 2 and 4-5 (0.9-0.99 R_ph on the saturated state, the dt rows, the L2
band) were confirmed by the USER on 2026-10-17.

**Floor.** The standard floor rule is kept: column floor decade(4e-6 rho_ph) (1e-15 for the BSG),
dfloor 1e-4 x, he_sponge_dmax 2 x the column floor. The oscillation AMPLITUDE depends on the floor
level, because the floor plateau sets the top of the oscillating layer (the oscillation fills the
hydrostatic atmosphere below the plateau and dies 0.02-0.03 R_ph into it): BSG g480 band rms
9.2 / 6.0 / 3.3 km/s for floors 1e-16 / 1e-15 / 1e-14, while the period changes by <= 3 % and the
star below the photosphere is unchanged (RESULTS sect. 9). (a') comparisons are therefore valid only
with the same floor rule on both runs; gate_aprime.py FAILs item 2 otherwise.

Items 3-5 are checked on both runs. BSG production must use a g480-class grid (the P grid,
~1 cell per H_rho, fails items 1 and 2); g480 itself has 4.24 cells per H_rho with the seam-fixed IC,
so the candidate BSG production grid is the refined r496 (plateau dr 0.0319 Rsun, 4.46 cells;
/raven/ptmp/jinma/bsgres_1017/decks/bsgres_r496.athinput, pair r808). Results: regress/RESULTS.md sect. 9.

## Files

- `make_star_ic.py`: the driver.
- `regress/gate_aprime.py`: the adopted atmosphere gate (a') (sect. 8a); `regress/gate_atm.py`:
  the old gate (a) table (diagnostic).
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
  - `hegiant_14p5M`: MESA-anchored startup (sect. 2d);
  - `hegiant_14p5M_cap`: the option-A gamma-cap attempt (superseded);
  - `hegiant_14p5M_global`: validation;
  - `hegiant_14p5M_hg145compat`;
  - `hegiant_15M`;
  - `hegiant_M2pt754`: reference.
