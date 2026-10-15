# tools/ic_star: the standard initial condition for he_star_m1 envelope runs

ONE pipeline for every stellar-envelope run (BSG, He giants, AG Car / LBVs, future stars):

```
star spec (YAML/JSON) --make_star_ic.py--> tables (cached by composition hash)
                                         -> IC column  ic_<name>.txt   (+ provenance header)
                                         -> deck block deck_<name>.athinput
                                         -> diagnostics diag_<name>.txt / .png, prov_<name>.json
```

```bash
cd tools/ic_star
export OPENBLAS_NUM_THREADS=1        # login nodes: OpenBLAS thread limit
PY=/u/jinma/miniconda3/envs/magritte/bin/python3      # numpy 1.23, scipy 1.10, PyYAML, mpl
$PY make_star_ic.py specs/bsg.yaml --out /raven/ptmp/jinma/ic_star_runs/bsg
# regression (AG Car A deep70, the E1 production IC):
$PY make_star_ic.py specs/agcar_A_deep70.yaml --out /raven/ptmp/jinma/ic_star_runs/agcar_A_deep70 \
    --compare /raven/ptmp/jinma/agcmodel_1012/ic/ic_agcar_A_ge_deep70.txt
```

A run takes 5-25 min on a shared login node (LSODA with a Python right-hand side; the
`hp_below_cz` rule adds one exploratory integration).

## 1. The method (the AG Car method, chosen by the user as the standard)

`icstd/column.py`, = `make_ic_agcar_ge.py` (Viper lbv_1008/agcar/geos, the AG Car A/B ICs) and
`make_ic_agcar_deep.py` (agcmodel-1012, the E1 deep wall), same equations:

1. R_ph (tau_R = 2/3) given, or sqrt(L / 4 pi sigma Teff^4). T(tau) grey Hopf (Mihalas fit).
2. Above R_ph: isothermal atmosphere at T_ph = T_Hopf(2/3), gravity only, rho_ph SHOT so that
   its tau(R_ph) = 2/3, down to the column floor; constant floor above, to the column top.
   E above R_ph is diluted, E = (F/c)(1 + (c a T_ph^4/F_ph - 1)(R_ph/r)^2). In the run,
   `he_ic_rad = sph1d` + `he_ic_rad_gas = re_atm` + `he_ic_rad_hse` replace this layer by the RE
   atmosphere of the code's own transport, and `he_ic_equil` makes the discrete state static.
3. R_ph inward to tau = 3: Hopf T(tau) + gas HSE with the radiation force of F.
4. tau = 3 inward to the column bottom: dP_tot/dr = -rho GM/r^2, dlnT/dlnP_tot = nabla, where
   nabla is local Bohm-Vitense MLT (alpha_MLT 1.5) with MLT++ (Paxton+2013 sect. 7.2, MESA
   adjust_gradT_excess: nabla <- af nabla + (1 - af) nabla_ad where nabla - nabla_ad > f1,
   af = f2 + (1 - f2)(1 - pp_alpha); f1 1e-4, f2 1e-3; pp_alpha 1 gives af = f2), and
   nabla = nabla_rad where radiative. F_r = F nabla/nabla_rad, F_MLT/F = 1 - nabla/nabla_rad,
   F = L/(4 pi r^2). L is constant and GM is a point mass (`he_gm`).
5. EOS: the code's general-EOS table dump (bicubic spline, clamped to the grid) + radiation.
   kappa_R: the run's Rosseland table, bilinear in (log T, log rho), extended below its grid as
   `he_opac_logd_min` does.
6. Rows: log grid dln r = 2e-4 from the column bottom up to the tau-3 join, with x1min as an exact
   row; Hopf and atmosphere rows every 2e-4 R_ph; floor rows every 50 x that. Columns:
   `r rho eint F_r(radiative) E T F_MLT/F` (he_ic_cols = 5, 7 columns; mlt_flux_frozen reads
   column 7).

Schemes (`grid.scheme`): `agcar_ge` is THE STANDARD. `mltpp_inward` exists only to reproduce
`make_ic_agcar_mltpp.py --mode inward` files bit for bit (hg145 ic2). It uses the same equations
with a different row placement and rtol 1e-9, and MLT++ in the algebraically identical 'f2'
form. `r_in.legacy_deep` reproduces the two-stage AG Car A deep70 file. Do not use either for
new stars.

## 2. Step by step for a new star

1. Write `specs/<name>.yaml` (see section 3; copy `bsg.yaml`). Give L, M (or GM), Teff or R_ph,
   X Y Z (+ GS98), and optionally a MESA profile.
2. Run `make_star_ic.py`. Tables are built once per composition (section 4).
3. Read `diag_<name>.txt` and check:
   - coverage: column inside the EOS and opacity grids, 0 clamped calls, opacity source map;
   - the wall: H_p buffer >= 4-5 below the Fe-CZ base, buffer RADIATIVE
     (max nabla_rad/nabla_ad < 1), L flat at x1min (MESA);
   - L consistency: 4 pi r^2 (F_r + F_MLT) = L, and F_diff/F_r on the radiative rows;
   - background: floor-layer extent, sponge ACTIVE (floor < dmax), tfloor < T_rad(x1max);
   - regime flags;
   - MESA comparison: T, rho ratios and envelope mass. A large mismatch means the inward
     MLT++ column found another branch (see section 7). Stop there and ask; do not patch the
     method per star.
4. Paste `deck_<name>.athinput` into a copy of the current production deck (E1:
   /raven/ptmp/jinma/overlap_1015/decks/e1_prod_rol2u.athinput). Each key replaces the line of
   the same name. Add the grid (nx1, f_stretch_*), tlim (long-tlim rule), vet_gd_rotate_dt
   (~20 dt) and the seed amplitude/k band from the production deck.
5. Before any production submit, run the background checklist on the full-size deck
   (check-background-before-prod).

## 3. Spec reference (defaults = AG Car E1 values; `icstd/spec.py` DEFAULTS)

| block | key | meaning / default |
|---|---|---|
| | name | output tag |
| star | L / logL | erg/s, or log L/Lsun (Lsun 3.828e33), or `mesa` (photosphere_L) |
| star | M / GM | Msun (GM = M x 1.3271244e26), or GM [cgs], or `mesa` (G_code x M_MESA) |
| star | Teff / R_ph_rsun / R_ph_cm | R_ph given, or from L and Teff; `mesa` = photosphere_r |
| composition | X Y Z mixture a_metal | GS98 (TOPS mixture from planck_tools/mixspec.py), or `custom` + `tops_mixture` (number fractions); a_metal = eos_a_metal (16 = E1 inert metals; He giants 1.77-1.80 fully ionized) |
| mesa_profile | path, format | `mesa` (profile*.data) or `cgs6` (r T rho m kappa F_r); comparison + optional L/M/R source |
| r_in | rule | `hp_below_cz` (default; n_hp 4.5 below the Fe-CZ base) or `explicit` (r_rsun / r_cm) = x1min |
| r_in | rlo_fac / bottom_rsun | column bottom = 0.95 x1min (ghost margin) unless given |
| r_in | cz_fmlt, explore_tmax, l_flat_tol | 0.01 (CZ definition), 3e6 K (explore depth), 0.01 (L flatness vs MESA) |
| r_top | rtop_fac, col_top_fac | x1max = 3.0 R_ph (E1); column top 1.15 x1max |
| mlt | alpha, f1, f2, pp_alpha | 1.5, 1e-4, 1e-3, 1.0 (AG Car A used pp_alpha 0.970755, the MESA (beta, lambda) fixed point) |
| grid | dlnr, drn_atm_rph, fl_step, tauj, rtol | 2e-4, 2e-4, 50, 3, 1e-12 |
| grid | rho_lo, rho_hi | the rho_ph shooting bracket (1e-14 .. 1e-7) |
| opacity | tops_hi, tops_lo | raw TOPS files (else fetched from LANL); `rosseland`/`planck` = explicit tables (no build); `ic_rosseland` = a different table for the column only (AG Car A legacy) |
| eos | dump or logd_max, logt_max, dlog | explicit dump, or the recipe (-21..-1, logt_min from tfloor, ..6.8, 0.02) |
| background | see section 5 | rules; any key may be given explicitly (rho_floor, tfloor, sponge_rate) |
| ramp | start_turnovers, time_turnovers | 2 and 10 Fe-peak turnovers (E1) |
| config | cache, eos_binary | /raven/ptmp/jinma/ic_star_cache, bubbles_1012 athena_cpu_m4 |

## 4. Tables: one recipe, any composition (`icstd/tables.py`)

- **Opacity = the E1 ext2 recipe** (Viper lbv_1008/agcar/tables_ext/scripts/build_ext.py fergR):
  - high T: TOPS/ATOMIC Rosseland and Planck for the mixture. 'hi' covers 0.0005-1 keV x
    rho 1e-14..1 (71 points); 'lo' covers rho 1e-21..1e-14 (36 points) and is optional
    (without it the extension starts at log rho -14);
  - low T: Ferguson+2005 GS98, linear in X and in log Z between the bracketing tables;
  - blended over log T 4.0-4.2, with the Saha-es + power-law low-rho extension;
  - grid log T 2.6-8.0 / log rho -21..0.
  - Cache: `<cache>/tables/opac_<hash>/` with mixture.txt, provenance.json and build_fields.npz
    (TOPS floors, Ferguson box, used for the coverage map). The hash covers X, Z, the mixture
    string and the recipe version.
  - Raw inputs: `<cache>/raw/ferguson05/{ross,planck}` (copied from Viper bench/m1_opac) and
    `<cache>/raw/tops/<hash>_{hi,lo}.dat`. TOPS is copied from `opacity.tops_hi/lo`, or fetched
    with data/stellar_opac/planck_tools/fetch_tops.py. **The LANL service was down on 10-15**
    (HTTP 500), so new mixtures need a retry.
  - Check: X 0.36 Z 0.02 rebuilt by this code equals the E1 tables as data. The rc_1009
    rosseland_ext2 (b28e97c5) and planck_ext2 (dc482452) match in every value. E1's
    planck_ext2orion (493f9d01) differs in 4462 of 91357 values: that is the Orion low-rho
    stitch, which is not part of the recipe.
- **EOS = the code's own dump**: a 1-cell `shock_tube` deck with eos_xh / eos_yhe / eos_a_metal
  and the range keys, run with `config.eos_binary`. Cache: `<cache>/tables/eos_<hash>/`.
- Coverage report (diag): the column vs the EOS grid and the deck's EOS range (he_star_m1 fatals
  outside), the builder's clamped calls, the column vs the Rosseland grid, and the source of
  every row (TOPS data / TOPS extension / blend / Ferguson data / Ferguson extension).

## 5. Rules in units of the star (hot-star calibrated; they reproduce the E1 values)

| deck key | rule | AG Car A (E1) |
|---|---|---|
| x1min | `hp_below_cz`: n_hp H_p below the base of the Fe CZ (F_MLT/F > 0.01 around the max of kappa_R at 1e5-1e6 K) | 83 Rsun = 3.4 H_p (explicit); the new default is 4.5 |
| x1max | rtop_fac R_ph | 3.0 |
| column floor | decade(4e-6 rho_ph) | 1e-16 (rho_ph 2.56e-11) |
| he_sponge_dmax | 2 x column floor | 2e-16 |
| dfloor | 1e-4 x column floor | 1e-20 |
| he_sponge_rate | decade(80 g/c_s,iso at R_ph) | 1e-3 /s (g/c_s 1.26e-5 /s) |
| he_sponge_r0 | 0.9 x1max | 0.9 |
| tfloor_kelvin | 100-K floor of 0.76 T_rad(x1max), T_rad = (E_col/a)^1/4 | 3000 K (T_rad 3956 K) |
| eos_logt_min | floor(10 log10 tfloor)/10 | 3.4 |
| pfloor | decade(6e-8 P_gas(floor, tfloor)) | 1e-12 |
| vceil | decade(4.3 v_esc(R_ph)) | 1e8 |
| wb_rmax | 1.1 R_ph | 1.1 |
| he_ic_balance_rmax | 0.8808 R_ph (inside the MLT++ adiabat, below the tau-3 join) | 342 Rsun |
| mlt_ramp_start / _time | 2 / 10 x t_Fe, t_Fe = H_p/v_MLT (textbook bubble speed) at the Fe opacity peak | 4.032e6 / 2.016e7 s |
| implicit_flux_x1min | column F_r(x1min), log-linear in r = he_star_m1's fine grid; the code fatals if the two differ by > 1e-3 (he_star_m1.cpp ~1552) | 8.1426e12 |
| he_seed_rlo / rhi | 1.007 x CZ base / 0.9786 R_ph | 146 / 380 Rsun |

**Caveats for big, cool stars (red giants, RSGs).** These rules are calibrated on hot,
compact-ish envelopes (AG Car Teff 9000 K, BSG 15.5 kK, He giants 7.6-14 kK) and are probably
NOT right for RGs/RSGs:
- **Column floor and dmax** follow rho_ph. An RSG photosphere is ~1e-9..1e-10 g/cc with an
  extended, possibly dynamically supported atmosphere, so a floor 4e-6 rho_ph may cut through gas
  that matters. Tie the floor to the RE/extended atmosphere at a chosen tau or radius instead.
- **Sponge rate** scales with g/c_s at R_ph. An RSG has g ~ 1 cm/s^2, so the absolute 1e-3 /s
  would freeze its atmosphere. Check that the rule's value stays a small multiple of the local
  dynamical rate wherever floor gas lives.
- **tfloor**: 3000 K is close to an RSG Teff. The rule ties tfloor to T_rad(x1max), but for
  Teff < ~5000 K it must stay below the coolest physical gas: H2/molecules matter there, the
  EOS dumps have eos_h2 = false, and Ferguson's tables include grains below ~2000 K.
- **x1max**: RSG extended atmospheres reach several R_ph, so 3 R_ph may be too small. Choose
  rtop_fac from the atmosphere, not from E1.
- The grey Hopf/RE atmosphere and the isothermal cap are crude for molecular, strongly non-grey
  atmospheres.

The diagnostics FLAG these regimes: Teff < 5000 K, T or tfloor in the Ferguson/grain range, an
extended atmosphere (H_ph/R_ph > 0.03), super-Eddington layers, and x1max > 5 R_ph. These are
flags only; there is no RSG implementation.

## 6. Keys every production deck must name (the deck block)

<mesh> x1min, x1max (+ the grid fit); <hydro> eos general/table, eos_h2 false,
eos_ionization true, eos_metal_ionization false, eos_radiation false, eos_xh, eos_yhe,
eos_a_metal, eos_logd_min -21, eos_logd_max, eos_logt_min, eos_logt_max, eos_dlog, wb_rmax,
dfloor, tfloor_kelvin, pfloor, dfloor_keep_velocity, vceil; <rad_m1> implicit_bc_x1min flux,
**implicit_flux_x1min = F_r(x1min) of the column**, implicit_bc_x1max marshak,
implicit_flux_x1max 0; <problem> he_gm, he_ic_file, he_opac_table, he_planck_table,
he_opac_logd_min, he_nfine 262144, he_tau_int 1, he_ic_cols 5, **the IC chain**: he_ic_rad sph1d,
he_ic_rad_gas re_atm, he_ic_rad_hse true, he_ic_rad_atm_force module, he_ic_rad_atm_post false,
he_ic_eint_from_t true, he_ic_balance true, he_ic_balance_rmax, he_ic_equil true (+ dtfac 1,
relax_n 40, relax_tol 1e-4, maxit 1, tol 1e-9, last_march true), **the scaffold**:
mlt_flux_frozen true, mlt_closure frozen, mlt_ramp_start, mlt_ramp_time, **walls/top**:
he_wall_zero_flux, he_wall_closed, he_bc_inner wall, he_bc_outer outflow, he_sponge_rate,
he_sponge_r0, he_sponge_mode radial, he_sponge_dmax, he_seed_rlo/rhi, he_grid_dump.

## 7. Known limits of the method (read before using it on a new class of star)

- **Super-Eddington convective envelopes (He giants).** For 14p5M, hg145 ic2 with the
  `mltpp_inward` scheme and the interim table gives Gamma_L up to ~5. MLT++ pins nabla to
  nabla_ad from the surface down, and the inward adiabat lands on a far too hot and dense
  branch: T(8 Rsun) 3.9e6 vs MESA 9.0e5 K, and an envelope mass of 3e4 Msun. The standard
  scheme is the same model and finds the same branch (diag of specs/hegiant_14p5M.yaml). The old
  He giant M2pt754 was built from MESA's T(r) for the same reason; it is not in thermal
  equilibrium (L(3 Rsun) = 418 L_ph). Whether the AG Car method can be used for He giants at
  all is a physics decision for the user (e.g. MLT++ controls), not a per-star patch.
- The column is in thermal equilibrium with the table it was built with (F_r + F_MLT = L). A
  different table in the deck (AG Car A: column on plain TOPS, deck on ext2) moves heat until
  the run relaxes. The standard builds the column with the deck's table.
- L is constant: no nuclear or thermal sources in the envelope. With a MESA profile, the wall
  must sit where L_MESA(r) is flat (checked, 1 %).

## 8. Regression and reproduction checks

See `regress/README` for the commands. Results (10-15):
- **AG Car A deep70** (specs/agcar_A_deep70.yaml vs
  /raven/ptmp/jinma/agcmodel_1012/ic/ic_agcar_A_ge_deep70.txt, md5 of its data in diag): see
  regress/RESULTS.md.
- **hg145 ic2** (specs/hegiant_14p5M_hg145compat.yaml vs
  /raven/ptmp/jinma/hegiant145_1015/ic2/ic_hg145_mltpp.txt): see regress/RESULTS.md.
- Tables: E1 ext2 Rosseland and Planck reproduced as data (section 4). The EOS dump recipe
  reproduces the AG Car dump b5e87701 (regress/RESULTS.md).

## Files

- `make_star_ic.py`: the driver.
- `icstd/micro.py`: EOS, opacity reader, MLT/MLT++.
- `icstd/column.py`: the column.
- `icstd/tables.py`, `icstd/saha_es.py`: the tables.
- `icstd/analyze.py`: rules and diagnostics.
- `icstd/deck.py`: the deck block.
- `icstd/diag.py`: the PNG.
- `icstd/mesa.py`: MESA readers.
- `specs/`:
  - `agcar_A_deep70`: the regression;
  - `bsg`: Ma+2026, X 0.7 Z 0.008;
  - `hegiant_14p5M`: r_in 8 Rsun;
  - `hegiant_14p5M_hg145compat`: hg145 ic2, bit for bit;
  - `hegiant_15M`;
  - `hegiant_M2pt754`: reference only.
- `regress/`: the check commands and results.
