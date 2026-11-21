# NOTE 2026-11-21hz Caltech -> viper + orion: R136a1 ic_star spec / top mode (answer to b43923f8f)

Read-only look at `/resnick/groups/carnegie_poc/jingze/r136a1_1104/` on Caltech. Verbatim copies of the spec, the
resolved spec (provenance), both logs, the run scripts and the seam script are in
`docs/handover/caltech_results/R136A1_IC/` (files named below).

**Short answer: R136a1 used the DEFAULT `top: hopf` (no `top`, `top_fit` or `r_max_rsun` key in the spec). The tau-3 join
was reached INSIDE MESA's T(P) range (no extension, no two-sided fit), and the outward Hopf layer reached tau 2/3 although it is
super-Eddington all the way (Gamma ~1.9 at r_j -> 5.4 at R_ph). It is only 0.0034 Rsun thick (rho ~1e-8, so the
remaining tau 3 -> 2/3 is used up almost at once).**

Identity check: `md5sum r136a1_1104/ic/ic_r136a1_B.txt` = `ec3c8a81dc1bfe6101bb84c88807f53d`, the same as
`r136a1_1104/ic/BB/col_3.txt` (the seam loop output). Chain: `ic/runAB.sh` (make_star_ic.py, spec `r136a1B_noseam.json`,
seam off) -> `ic/AB/ic_r136a1B_noseam.txt` -> `ic/runBB.sh` (`tools/seam_only_B.py`) -> `ic/BB/col_3.txt` -> copied to
`ic/ic_r136a1_B.txt`.

Note on the header: it says `tools/seam_only.py`, but the production column came from **`tools/seam_only_B.py`**. That file is a copy
of seam_only.py with the AB/BB paths, the B deck, nx1 831 and the explicit wall flux. The header string was not updated
in the copy. `seam_only.py` itself built the earlier first-look column `ic/ic_r136a1.txt` (ic/A -> ic/B, x1max 5.2178e12, nx1 512).

## 1. Spec

File: `/resnick/groups/carnegie_poc/jingze/r136a1_1104/ic/r136a1B_noseam.json` (copy: `R136A1_IC/r136a1B_noseam.json`),
full text:

```json
{
 "_comment": "r136a1_1104 (USER 11-04): R136a1, Keszthelyi+2025 TS (MLT++) v200 281.8 Msun midMS profile, tools/ic_star mode mesa (RSG E / HG / bsgsp way). Pass A (seam off) -> grid design; pass B (seam on, deck = the r136a1 deck) = production column. Opacity: Rosseland = MESA a09 X0.58 Z0.0064 (TOPS same-composition outside MESA's log R box), Planck = TOPS X0.58 Z0.0064 AGSS09 (LANL 11-04). x1min = r(log T 6.0) of the profile, x1max = 3 R_ph(MESA). DEEP (USER 11-04 late): x1min 19.45 Rsun (log T 6.50, MESA radiative, below the MLT++ envelope CZ base 20.6 Rsun), x1max 100 Rsun, grid leanB (USER: variant B).",
 "name": "r136a1B_noseam",
 "star": {"L": "mesa", "R_ph": "mesa", "M": "mesa"},
 "composition": {"X": 0.58, "Y": 0.4136, "Z": 0.0064,
  "mixture": "AGSS09 (TOPS number fractions in tops/mix_x0.58_z0.0064_a09.txt)", "a_metal": 16.0},
 "mesa_profile": {
  "path": "/resnick/groups/carnegie_poc/jingze/r136a1_1104/model/r136a1_ke25_282_v200_mltpp_midMS.data",
  "format": "mesa"},
 "opacity": {
  "rosseland": "/resnick/groups/carnegie_poc/jingze/r136a1_1104/opac/rosseland_mesakap_a09_x0.58_z0.0064_topsbase_T7.txt",
  "planck": "/resnick/groups/carnegie_poc/jingze/r136a1_1104/opac/planck_tops_x0.58_z0.0064_a09_T7.txt",
  "ld_ext": -21.0},
 "eos": {"logd_min": -21.0, "logd_max": -1.0, "logt_min": 3.3, "logt_max": 6.8, "dlog": 0.02,
  "binary": "/resnick/groups/carnegie_poc/jingze/builds_rstfixx/bin/athena_cpu_built_in_pgens_prod_18bfdeb9_fmadoff"},
 "r_in": {"rule": "explicit", "r_rsun": 19.45, "bottom_rsun": 18.5},
 "r_top": {"x1max_cm": 6957000000000.0, "col_top_fac": 1.15},
 "startup": {"re_zone": "hopf"},
 "taper": {"enable": false},
 "seam": {"enable": false},
 "config": {"cache": "/resnick/groups/carnegie_poc/jingze/r136a1_1104/ic_cache",
  "eos_binary": "/resnick/groups/carnegie_poc/jingze/builds_rstfixx/bin/athena_cpu_built_in_pgens_prod_18bfdeb9_fmadoff"}
}
```
(Whitespace compacted here; the copied file is byte-for-byte.)

The keys you asked about (facts from the spec, plus the code defaults at the build commit 83d97f8c9, make_star_ic.py lines 276-277):
- `mesa_profile.top`: **absent -> default `hopf`** (`s['mesa_profile'].get('top', 'hopf')`).
- `mesa_profile.top_fit`: absent -> default `two_sided`. **It was never used**: the two-sided branch only runs when the
  MESA profile ends before tau_j, and here the tau_j event fired first (section 2).
- `mesa_profile.r_max_rsun`: **absent (None) -> no cut.**
- `taper`: `enable false` (resolved: tau_lo 3, tau_hi 30, unused).
- `seam`: `enable false` in make_star_ic. The seam was run separately with seam_only_B.py (section 3).
- `startup.re_zone`: `hopf` -> deck `he_ic_rad_taua 1`, `he_ic_rad_taub 3`.
- `grid.tauj`: not in the spec -> default **3.0**. Resolved grid block (prov json): `dlnr 0.0002, drn_atm_rph 0.0002,
  fl_step 50, tauj 3.0, rtol 1e-12, rho_lo 1e-14, rho_hi 1e-07, atm_rad 0, scheme agcar_ge`.
- opacity: the block above, i.e. explicit tables with no build. Rosseland = MESA kap a09 (OPAL + Ferguson) X0.58 Z0.0064, using TOPS of the same
  composition outside MESA's log R box (md5 c5f783e0). Planck = TOPS X0.58 Z0.0064 AGSS09 (md5 06dbbaea). `ld_ext -21`.
- MLT (resolved defaults): alpha 2, Cox, alt_scale_height, `mltpp mesa_gated`, `use_mltpp false`.
- Full resolved spec with all defaults: `R136A1_IC/prov_r136a1B_noseam.json`.

## 2. make_star_ic.py log

File: `r136a1_1104/ic/runAB.log` (copy `R136A1_IC/runAB.log`, full). The run was tools/ic_star at git 83d97f8c9c6c4376 and
MESA profile md5 0f452a80. Relevant lines, verbatim:

```
star r136a1B_noseam: L 2.60446120e+40 erg/s (log L/Lsun 6.8327; MESA photosphere_L), GM 3.10569999e+28 (G_code x MESA total mass), R_ph 1.73269051e+12 cm = 24.9057 Rsun (MESA photosphere_r), Teff 59069.54 K, T_ph(Hopf, tau 2/3) 59403.16 K
x1min 1.3531365000e+12 (19.4500 Rsun, explicit); column bottom 1.2870450000e+12 (18.5000 Rsun); x1max 6.9570000000e+12 (4.0151 R_ph); column top 8.000550e+12
  mode mesa: tau-3 join at r_j 1.74174603e+12 (1.0052 R_MESA), T 76094.28 K, rho 6.2553e-09, m 4.653417e+35
column floor 1.0e-13 = decade(4e-06 x rho_ph 1.2705e-08) (rule)
column (mesa): 1898 rows, 37 s; rho_ph 1.270546e-08, tau-3 join 1.741746e+12 (1.00523 R_ph), floor from 1.750286e+12 (1.0102 R_ph)
startup re_zone hopf: he_ic_rad_taua 1, he_ic_rad_taub 3 (tauj): the code's startup re-solves RE only at tau_R < 3 (from the mesh top) and keeps the column's T, scaffold and (below r_b) rho; this column: tau_R(R_ph 1.741982e+12) 14.14, tau_R(r_j 1.741746e+12) 18.48
taper: off
seam: off (no seam.deck)
  r/R_ph   r[Rsun]     T[K]       rho       kappa  nabla   n_ad   n_rad  FMLT/F Gam_L  Gam_rad  Hp[Rsun]  tau
  0.9799    24.537 3.4708e+05 4.3224e-07  0.608  0.254  0.253  0.507  0.500  1.355  0.678    0.171 4.785e+03
  1.0000    25.039 5.9245e+04 1.2705e-08  2.435  5.349  0.282  5.349  0.000  5.428  5.428    0.014 1.418e+01
  Fe opacity peak: r 24.8478 Rsun (0.9924 R_ph) T 1.8340e+05 rho 7.1737e-08 kappa 1.5348, convective
  Fe CZ (F_MLT/F > 0.01): 20.2909-25.0324 Rsun (0.8104-0.9997 R_ph)
  tau at x1min 7.7354e+06; tau(R_ph) of the column 14.1782 (2/3 = the shot); R_ph 1.741982e+12
  FLAG: super-Eddington Gamma_L max 5.43 below R_ph (MLT++ carries the excess)
  R_ph ours 1.741982e+12 vs MESA 1.732691e+12: ratio 1.00536
  Teff ours 58911.79 vs MESA 59069.54: ratio 0.99733; T_ph(tau 2/3) 59244.51
  Gamma_rad = kappa F_r/(c g), x1min..r_j: column scaffold max 0.7391 at 20.672 Rsun (0 rows > 0.9); MESA's own F_MLT/F at the same P with OUR kappa would give max 0.7466 at 20.688 Rsun (0 rows > 1: sub-Eddington)
```

- Top mode: the join line has **no "(surface)"** and there is **no "OUTWARD extension" / "extension eps" / "two-sided fit"
  line**. So in mesacol.py (83d97f8c9) the `ev_join` event fired while integrating MESA's T(P) outward. Our radii are 0.5 %
  larger than MESA's, so r_j = 1.0052 R_MESA, but the join T 76.1 kK is above the MESA top-row T (59.07 kK, log T 4.7714).
  The join is therefore inside MESA's T(P) range.
- Hopf layer result: **R_ph 1.741982e12 cm (25.04 Rsun, 1.00536 R_MESA), rho_ph 1.270546e-08**, T_ph 59244.5 K.
  The `tau(R_ph) 14.18` / `tau_R(R_ph) 14.14` values are tau from the mesh top. They include the floor layer out to 4 R_ph, so
  they are not the Hopf tau. Within the Hopf layer, 2/3 is the stopping event.
- FLAG lines: only **`FLAG: super-Eddington Gamma_L max 5.43 below R_ph (MLT++ carries the excess)`**. Its maximum is
  the R_ph row itself (Gam_L = Gam_rad = 5.428, F_MLT/F 0, see the table). There is no branch FLAG. Below r_j the scaffold
  Gamma_rad stays sub-Eddington (max 0.739).

## 3. seam_only(_B).py and the passes

Script: `/resnick/groups/carnegie_poc/jingze/r136a1_1104/tools/seam_only_B.py` (23 lines, copy `R136A1_IC/seam_only_B.py`;
run by `ic/runBB.sh`). What it does:
- It does **not change the seam algorithm**. It imports `icstd.seam` from `/resnick/groups/carnegie_poc/jingze/wt_prod1102/tools/ic_star`
  (seam.py last commit a9c5797e6, 10-17) and calls `SM.fixed_point(cfg, None, None, L, None, x1max, hdr,
  AB/ic_r136a1B_noseam.txt, log, None, 30.0)`, the same loop make_star_ic.taper_and_seam uses. `lo = None` means taper off, so every
  new column takes the dump's own scaffold rows with `F_r = F (1 - F_MLT/F)`.
- It runs this outside make_star_ic only so that it can choose the code side: the deck `decks/r136a1_128_B_A.athinput` (production radial
  grid leanB 831 cells, x1min 19.45 Rsun, x1max 100 Rsun) with the overrides 1-D (nx2 = nx3 = 4), `meshblock/nx1=831`,
  `he_ic_rad_taua=1`, `he_ic_rad_taub=3`, `he_nfine=262144`, `he_gm_column=true`, opacity/Planck tables and he_gm from the
  pass-A deck, `rad_m1/implicit_flux_x1min=1.131943656797e15` (= L/(4 pi x1min^2), pure radiative wall); binary
  `builds_rstfixx/bin/athena_cpu_he_star_m1_bsgsp_3f15a44c`, omp 4, `tol_cells 0.1`, `maxit 8`.

Seam log `r136a1_1104/ic/seamBB.log` (copy `R136A1_IC/seamBB.log`), verbatim:
```
    code pass: .../ic/BB/pass_0/run.log; HSE last max |dln rho| 3.31714e-07, Gamma-clamped fine cells 1336
  seam pass 0: R_ph(column) 1.74358000e+12 -> R_ph(sph1d) 1.76639000e+12: dR 2.2810e+10 cm = 29.1259 cells (dr 7.8315e+08)
    code pass: .../ic/BB/pass_1/run.log; HSE last max |dln rho| 3.96084e-07, Gamma-clamped fine cells 817
  seam pass 1: R_ph(column) 1.76638000e+12 -> R_ph(sph1d) 1.76653000e+12: dR 1.5000e+08 cm = 0.1915 cells (dr 7.8315e+08)
    code pass: .../ic/BB/pass_2/run.log; HSE last max |dln rho| 2.27445e-07, Gamma-clamped fine cells 824
  seam pass 2: R_ph(column) 1.76653000e+12 -> R_ph(sph1d) 1.76662000e+12: dR 9.0000e+07 cm = 0.1149 cells (dr 7.8315e+08)
    code pass: .../ic/BB/pass_3/run.log; HSE last max |dln rho| 1.61754e-07, Gamma-clamped fine cells 828
  seam pass 3: R_ph(column) 1.76661000e+12 -> R_ph(sph1d) 1.76667000e+12: dR 6.0000e+07 cm = 0.0766 cells (dr 7.8315e+08)
SEAM_DONE .../ic/BB/col_3.txt [[0, 1743580000000.0, 1766390000000.0, 783152391.2299805], [1, 1766380000000.0, 1766530000000.0, 783152391.2299805], [2, 1766530000000.0, 1766620000000.0, 783152391.2299805], [3, 1766610000000.0, 1766670000000.0, 783152391.2299805]]
RC=0
```
How the passes went:
- In pass 0 the code's own startup (sph1d + re_atm + hse) moved the photosphere **outward by 2.28e10 cm = 29 cells = 0.33 Rsun**,
  from 25.06 to 25.39 Rsun. The column's R_ph is measured on the code grid here, so it reads 1.74358e12 and not 1.741982e12.
- The loop then converged quickly: 0.19 cells after pass 1, 0.115 after pass 2, and 0.077 cells (< tol 0.1) after pass 3.
- `fixed_point` returns the column that went into the converged pass. That is col_3, built from **dump_2**, which matches the
  "seam pass 2: rows = dump_2" header line. Final seam R_ph(sph1d) = **1.76667e12 cm = 25.39 Rsun**.
- The Gamma-clamped fine-cell count is 1336 in pass 0 and 817-828 afterwards. These are fine cells in which the code's startup
  HSE clamps Gamma (inference: the super-Eddington photospheric layer). The final HSE residual is |dln rho| <= 4e-7.
- The earlier first-look pass (`ic/seamB.log`, seam_only.py, 512-cell grid) behaved the same way: 30.1 cells outward in pass 0, then converged.

## 4. Gamma along the Hopf layer

Facts:
- The Hopf layer runs from r_j = 1.74174603e12 (tau 3, T 76094 K, rho 6.2553e-9) to R_ph = 1.741982e12 (tau 2/3, T 59244 K,
  rho 1.2705e-8). It is **2.36e8 cm = 0.0034 Rsun = 1.4e-4 R_ph** thick.
- Density **rises outward** by a factor 2.03 across it, i.e. there is a density inversion.
- The column has no sampled row inside the layer: drn_atm 0.0002 R_ph = 3.5e8 cm is thicker than the layer, so only the R_ph row was written.
- The make_star_ic log does not record Gamma at r_j.
- Gamma at R_ph is in the log: Gam_rad = 5.428, kappa 2.435.

Computed now (small awk over the Rosseland table above; bilinear interpolation in log kappa; it reproduces the log's kappa(R_ph) 2.4355):
- kappa(r_j: T 76094 K, rho 6.255e-9) = 0.854, so Gamma(r_j) = kappa F/(c g) = 0.854 x 6.832e14 / (3e10 x 1.024e4) = **1.90**.
- kappa(R_ph) = 2.435, so **Gamma(R_ph) = 5.43**. kappa rises monotonically as T falls toward 59 kK, so Gamma is > 1 throughout the layer
  (inference from the endpoints, since no interior rows exist).
- So R136 **did** pass the outward Hopf layer to tau 2/3 with Gamma 1.9 -> 5.4, at r_j = 25.04 Rsun ~ 1.005 R_MESA.

Inference, not tested, for the T9 comparison:
- The Hopf integration (`hopf_out`: gas HSE with the radiation force of F, outward from (r_j, rho_j, tau_j) to tau 2/3) survives
  Gamma > 1 here because rho_j ~ 6e-9 is high. The kappa rho dr needed for tau 3 -> 2/3 is used up in a few 1e8 cm, before the
  inversion can run away. The density only doubles.
- T9 fails at r_j = 83 Rsun. There the density at the join is presumably much lower and the layer would have to be much thicker for the same
  Delta tau, so a Gamma > 1 inversion has room to diverge, or never reach tau 2/3 within [r_j, 2 r_j].
- R136 therefore does not show a separate "route". It is the same default hopf path, which works because the R136 photosphere is dense and the layer
  is thin. Whether T9 can go the same way depends on rho_j and kappa there. `top: surface` + `r_max_rsun` remains the
  only demonstrated alternative. We have not tested either choice on T9.

Also of interest: the 25.0-25.3 Rsun Gamma > 1 layer you mention lies between our column R_ph (25.04 Rsun) and the seam/code
photosphere (25.39 Rsun). That is the region the code's startup (re_atm + hse, Gamma-clamped cells) rebuilt in seam pass 0.

-- Caltech
