# wasp121_0925: the dhj model moved to WASP-121b

This directory holds the WASP-121b inputs for `deep_hot_jupiter_rt`, the literature behind
them, where the observed radius sits, and the setup tool that derives the keys.

- **Inputs:**
  - `sparc_w121.athinput` is a drop-in for `sparc_0925/sparc.athinput`.
  - `prod_w121.athinput` is the nx1 = 256 production grid.
  - `grid_w121.env` holds the command-line grid overrides.
- **Binary:** needs rt-integration a6eaac93 plus **branch dhj-albedo (cbdb9028)**, which adds
  `problem/albedo`. A binary without it silently ignores the key and uses the fit albedo
  (about 0.009).
- **Tool:** `inputs/production/tools/planet_setup.py` on **branch planet-setup (35d54d33)**, in
  worktree `/viper/ptmp2/jinma/wt_planet`. Nothing is merged or pushed.

## 1. System parameters

### Candidate values

Citation counts are Semantic Scholar counts on 2026-09-25.

| source (cites) | M_p [M_J] | R_p [R_J] (band) | a [AU] | Teff [K] | R_* [Rsun] | M_* | [Fe/H]_* | Teq quoted |
|---|---|---|---|---|---|---|---|---|
| Delrez+2016 MNRAS 458, 4025 (147) | 1.183 +-0.064 | 1.807 sph. / 1.865 +-0.044 Roche-corrected (TRAPPIST z', Euler r', B) | 0.02544 | 6459 +-140 | 1.458 +-0.030 | 1.353 | +0.13 | 2358 +-52 (A = 0, f = 1/4) |
| Evans+2018 AJ 156, 283 (119) | 1.18 (D16) | ~1.7 (HST white light 0.29-1.7 um) | D16 | D16 | D16 | D16 | D16 | ">2500" |
| Bourrier+2020 A&A 635, A205 (161) | 1.157 +-0.070 | 1.753 +-0.036 (TESS 0.6-1.0 um) | 0.02596 | D16 | D16 | 1.358 | - | - |
| Patel & Espinoza 2022 AJ 163, 228 (46) | - | Rp/R* 0.1217 (TESS) | a/R* 3.81 | D16 | D16 | | | |
| **Sing+2024 AJ 168, 231 (14), NASA Exoplanet Archive default** | **1.170 +-0.043** (dynamical) | **1.7420 +-0.0060** (Rp/R* 0.122551 from Mikal-Evans+2023, **JWST NIRSpec NRS1 2.7-3.7 um**) | **0.02571 +-0.00010** | **6628 +-66** | **1.461 +-0.005** | 1.330 | +0.17 | 2409 +-24 |

- Bourrier+2020 measured **P = 1.27492504 d**, and Sing+2024 uses it. The orbit is circular:
  e < 0.0078 at 3 sigma (Bourrier).
- The planet is **tidally locked**. Delrez, Bourrier, Mikal-Evans+2020 and Sing+2024 all
  assume synchronous rotation, so P_rot = P_orb.

### Adopted set: Sing+2024

It is the latest homogeneous analysis, and it is the archive default. Its M_p, R_p, a,
Teff and R_* come from one self-consistent solution.

- **Teq check:** Teff sqrt(R_*/2a) = **2409.3 K**, which reproduces their 2409 K. This is the
  zero-albedo, full-redistribution convention, the same one the code uses (T_irr = sqrt2 Teq).
- **Most-cited set (Delrez+2016):** Teq 2357.9 K, which is -51 K (-2.1 %).
  - Its R_p of 1.865 R_J is the Roche-corrected z' radius, +7.1 % over ours.
  - g(R_p) is 843 against 956 (-12 %).
  - Through the tool it gives x1min 1.1862e10 (+6.3 %) and grav 1065.1 (-10.5 %).
- **Bourrier+2020:** Teq 2334 K, x1min 1.1159e10 (-0.01 %), grav 1177.1 (-1.1 %).
- These sensitivity runs are in `param_sets.txt`. They reuse the adopted IC and each set's own
  p_ref estimate.

## 2. Metallicity

### Literature

The host star is [Fe/H]_* = +0.13 (Delrez) to +0.17 (Sing). Retrieved atmospheric
metallicities of WASP-121b:

| source | atmospheric metallicity |
|---|---|
| Evans+2018 | 20x solar (10-30x), equilibrium-chemistry grid |
| Mikal-Evans+2019 | [M/H] = 1.09 (+0.57 / -0.69) |
| Mikal-Evans+2022 | 0.76 (+0.30 / -0.62) dayside; they then fixed 5x solar |
| Smith+2024 (IGRINS) | about stellar |
| Pelletier+2025 (NIRISS eclipse) | refractories [R/H] = 1.17 +-0.2 |
| Evans-Soma+2025 | C/H about 24x, O/H about 12x stellar |
| Kahle+2026 (0.6-12 um, preprint) | C/H about 3.4x, O/H about 2.9x, refractories about 10x stellar |

C/O is 0.57 (Kahle) to 0.92 (Evans-Soma).

**The most defensible single value is about 3-10x solar.** Refractories, which carry the
optical opacity, are near 10x. Volatiles are about 3-5x.

### Adopted: 10x solar, [M/H] = +1

This is the atmospheric value, as a GCM should use. Premixed k-tables exist only at 1x,
10x, 100x and 1000x (upstream Exo-FMS_column_ck), so there is no 5x table. Every place
metallicity enters uses the same value:

- **`problem/met = 1.0`.** The EOS metal donors follow it through eos_metal_mh. `rad_met` is
  set to 1.0 too; it is unused under route B.
- **EOS composition:** `eos_xh = 0.6588` and `eos_yhe = 0.2218`, giving Z = 0.1194. This keeps
  the solar He/H and puts 10x the solar metals-to-H (Z/X = 0.1815).
  - EOS dump: `met10/dump/`, from the deepconv_0925 dump binary with these keys.
- **ck tables:** upstream `Premixed_10x_g8_11.txt` and `FastChem_ck_10x_int.txt`, plus a hiT
  extension.
  - The extension is made by the repo's `gen_hitemp.py`, unchanged, using CK_DATA pointing at
    the 10x files and FastChem with Asplund09 metals +1 dex (`met10/gen10.sh`, `gen10.log`).
  - Continuity at 6100 K: the CE ratio upstream/FastChem is 0.935-1.005 (1x: 0.94-1.006).
  - The largest band-mean k jump is 5900 to 6100 K, 0.64 in ln (1x: 0.50). Above 6100 K it is
    the same carrier taper as 1x.
  - Installed in `ckdata10/`: `ck/`, `CE_tables/`, and symlinks to the repo's `cia/`, `ray/` and
    `sw_flux/`.
- **RCE:** `rce_w121.py` uses the same 10x tables and 10x EOS.
- **Tint:** Thorngren's Tint is independent of metallicity. It is 567.37 K from Teq 2409.3 K.

### What 10x changes (global-mean RCE, T in K at p in bar)

| profile | 1e-8 | 1e-6 | 1e-4 | 1e-3 | 1e-2 | 0.1 | 1 | 10 | 100 | 250 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1x, A fit (0.0094) | 3862 | 3767 | 3266 | 2733 | 2184 | 2105 | 2125 | 2632 | 3839 | 4383 |
| 1x, A 0.277 | 3735 | 3614 | 3124 | 2312 | 2034 | 1958 | 1984 | 2533 | 3763 | 4306 |
| **10x, A 0.277 (adopted)** | 3714 | 3190 | 2750 | 2034 | 1978 | 2018 | 2210 | 3199 | 4324 | 4898 |

- At 10x the deep atmosphere is 560-670 K hotter from 10 to 250 bar.
- The upper layers are up to 420 K colder.
- The radiative-convective boundary moves from about 7.5 bar (1x) to about 1.5 bar (10x).
- p_ref moves from 7.0e-4 to 1.8e-4 bar (section 5).

## 3. Stellar SED

The code takes the SED either from `ck_star_teff > 0` (a blackbody) or from `ck_star_teff <= 0`
plus `ck_swflux` (band fluxes under `ck_data_dir/sw_flux`).

- The repo file `sw_band_flux_W121_11.txt` is the Exo-FMS WASP-121 model spectrum (Lee+2021).
  PROVENANCE marks it "examples only".
- Its shape departs from any blackbody:
  - 0.26-0.42 um carries 0.150 of the flux, against 0.208-0.220 for a 6460-6628 K
    blackbody, because of line blanketing and the Balmer jump.
  - 0.42-0.61 um carries 0.297 against 0.262.
  - The best-fit blackbody would be 6240 K.
- Changing the blackbody from 6460 to 6628 K moves each band by at most 2 %. That is far
  smaller than the 28 % model-versus-blackbody difference in the UV.
- **Adopted: `ck_star_teff = 0`, `ck_swflux = sw_band_flux_W121_11.txt`,** in both the input and
  the RCE (`RCE_TSTAR=0`). It is renormalised to sigma T_irr^4 as before. I did not regenerate
  the band fluxes from a PHOENIX spectrum.

## 4. Bond albedo

### Literature

| source | A_B | band / method |
|---|---|---|
| **Splinter+2025 (AJ 170, 323)** | **0.277 +-0.016** | JWST NIRISS/SOSS full phase curve, 0.6-2.85 um (50-83 % of the bolometric flux); T_day 2717 K, T_night 1562 K; Cowan & Agol |
| Splinter+2025, other nightside treatment | 0.307 | same data |
| Splinter+2025, energy-balance-model grid | about 0.31 | same data |
| Mikal-Evans+2022 | 0.14 +-0.08 | WFC3 1.1-1.7 um |
| Morello+2023 | 0.37 / 0.32 | Spitzer 3.6 / 4.5 um |
| Morello+2023, from TESS | 0.05 +0.18/-0.22 | TESS |
| Davenport+2025 | 0.23 / 0.002 | Spitzer; the two bands disagree |

- The brief's "Frazier+2026" value is Splinter's number. Frazier's GCMs have A_B <= 0.05 and
  over-predict the NIRISS emission, which supports a high A_B.
- **Adopted: A_B = 0.277 (Splinter+2025).** It is the only energy budget that covers most of
  the bolometric dayside and nightside flux. The systematic range is 0.25-0.33.
- The code's Parmentier+2015 fit gives 0.0094 here.

### Code (branch dhj-albedo, cbdb9028)

- `problem/albedo`, when present, replaces the get_albedo value at both call sites:
  - the RT pass (two_stream_rt.hpp: every (1-albedo) Fstar deposit, grey and ck, and the
    picket-fence Teff);
  - the picket-fence IC (`get_picket_fence_Ttau_coeff`).
- When the key is absent, behaviour and the parameter dump are unchanged.
- Gates (`gate_alb/`):
  - Unset is bitwise against a6eaac93 on both the ic_profile and the picket-fence IC paths
    (6/6 files each).
  - With albedo = 0.3, Q_sw scales by 0.7059715752 = 0.7/(1 - A_fit) in every one of 61
    cells, and the column total by the same factor. T is identical.
  - The picket-fence IC changes with the albedo, as intended.

### RCE and Teq

- The RCE uses the same (1 - 0.277) absorbed flux (`RCE_ALB`).
- Teq stays the zero-albedo literature value, 2409.3 K. Tint(Teq) is unchanged.

## 5. Where R_p sits: p_ref

### Method (`transit.py`)

- The model's own ck opacities are used: line + continuum + Rayleigh + CIA + H-, through the
  rr_lib transcription of correlated_k.hpp, with the 10x tables.
- They are evaluated on the adopted IC column.
- For each chain, the chord optical depth is computed through exact spherical shell chords.
- Band transmission is the sum over g of w_g exp(-tau).
- Per-band effective radius: R_b^2 = r0^2 + 2 int (1 - t) y dy.
- The band radii are combined as photon-weighted transit depths.
- p_ref = p(R).

### Result

- **p_ref = 1.8e-4 bar in the adopted radius band (NIRSpec NRS1, 2.7-3.7 um; ck bands 2.5-3.5
  and 3.5-4.4 um).**
  - The tau_chord = 0.56 level of the band-mean transmission gives the same value, 1.80e-4
    (Lecavelier des Etangs+2008).
  - TESS 0.6-1.0 um: 1.2e-4. WFC3 G141: 3.3e-4. Band 8 (0.61-0.85 um): 9e-5.
- **Uncertainty (limb):** the profile is a global mean. Shifting T above 0.1 bar by -300 / +300 K
  gives p_ref(NRS1) = 1.33e-4 / 2.37e-4 bar, about 0.3 H or 0.4 % in radius.
  - At 1x the same calculation gives 7.0e-4 bar (TESS 7.4e-4). Metallicity is the dominant
    uncertainty.
- **Literature:**
  - Evans+2018 set its retrieval reference radius at 1 mbar, and put the near-UV radius at
    about 20 mbar.
  - Parmentier+2018 fitted the 1 bar radius, with tau = 0.56 limb photospheres at 20-40 mbar
    at 1.4 um.
  - Gapp+2025 say transmission probes about 1e-3 bar.
  - Our p_ref is lower because at 10x metallicity the 3 um band (H2O, CO; about 10x the
    solar column opacity) and the optical band (Fe, TiO/VO, H-) become opaque higher up.
  - Our 1x value, 0.7 mbar, sits at the low end of the literature's 1-10 mbar.

### Oblateness (new with WASP-121b)

- P_rot = 1.27 d makes the rot_potential centrifugal term 7.7 times larger than for the old
  planet.
- On the p_ref isobar the equator bulges 2.2 % above the pole.
- The tool therefore sets the polar radius so that the terminator's area-equivalent radius,
  sqrt(<r^2>), equals R_p:
  - pole 1.232081e10;
  - equator 1.259071e10;
  - area-equivalent 1.245391e10 = 1.742 R_J.

## 6. Derived keys (`setup_w121.txt`)

**<problem> and <mesh>:**

| key | value | note |
|---|---|---|
| ap = x1min | 1.116018e10 cm | p = 250 bar on the polar column |
| grav | 1190.073 | g(R_p) = 955.7 |
| Teq | 2409.3 | |
| omega | 5.704026e-5 | P = 110153.5 s |
| x1max | 1.530921e10 | p = 1.43e-9 bar on the IC column |

- **x1max rule.** The old x1max, 2.0556e10, was measured on hot 3-D runs. On a cold RCE IC it
  puts about 20 % of cells on the density floor (sparc smoke job 11971416). x1max is now set
  3 % of the height inside r(1e-9 bar) of the IC column.
- **Bounds.** A column 40 % colder above 0.01 bar (a nightside) would have p(x1max) = 5e-15 bar.
  A column 30 % hotter would have 4.3e-8 bar. No single x1max keeps every column at 1e-9 bar.
  - The old measured-extent rule, scaled, would give 1.583e10.
- **Floors.**
  - dfloor 5e-14 would floor the top two decades: the IC density at x1max is 6.4e-15. It is
    lowered to **1e-16**.
  - pfloor 1e-5 barye (1e-11 bar) and tfloor 200 K are kept.
- **The ic_profile** (`ic_w121.txt`) is the final RCE (`rceC10`, iterated to the derived ap and
  grav) plus an isothermal top row at 1e-10 bar. It covers 1e-10 to 283 bar.

## 7. Radial grids (`planet_setup.py --grids`, basis = the IC column, H = p/(rho g))

| grid | nx1 | rule | min cells/H by band (1e-6..1e-3 / 1e-3..0.1 / 0.1..10 / 10..300 bar) | dr [cm] | dt est. |
|---|---|---|---|---|---|
| SPARC-like (`G_SPARC`) | **74** | >= 3/H where p > 1e-6 bar; need 72.9 | 3.07 / 3.19 / 3.23 / 3.18 | 1.84e7 .. 5.04e8 | 14.3 s |
| production (`G_PROD`) | **256** | plan A (5/H 1e-6..300 bar), 8c | 9.36 / 9.96 / 9.73 / 9.62 (f = 1.87 of plan A) | 5.8e6 .. 7.1e7 | 4.6 s |

- **dt** = CFL 0.3 min(dr/(c_s + 1 km/s)), radial only.
  - On the old planet the same estimator gives 22.7 s for the sparc 66 grid (calibrated
    estimate 25.8 s) and 8.2 s for A 8c at 256 (design256: 7.15 s).
- **Plan A at 256 is over-met.** The domain spans fewer scale heights than the old one, so
  about 137 cells would meet plan A.
- **Cold-column bound.** A 0.6x-T nightside column gets 1.8-1.9 cells/H in 1e-6..0.1 bar on the
  sparc grid.
- **The old grids on the new IC.** On the old planet's cold RCE IC, the old sparc 66 grid gives
  only 0.95 cells/H at p > 1e-6 bar. It was designed on the hot 3-D states.
- **Regression.** The template basis (the old 3-D columns mapped at equal pressure) reproduces
  both old grids exactly with the old planet. For the new planet it produced artifacts where
  the old and new inversion bases differ, so it is not used.

## 8. Using the inputs (run.sub / smoke.sub)

- **Input and grid:** `INP=/viper/ptmp2/jinma/wasp121_0925/sparc_w121.athinput`, and
  `source /viper/ptmp2/jinma/wasp121_0925/grid_w121.env`.
  - That defines G_SPARC and G_PROD, and sets G to G_SPARC.
  - Use `$G` in place of run.sub's hard-coded `G=` block.
  - It includes x1min and x1max, which the old block did not need to set.
- **Rotation period:** set `PR=1.101535e5` in place of 3.05e5, so that TROT counts rotations.
- **Output cadences** in `O=` are hard-coded for P = 3.05e5 and must scale to the new period:
  - `output1/dt=1.10154e3` (0.01 rot);
  - `output3/dt=5.50768e4` and `output4/dt=5.50768e4` (0.5 rot).
  - The file itself has 0.1 rot (hst, log), 2 rot (bin) and 0.5 rot (rst).
- **Binary:** needs dhj-albedo, merged onto whatever athena.gpu is built from.
- **Keys:** `keys_changed.txt` lists every changed key, old -> new, with its source. It covers
  34 keys, 2 of them new: `ck_swflux` and `albedo`. Route B, nq2, hiT, sponges (pressure-based)
  and every other key are verbatim.

## 9. Gates (all CPU; `gate_w121/`, `gate_alb/`, `regress_final.txt`)

- **Tool regression** (the current planet from its own parameters, R_p taken from the
  ckrce_0925 RCE radii at 1e-3, 1e-2 and 1 bar):
  - x1min within 1.0e-4, grav within 2.0e-4, x1max within 2.2e-4;
  - IC-fill emulation r(p_ref) within 4e-5;
  - with `--grids --grid_basis template`, the sparc grid comes back at nx1 66 with max
    |dc| = 5e-5, and A 8c @256 with max |dc| = 1.5e-4.
  - REGRESSION PASS.
- **Parse:** `athena -n` passes on both inputs (binary athena_cpu_albedo_cbdb9028).
- **Fresh start** of sparc_w121 on 74 x 16 x 16 per panel:
  - the counters show **no dfloor, efloor or tfloor event at cycles 0 and 1**. Cycle 2 has 140
    dfloor and 104 efloor events in the top cells once the dynamics start;
  - IC density minimum 1.33e-14 (dfloor 1e-16);
  - IC T matches the table to 1.7 %.
  - For comparison, the old sparc input on the same mesh has 51840 dfloor and 6912 tfloor
    events at cycle 0.
  - prod_w121 runs 1 cycle cleanly.
- **IC fill** (hydro_w bin at t = 0, p from rho and e via the 10x EOS grid):
  - over 1536 columns, r(p_ref) ranges from 1.232375e10 to 1.259518e10;
  - the tool's pole and equator are 1.232081e10 and 1.259071e10, so the error is +2.4e-4 and
    +3.6e-4;
  - the terminator mean equals R_p to better than 1e-3. The gate required < 1 %.

## 10. Caveats

- **Profile extended (09-25 follow-up).** `ic_w121.txt` now has 269 uniform rows
  (0.0499 dex) from 9.5e-11 to 2240 bar: isothermal at the top, and the exact 10x-EOS
  adiabat below 283 bar. The RCE bottom slope is 0.1411 against nabla_ad 0.1407.
  - Fresh start (`gate_w121/run_gz`, ghost zones written): no floor events at cycles 0-1.
  - First active cells: the pole-like cell sits exactly on the table (215 bar, 4797 K). The
    equator-like cell is 631 bar at 5068 K, against 5603 K on the adiabat (-9.6 %).
  - Ghosts reach 1235 bar at 5182 K against 6229 K (-17 %).
  - Cause: the code, not the table. get_init_eos fills every cell below z_eff = 0 with the
    isothermal T(250 bar) and an exponential p, so the table is never read past its end.
  - Putting those cells on the adiabat needs a small pgen change: start get_init_eos_arr
    below z = 0. It is not made.
- **Inner wall at the equator.** The rot_potential bulge puts the equatorial inner-wall cells
  at up to about 630 bar (the pole is at 250 bar). Below z = 0 the IC is the pgen's isothermal
  extension. The bottom sponge (50-100 bar) is pressure-based and unaffected.
- **ck convergence not tested.** The 10x hiT tables were not run on GPU. The ck-implicit
  convergence from this fresh start was not tested: the smoke job saw 14-19 of 75
  NOT-CONVERGED calls on the old planet.
- **The IC is a global mean.** The deep (4300-4900 K at 100-250 bar) sits on the 1-D RCE
  adiabat.

## Files

- `rce_w121.py`: rce.py parametrised by environment variables: RCE_TEQ, RCE_G0, RCE_AP,
  RCE_TSTAR (0 = SED file), RCE_MET, RCE_CKDIR, RCE_EOSDIR, RCE_ALB, RCE_RESTOL.
  - Its regression against ckrce_nq2_B is 6.5e-5 (`rce_regr/`).
- `rr_lib_met.py`: rr_lib with switches for the metallicity tables and the SED file.
- `transit.py`; results in `transit_C10.txt`.
- `make_inputs.py`, `keys_changed.txt`.
- `rce*/`: RCE runs. `rceC10` is the adopted one; rce1, rceA1, rceA10 and rceB10 are the
  comparison and iteration runs.
- `met10/`: 10x tables, FastChem input and the EOS dump.
- `ckdata10/`: the ck_data_dir.
- `setup_w121.{txt,npz}`, `grid_w121.env`, `param_sets.txt`.
- `athena_cpu_a6eaac93`, `athena_cpu_albedo_cbdb9028`: the gate binaries. The build trees are
  deleted.

## 11. Follow-up 09-25 (afternoon)

### IC column below the anchor

**Branch dhj-ic-deep 23b3fb6a** (from rt-integration 038148c6). When `ic_profile` is set, the
IC column is integrated down from 250 bar to 1.1 x the lowest z_eff of any cell or inner
ghost.

- The WASP-121b column reaches z = -2.75e8 cm.
- Table lookups are clamped, so nothing is read past the table's end.
- Gates (`gate_icdeep/`):
  - The picket-fence path is bitwise (6/6 files).
  - Deepest cells and inner ghosts are on the table adiabat to < 1e-4:
    - active cells: pole 215 bar / 4797 K, equator 593 bar / 5552 K;
    - inner ghosts up to 1037 bar / 6055 K.
  - No floor events at cycles 0-1.
  - r(p_ref) at the tool's pole and equator to 2.4e-4 and 3.7e-4.
- The CPU binary is `athena_cpu_icdeep_23b3fb6a`.

### sparc_0925 binary and scripts

- `athena.gpu` = 23b3fb6a, md5 e2d1a403 (the previous binary is kept as
  `athena.gpu.a6eaac93`; details in BUILD_COMMIT.txt).
- `run.sub` and `smoke.sub` now default to `PLANET=w121`:
  - INP is `sparc_w121.athinput`;
  - G comes from `GRIDENV` (default `grid_w121.env`);
  - PR = 1.101535e5;
  - output dt 1.10154e3 / 5.50768e4 / 5.50768e4.
- `PLANET=old` restores the previous planet.
- `DRY=1` prints each arm's command line instead of running it.
- The pre-edit copies are `*.pre_w121`.

### Top variants

The variants lower the model top. Both use the same >= 3 cells/H rule below 1e-6 bar and the
same keys otherwise.

| variant | p_top | x1max [cm] | p(x1max) [bar] | nx1 | dt est. | files |
|---|---|---|---|---|---|---|
| baseline | 1e-9 | 1.530921e10 | 1.43e-9 | 74 | 14.3 s | sparc_w121.athinput, grid_w121.env |
| top8 | 1e-8 | 1.455079e10 | 1.39e-8 | 72 | 14.5 s | sparc_w121_top8.athinput, grid_w121_top8.env |
| top7 | 1e-7 | 1.387406e10 | 1.35e-7 | 72 | 14.0 s | sparc_w121_top7.athinput, grid_w121_top7.env |

- Both variants parse. CPU fresh start: no floor events through cycle 2.
- Lowering the top saves only 2 cells: the free top uses few, coarse cells.

### Top sponge in each variant

The top sponge is fixed in pressure. fdrag ramps linearly in ln p from 0 at 1e-6 bar to 1 at
1e-7 bar, and the rate is fdrag/1000 s.

- **Baseline:** a full-strength layer 1.9 dex thick (1e-7 to 1.4e-9 bar) sits above the ramp.
- **top8:** the full-strength layer is 0.85 dex (1e-7 to 1.4e-8 bar).
- **top7:** there is no full-strength layer. The whole top decade is ramp, reaching only
  fdrag = 0.87 at the lid.
- Moving the window with the top (to 10-100 x p_top) would put the top7 sponge at 1.35e-6 to
  1.35e-5 bar, inside the p > 1e-6 bar data region. So the window should stay at 1e-7 to
  1e-6 bar.
- top7 therefore tests a lid with only a partial sponge. top8 keeps a full-strength layer
  about one decade thick.

## ck Newton settings for the 10x inputs (09-25, user-approved)
- sparc/prod inputs now set:
  - `ck_impl_rsec = 20`
  - `ck_impl_maxit = 12`
  - `ck_impl_tol = 1.0e-7`
  - `ck_impl_floorbound = true`
  - `ck_impl_kkt_demax = true`
- Needs rt-integration >= 6e897934 (ck-newton10x merge).
- **tol 1e-7 is a USER-APPROVED EXCEPTION for 10x metallicity.**
  - Why: the 10x day-side upper layers leave a few columns converging linearly. With rsec 20 + maxit 12 the measured
    max residual is 1.6e-7 and the median is 1e-8 (wasp121_0925/SMOKE_DIAG.md section 6).
  - The default and the 1x inputs keep 1e-8.
- Backups of the previous inputs: *.pre_rsec.
