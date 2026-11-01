# wasp178b: hydro-first deep-hot-Jupiter setup for WASP-178b = KELT-26b (dhj_planets_1101)

WASP-178b sits mid-trend in Seidel et al. 2026 (Nature Astronomy, DOI 10.1038/s41550-026-02870-1):
Teq 2398 +- 77 K (their zero-albedo recomputation, ED Table 1, ref. [73] = Damasceno+2024),
net terminator Fe I blueshift **4.79 +- 1.19 km/s** (2 ESPRESSO nights; ../seidel2026.md).
Host: the hottest star of the sample (A1, Teff 9200-9360 K), so the strongest short-UV flux.
This directory is the setup package for a Caltech hydro run with the WASP-121b hydro fresh
deck (`/viper/ptmp2/jinma/w121hfresh_1101/deck/hfresh.athinput`, fork/fin51-r-fix 182477553;
hydro, lhllc, `<hydro>/am_vector = true`, `problem/am_transport = false`, C32 x nx1 112).
Pipeline = the TOI-1518b worker's chain (`../toi1518b`, read-only) = the w121prod_0929/met3
chain: face-0 route-B ck RCE -> face-5 exact RCE -> TESS transit radius -> IC ->
planet_setup.py -> W121-cells/H radial grid -> deck.

## 0. What to run

- **Nominal deck (3x solar): `deck/wasp178b_hfresh.athinput`**; **1x arm: `deck/wasp178b_hfresh_1x_arm.athinput`**
  (the uniform-1x trend set; Caltech round 1 runs the 1x decks of all planets first). md5s in `MD5SUMS`.
  - Every changed key is marked `# wasp178b:` with the WASP-121b value; made by `scripts/make_deck.py
    deck/ovr_3x.json` / `deck/ovr_1x_arm.json` from `/viper/ptmp2/jinma/w121hfresh_1101/deck/hfresh.athinput`.
  - Everything else is hfresh verbatim: hydro, `rsolver = lhllc` (cb_sensor default), `<hydro>/am_vector = true`,
    `problem/am_transport = false` (required), cs_seam keys, top sponge 1e-7..1e-6 bar, no bottom sponge,
    ck T4/C2/E4, dfloor 1e-18 / pfloor 1e-7 / tfloor 200 K, eos_logd_min -19.
  - 3x nominal changes: grid (nx1 120, x1min 1.178972e10, x1max 1.898531e10, 8 coefficients + 2 bumps),
    Teq 2452.3, omega 2.174158e-05, grav 1512.9734, ap 1.178972e10, albedo 0.15, ck_swflux
    sw_band_flux_W178_11.txt, ic_profile ic/ic_wasp178b_3x_f5.txt, tlim 120 rot, output cadences in rotations
    (hst 0.01, log 0.1, bin 0.5, rst 1 rot = 288993.9 s); metallicity keys = hfresh (3x).
  - 1x arm: nx1 120, x1min 1.178808e10, x1max 1.932467e10, grav 1513.3951, ap 1.178808e10,
    ic/ic_wasp178b_1x_f5_arm.txt, and the 1x keys (eos_xh 0.7381, eos_yhe 0.2485, met 0.0, rad_met 0.0,
    ck_data_dir ck/CK1X, ck_table ck/CK1X/ck/Premixed_1x_g8_11_hiT2.txt, ck_ce_table
    CE_tables/FastChem_ck_1x_int_hiT2.txt).
- **Run directory:** `ck` and `ic` (copies of or symlinks to `ck/`, `ic/` here); relative table paths; no `-d`.
  `ck/CK3X`, `ck/CK1X` (and unused `ck/CK10X`) = symlinks to the WASP-121b package / repo tables plus a real
  `sw_flux/` with **`sw_band_flux_W178_11.txt`** (new): Caltech copies that file into its CK3X and CK1X
  `sw_flux/`.
- **Commands:**
  ```bash
  cd <rundir> || exit 1          # ck -> <pkg>/ck, ic -> <pkg>/ic
  srun -n 2 <athena> -i <pkg>/deck/wasp178b_hfresh.athinput time/nlim=200 > smoke.log 2>&1
  srun -n 2 <athena> -i <pkg>/deck/wasp178b_hfresh.athinput -t <hh:mm:ss> > run.log 2>&1
  ```
- **Binary:** any fork/fin51-r-fix >= 182477553 with PROBLEM=deep_hot_jupiter_rt. Roche tide keys: section 6.

## 1. System parameters

### Candidate sets

| source | Teff [K] | R_* [R_sun] | M_* [M_sun] | a/R_* | R_p [R_J] | M_p [M_J] | [Fe/H]_* | Teff/sqrt(2 a/R_*) |
|---|---|---|---|---|---|---|---|---|
| Hellier+2019 MNRAS 490, 1479 (discovery; NASA archive default) | 9350 +- 150 | 1.67 | 2.07 | 7.17 +- 0.21 | 1.81 +- 0.09 | 1.66 +- 0.12 (RV) | +0.21 +- 0.16 | 2469 |
| Rodriguez Martinez+2020 AJ 160, 111 (KELT-26b) | 8640 +500/-240 (SED) | - | 1.93 | - | 1.94 +0.060/-0.058 | 1.41 +0.43/-0.51 | - | - |
| Pagano+2024 A&A 682, A102 (CHEOPS) | (9350) | 1.772 +- 0.020 | 2.169 | - | - | - | - | - |
| Damasceno+2024 A&A 689, A54 (ESPRESSO; Seidel's ref. [73]) | 9350 (Hellier) | 1.772 (Pagano) | 2.169 (Pagano) | 7.00 +0.03/-0.13 | 1.796 +- 0.022 | 1.66 (Hellier) | 0.21 (Hellier) | 2499 |
| **Lothringer+2025 AJ 169, 274 (HST+JWST+TESS, new stellar characterisation)** | **9200 +200/-170** | **1.70 +- 0.01** | **2.02 +0.04/-0.05** | **7.037 +- 0.015** | **1.77 +- 0.01** | 1.66 (Hellier) | +0.03 +0.10/-0.12 | **2452** |

### Adopted: Lothringer+2025 (star, a/R_*, R_p), Hellier+2019 (M_p), Damasceno+2024 (P)

- It is the newest homogeneous solution (star re-derived, joint transit fit), and it is
  self-consistent: R_p/R_* = 1.77 R_J / 1.70 R_sun = 0.10699 matches the measured transit depth
  (Damasceno R_p/R_* = 0.10716), and Kepler a(P, M_* + M_p) = 0.05535 AU vs a/R_* R_* = 0.055633 AU
  (0.5 %). Damasceno's set mixes sources: its R_p/R_* x R_* gives 1.848 R_J, not the quoted
  1.796 R_J.
- P = 3.3448371 d (Damasceno+2024 ephemeris); circular (e fixed 0 everywhere); tidally locked.
- M_p = 1.66 +- 0.12 M_J (Hellier RV, 7 %; adopted by all later papers). With Lothringer's
  M_* 2.02 instead of Hellier's 2.07 the RV mass would be 1.63 M_J (-1.6 %, ignored).
- g(R_p) = 1313.3 cm/s^2 (log g 3.118; Damasceno 3.106 +- 0.033).
- **Teq (code) = 9200/sqrt(2 x 7.037) = 2452.3 K.** Seidel's 2398 K (cited from Damasceno+2024)
  is NOT reproduced by Teff/sqrt(2 a/R_*) of any published set (Damasceno gives 2499 K,
  Hellier 2469 K, Lothringer 2452 K); 2398 K would need a/R_* = 7.4 at 9200 K or Teff 9000 K at
  7.037. The model uses 2452.3 K (+2.3 % vs the paper, the same sign as WASP-121b, +2.2 %, and
  TOI-1518b, +2.2 %, in the code convention). **Flag:** Teff is the dominant uncertainty: the
  KELT SED value 8640 K would give Teq 2303 K.
- omega = 2.174158e-05 s^-1 (v_eq(R_p) = 2.75 km/s; slowest rotator of the sample, P = 2.62x WASP-121b's).
- Stellar v sin i 10.5 km/s (Damasceno) with a polar orbit (lambda 105.7 deg): irrelevant for the
  synchronous planet model.

## 2. Metallicity (own section; user decision 11-01: per-planet nominal from the evidence)

**Decision (user, 11-01): WASP-178b = 3x solar NOMINAL, with a 1x arm** (the uniform-1x trend set
used for the other Seidel planets). Both decks are built and CPU-gated.

Evidence ranked as asked: (1) the planet's own retrieved metallicity; (2) mass-metallicity
expectation; (3) host [Fe/H] as context only.

### (1) Retrieved atmospheric metallicity of WASP-178b

| source | data | result |
|---|---|---|
| Saha & Jenkins 2025 (arXiv:2510.11479) | reanalysis of JWST NIRSpec G395H (GO-2055) transmission | **11.44 +12.54/-6.94 x solar**, super-solar C/O (abstract; full paper not read) |
| Lothringer+2025 AJ 169, 274 (arXiv:2503.15472), Table 3 | HST UVIS/WFC3 + JWST NIRSpec | evidence-favoured pRT JWST+HST chem.-eq. 3TP model (weight 0.96): **[O/H] = +1.49 +inf/-1.01** (>= 3x solar), **[C/H] = -2.04 +- 0.2**, **[Si/H] = -0.13 +- 0.13** (0.74x), C/O ~0.01; ATMO JWST-only [O/H] +1.64, [C/H] +1.47, [Si/H] +1.24; masked-CO2 test [Z/H] = +1.61 +- 0.28. Abstract: "solar to super-solar [O/H] and [Si/H], sub-solar [C/H]". |
| Cont+2024 A&A 688, A206 (arXiv:2406.08166) | CRIRES+ K-band dayside + TESS/CHEOPS eclipses | **metallicity 1.4 +- 1.6 dex** ("solar to super-solar"), C/O 0.6 +- 0.2, A_g < 0.23 |
| Lothringer+2022 Nature 604, 49 | HST NUV transmission | SiO/Mg/Fe/Fe+ absorption; "solar or slightly super-solar" (as summarised in arXiv:2506.01132) |
| Damasceno+2024 A&A 689, A54 | ESPRESSO | Na, H-alpha/beta, Mg, Fe I, Fe II detections; no bulk value |

All three bulk/volatile constraints sit at or above solar (central values 11x, 25x, >=3x for O),
with ~1 dex errors; Si is ~solar and C is sub-solar in the Lothringer model (but C/O is solar in
emission and super-solar in Saha's reanalysis). Uniformly scaled tables cannot represent the
element pattern.

### (2) Mass-metallicity expectation (1.66 M_J)
- Thorngren+2016 (ApJ 831, 64) BULK: Z_p/Z_* = 9.7 M^-0.45 = 7.8x stellar = ~8x solar: an upper
  bound on the atmosphere (Thorngren & Fortney 2019).
- Welbanks+2019 (ApJL 887, L20): hot-Jupiter H2O generally sub-stellar (scatter ~1 dex).

### (3) Host star
[Fe/H]_* = +0.03 +0.10/-0.12 (Lothringer+2025; Hellier+2019: +0.21 +- 0.16).

### Choice and how it compares with 1x
- **Nominal 3x solar ([M/H] +0.477), 0.45 dex above the host.** It is the lowest value
  consistent with the volatile constraints (O >= 3x; Saha 11x, -0.4 dex at 1 sigma ~4.5x),
  below the bulk upper bound (~8x), and equal to the WASP-121b production metallicity.
- **1x arm**: what the uniform-1x trend set uses (Seidel 2026 assumes ~solar/stellar and Gandhi+2023
  find UHJ Fe ~ stellar); it is 1.3-1.5 sigma below the JWST bulk values but matches the
  Lothringer Si.
- A 10x arm (CK10X tables exist) was started and stopped after the decision (no deck).

| arm | ck table (relative to the run dir) | CE table (relative to ck_data_dir) | EOS dump | deck EOS keys | available? |
|---|---|---|---|---|---|
| **3x (NOMINAL)** | `ck/CK3X/ck/Premixed_3x_g8_11_hiT2.txt` (-> w121_s10fresh_1027/pkg) | `CE_tables/FastChem_ck_3x_int_hiT2.txt` | `/viper/ptmp2/jinma/w121prod_0929/met3/dump` | as hfresh (eos_xh 0.7188, eos_yhe 0.2420, met 0.4771, rad_met 0.4771) | yes |
| 1x (arm) | `ck/CK1X/ck/Premixed_1x_g8_11_hiT2.txt` (-> repo data/exo_fms_ck) | `CE_tables/FastChem_ck_1x_int_hiT2.txt` | `/viper/ptmp2/jinma/deepconv_0925/dump` | eos_xh 0.7381, eos_yhe 0.2485, met 0, rad_met 0 | yes |
| 10x | `ck/CK10X/...` (-> wasp121_0925/ckdata10) | 10x | `wasp121_0925/met10/dump` | 0.6588, 0.2218, met 1 | yes (unused) |
| 30x / 100x; non-scaled C-poor/O-rich | - | - | - | - | MISSING (generation + EOS dump needed) |

## 3. Host SED, 11-band stellar flux, Bond albedo, T_int

### SED -> `sw_band_flux_W178_11.txt`
- `sed/make_sed.py` (TOI-1518b script generalised to [M/H] in [-0.5, +0.5], ck04models
  ckp00/ckp02): ATLAS9 surface fluxes trilinear in Teff 9200 K (9000/9250 grids), log g 4.28,
  [M/H] +0.03 (Lothringer+2025 star). Same 11 band edges as `sw_band_flux_W121_11.txt`; < 0.26 um
  folded into band 10 (0.26-0.42 um), > 324.68 um into the last; ascending. Integral 0.9982
  sigma Teff^4. Command: `python3 make_sed.py 9200 4.28 sw_band_flux_W178_11.txt 0.03`
  (`sed/make_sed.out`). Superseded Damasceno-star version (9350 K, +0.21) in
  `sed/old_damasceno/` (below-0.26 fraction 17.2 %).
- `ck_star_teff = 0`, `ck_swflux = sw_band_flux_W178_11.txt` (in `ck/CK{1X,3X,10X}/sw_flux/`).

| band [um] | ATLAS9 9200 K | blackbody | model/BB |
|---|---|---|---|
| 0.26-0.42 | 0.4092 | 0.4556 | 0.898 |
| 0.42-0.61 | 0.3101 | 0.2466 | 1.257 |
| 0.61-0.85 | 0.1415 | 0.1466 | 0.965 |
| 0.85-1.32 | 0.0923 | 0.0978 | 0.943 |
| 1.32-2.02 | 0.0320 | 0.0357 | 0.895 |
| 2.02-2.50 | 0.0067 | 0.0078 | 0.870 |
| 2.50-3.50 | 0.0051 | 0.0060 | 0.845 |
| 3.50-4.40 | 0.0015 | 0.0018 | 0.828 |
| 4.40-8.70 | 0.0014 | 0.0017 | 0.814 |
| 8.70-20.00 | 0.0003 | 0.0003 | 1.192 |
| 20.00-324.68 | 0.0000 | 0.0000 | 1.341 |

- **Fraction of the stellar flux below 0.26 um folded into band 10: 16.5 %** (blackbody 13.9 %;
  the A-star Balmer jump moves flux from the near UV into 0.42-0.61 um, but the far UV is above a
  blackbody). Band 10 then carries 0.41 of the absorbed stellar power. More than KELT-20
  (14.4 %) and TOI-1518 (4.3 %); the hardest host of the sample.
- **Where band 10 deposits** (kelt20b_work/uvdep.py on the RCE column, `rce/uvdep.txt`, 3x nominal column): band 10 is 10 % deposited by 1.2e-4 bar, 50 % by 3.5e-3 bar,
  90 % by 5.3e-2 bar (all bands: 50 % at 2.1e-3 bar).
  Real < 0.26 um flux (Fe II, Mg, SiO continua, Lothringer+2022) would be absorbed far higher
  (~1e-5-1e-6 bar or above). With 16.5 % of the stellar flux deposited ~2 decades too deep, the
  model's upper atmosphere (1e-6-1e-4 bar) is likely too cool and the 1e-3-1e-2 bar layer too
  warm. **Flag for Caltech:** the KELT-20b worker built a UV pseudo-band test
  (kelt20b_work/rce_f5_uv.py); the same test would apply here, more strongly.

### Bond albedo: A_B = 0.15 (range 0-0.35)
- Pagano+2024 (A&A 682, A102; CHEOPS): geometric albedo 0.1 < A_g < 0.35, recirculation
  eps > 0.7, dayside brightness 2250-2800 K. Cont+2024: A_g < 0.23. Lambertian A_B = 1.5 A_g gives
  < 0.35 from Cont's limit and >= 0.15 from Pagano's lower bound (the CHEOPS band of a 2500 K
  planet still has thermal emission, so Pagano's lower bound is soft). Adopted 0.15.
- `problem/albedo = 0.15`; RCE_ALB = 0.15.

### T_int = 553.4 K
- Thorngren+2019 fit from Teq (as rce_w121.py and `two_stream_rt::get_Tint`): TINT 553.42 in the
  RCE log. WASP-121b: 567.4 K.

## 4. RCE, transit radius, IC (`rce/`, `ic/`)

Same chain and solvers as WASP-121b 3x / TOI-1518b (`chain.sh`, `finish_f5.sh`; read-only
`wasp121_0925/rce_w121.py`, `w121fresh_1001/ic/rce_f5.py` (F5_INLINE=1), `wasp121_0925/transit.py`,
`scripts/make_ic.py`, `scripts/planet_setup.py` = fork/fin51-r-fix planet_setup.py). Env
`scripts/rce_env.sh` (3x: `ck/CK3X`, EOS dump w121prod_0929/met3, RCE_TEQ 2452.3, RCE_ALB 0.15,
RCE_TSTAR 0 + `sw_band_flux_W178_11.txt`). The first parameter set (Damasceno star, Teq 2498.9)
was abandoned after a few iterations (`rce/old_damasceno_rceA`).

| step | dir | start | G0 / AP | result |
|---|---|---|---|---|
| face-0 pass A | `rce/r3x_A` | W121 3x RCE (w121prod_0929/met3/rceB3f) | 1510 / 1.18e10 (guess) | 40 it (first launch killed at ~05:22 by another worker's pkill, restarted from its checkpoint); keys ap 1.18123e10, grav 1507.195 |
| face-0 pass B | `rce/r3x_B` | pass A | 1507.195 / 1.18123e10 | plateau res 4.2e-6 (it 146), stopped; 2-it pass writes ckrce |
| face-5 | `rce/r3x_B/f5` | face-0 B | same | it 40 res 4.6e-4 (it 20 -> 40: < 10 K change at every level); finalized from `ckpt_it40.npz` (1 it). **Energy closure 0.9776**; outer max|dT/T| 9.5e-5 |
| transit / setup | `rce/r3x_B/{transit,setup}.txt` | | | **p_ref (TESS) = 6.651e-04 bar** (dT -300 K: 6.28e-4, +300 K: 7.41e-4) |

- **Keys: ap = x1min = 1.178972e10 cm, grav = 1512.9734 cm/s^2** (g(R_p) 1313.34); p_ref isobar
  pole 1.263966e10, equator 1.266855e10, area-equivalent 1.265408e10 = R_p. The RCE ran at
  1507.195 / 1.18123e10 (0.4 % / 0.2 % off; no third pass).
- IC `ic/ic_wasp178b_3x_f5.txt` (269 rows, 9.49e-11..2240.7 bar, 6158 K at the bottom). Column:

| p [bar] | r [cm] | T [K] | H [cm] |
|---|---|---|---|
| 100 | 1.18842e10 | 3856 | 9.57e7 |
| 10 | 1.20677e10 | 2658 | 6.47e7 |
| 1 | 1.21962e10 | 2141 | 5.29e7 |
| 0.1 | 1.23200e10 | 2173 | 5.51e7 |
| 0.01 | 1.24503e10 | 2205 | 5.82e7 |
| 1e-3 | 1.26015e10 | 2476 | 8.55e7 |
| 1e-4 | 1.29582e10 | 3596 | 1.852e8 |
| 1e-6 | 1.39581e10 | 4192 | 2.506e8 |
| 1e-9 | 1.59530e10 | 4239 | 3.325e8 |

- 1x arm (`arms/r1x`, `arms/r1x_C`): see section 5b.

## 5. Radial grid, domain top, sponge (`grid/nominal/`)

- `scripts/grid112.py` (W121 hfresh cells/H at equal pressure, 8 coefficients + 2 bumps, minimax):
  `PLANET_EOSD=/viper/ptmp2/jinma/w121prod_0929/met3/dump python3 scripts/grid112.py
  ic/ic_wasp178b_3x_f5.txt 1512.9734 1.178972e10 112 116 120 124` (via `post.sh`).
- x1max = r(W121 top pressure 6.12e-18 bar) on this column = **1.898531e10** (x1max/ap 1.61, W121 1.87;
  0.36 R_Hill: no Roche-lobe issue).
- max dr/target: nx1 112 1.077, 116 1.076, **120 1.0024**, 124 0.958. **Adopted nx1 = 120** (1.0024 is
  inside the W121 self-regression of the same fit, 1.0038 at 112, `../toi1518b/gate/grid112_regress_w121.txt`).
  Cells/H (min/median): 1e-6..1e-4 bar 3.26/4.27, 1e-4..1e-2 4.29/4.59, 1e-2..1 4.43/4.70, 1..10 4.04/4.13,
  10..300 4.23/4.30 (W121 2.98/3.99, 4.04/4.34, 4.10/4.39, 3.80/3.89, 4.01/4.03); 33 cells above 1e-6 bar.
- Sponge: pressure based (1e-6 -> 1e-7 bar), r 1.396e10..~1.45e10, inside x1max; unchanged.
- dt radial CFL estimate 8.41 s (W121 11.43 s, ratio 0.74).

## 5b. 1x arm (`arms/r1x`, `arms/r1x_C`, `grid/arm1x/`)

- Env `scripts/rce_env_1x.sh` (CK1X tables, deepconv_0925 1x EOS dump). Chain `run_arm.sh` +
  `arms/stop_at40.sh` + `finish_f5.sh`:
  face-0 60 it from the 3x face-0 checkpoint at the 3x pass-A keys (1507.195 / 1.18123e10) -> face-5,
  stopped at it 40 (res 6.6e-4), finalized (FIN_RESTOL 2e-3; one finalize run was killed by my own
  mis-targeted kill of the 10x arm and rerun) -> keys ap 1.179009e10 / grav 1512.8786 ->
  **pass C** (`arms/r1x_C`) face-5 at those keys: limit cycle at res ~8e-4 (it 40-49), finalized from
  its it-40 checkpoint. **Energy closure 0.941** (3x nominal 0.978; TOI-1518b 1x 0.967): the 1x column is
  less converged, as for TOI-1518b.
- **p_ref (TESS) = 1.028e-03 bar** (dT -300 K: 9.97e-4; +300 K: 1.51e-3); **ap = 1.178808e10,
  grav = 1513.3951** (0.02 % / 0.03 % from the 3x keys: metallicity barely moves the anchor, because the
  deeper 1x p_ref is compensated by the cooler upper column).
- Grid: x1max = r(6.12e-18 bar) = 1.932467e10; nx1 116 1.057, **120 1.0016**, 124 0.987 -> nx1 120.
  Cost estimate 4.0x W121 per rotation at nx1 124 (`grid/arm1x/cost_124.txt`), i.e. ~3.9x at 120.
- roche_anchor.py (1x): `--pref 1.028e-03 --profile ic/ic_wasp178b_1x_f5_arm.txt --eos
  /viper/ptmp2/jinma/deepconv_0925/dump/eos_table.txt --x1max 1.932467e10` (other inputs as 3x) ->
  ap 1.181312e10, grav 1506.9869, roche_a 8.322578e11; x1max/r_L1 0.373 (`grid/arm1x/roche_anchor.txt`).
- The 10x arm (`arms/r10x`) was stopped after the decision; no IC.

## 6. Tide and Roche potential (`grid/nominal/system.txt`, `grid/nominal/roche_anchor.txt`)

- a/R_* = 7.037 (a = 8.323e11 cm), q = 7.84e-4. Hill radius 5.322e10 cm = 4.21 R_p; Eggleton
  3.717e10 = 2.94 R_p; Roche L1 (roche_anchor.py) 5.180e10 cm.
- Tide / gravity (`scripts/system.py wasp178b <r>`):

| r | r/R_p | r/R_Hill | 2(r/a)^3/q | 3(r/a)^3/q | omega^2 r/g |
|---|---|---|---|---|---|
| x1min 1.1790e10 | 0.932 | 0.222 | 0.0072 | 0.011 | 0.0037 |
| R_p 1.2654e10 | 1.000 | 0.238 | 0.0090 | 0.013 | 0.0046 |
| r(1e-6 bar) 1.3958e10 | 1.103 | 0.262 | 0.012 | 0.018 | 0.0061 |
| x1max 1.8985e10 | 1.500 | 0.357 | 0.030 | 0.045 | 0.015 |

  **The tide is ~1 % of g in the photosphere and <= 4.5 % at the top: negligible** (a/R_* 7, the widest
  orbit of the sample). The Roche port is not needed for WASP-178b.
- **roche_anchor.py** (fin51-r-fix c8aac483c, copy in `scripts/`):
  ```
  python3 scripts/roche_anchor.py --mp 1.66 --rp 1.77 --pref 6.651e-04 --a_au 0.055633 \
      --porb 3.3448371 --teff 9200 --rstar 1.70 --mstar 2.02 --profile ic/ic_wasp178b_3x_f5.txt \
      --eos /viper/ptmp2/jinma/w121prod_0929/met3/dump/eos_table.txt --x1max 1.898531e10
  ```
  Output: **ap = x1min = 1.181477e10 cm, grav = 1506.5661, roche = true, roche_a = 8.322578e11 cm**;
  p_ref surface pole 1.768, terminator-eq 1.772, substellar 1.785, antistellar 1.784 R_J; inner wall
  pole 205 bar / 4275 K, terminator 250 bar, substellar 440 bar / 4778 K; L1 5.180e10 (7.25 R_J),
  x1max/r_L1 0.367 (inside the lobe). Kepler M_* 2.0517 M_sun (+1.6 %).
- **Anchor inputs (Caltech re-anchors round 1 with the tide ON), per metallicity deck:**

| input | 3x nominal | 1x arm | source |
|---|---|---|---|
| M_p | 1.66 M_J (G M_p 2.102996e23) | same | Hellier+2019 (RV) |
| R_p (terminator, at p_ref) | 1.77 R_J | same | Lothringer+2025 |
| p_ref (TESS band) | 6.651e-04 bar (rce/r3x_B) | 1.028e-03 bar (arms/r1x_C) | transit.txt |
| a | 0.055633 AU (= a/R_* 7.037 x R_* 1.70) | same | Lothringer+2025 |
| P_orb = P_rot | 3.3448371 d | same | Damasceno+2024 |
| Teff | 9200 K | same | Lothringer+2025 |
| R_* | 1.70 R_sun | same | Lothringer+2025 |
| M_* | 2.02 M_sun (Kepler check only) | same | Lothringer+2025 |
| profile / EOS | ic/ic_wasp178b_3x_f5.txt / met3 dump | ic/ic_wasp178b_1x_f5_arm.txt / deepconv_0925 dump | |
| **output ap = x1min / grav / roche_a** | **1.181477e10 / 1506.5661 / 8.322578e11** | **1.181312e10 / 1506.9869 / 8.322578e11** | roche_anchor.txt |
| Kepler check (implied vs literature M_*) | 2.0517 vs 2.02 (**+1.57 %**) | same | |

  Keys to switch the tide on: `problem/roche=true problem/roche_a=8.322578e11 problem/ap=1.181477e10
  problem/grav=1506.5661 mesh/x1min=1.181477e10` (decks stay on the 182477553 defaults, roche off).

## 7. Dynamo estimate (for later MHD)

B ~ F_int^(1/3) relative to Jupiter, T_int = 553.4 K:
- vs F_J 7.485 W/m^2 (Li+2018): F ratio 711, **B/B_J = 8.9 -> B_dip,eq ~ 37 G** (rms ~62 G);
- vs 5.44 W/m^2 (Hanel+1981): 978, 9.9 -> 42 G / 69 G.
- Seidel+2026 Fig. 3: B_max 1.6-2.7 G for WASP-178b: the scaling is ~15-20x above it.

## 8. Cost and rotations (`grid/nominal/cost_120.txt`, `gpu/`)

- Estimate (`scripts/cost.py`): dt = 13.79 x (8.41/11.43) = 10.15 s; 28470 cycles/rot; nx1 120 ->
  3.82x W121 per rotation, 14.1 min/rot.
- **Measured: apudev job 12528431 (vipa1326, same job and node as the WASP-76b smoke, binary
  athena_dhj_gpu72_18247755_hfresh, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1, the 3x nominal deck
  verbatim, 5.5 min wall; `gpu/ana_out.txt`): 10900 cycles, 0.328 rot, 0.0295 s/cycle, mean dt 8.67 s
  (last 8.67 s), 33331 cycles/rot, 16.4 min/rot on 1 node.** Against W121 hfresh (0.0277 s/cycle,
  13.79 s; other job/node): **4.44x W121 per rotation** (per cycle 1.06x, dt 0.63x, rotation 2.62x
  longer), **1.69x per unit simulated time.** The measured dt is 15 % below the radial-CFL estimate.
- Rotations: t_rad(100 bar) = 1.87e6 s = 6.5 rot (W121 2.42e6 s = 22 rot); 300 W121 rot = 13.6
  t_rad(100 bar) = 88 WASP-178b rotations. **Recommended 120 rot (18.5 t_rad(100 bar), 1.05x the
  W121 300-rot physical time); tlim in the deck = 120 rot.** Cost ~33 h on 1 viper node
  (W121 300 rot 18.5 h): 1.8x the W121 campaign. Check the 10-100 bar drift at 60 and 90 rot.
- Events (smoke, 4 log rows): dfloor 3.2e8, efloor 6.2e8, tfloor 2.8e8 (the top starts on the floor as
  W121's, same top pressure); vceil, eos_fail, c2p_it, fofc 0.
- 1x arm: same grid class; cost within a few % of the 3x deck (not measured on GPU).

## 9. Gates

CPU binary `/viper/ptmp2/jinma/dhj_planets_1101/kelt20b_work/athena_dhj_cpu_18247755` (fin51-r-fix
182477553, KELT-20b worker's build, md5 fac6d83e..., used read-only); `gate/gate.sh <deck> <name>`.

| gate | 3x nominal | 1x arm |
|---|---|---|
| parse (`gate/parse_*.log`) | rc 0 | rc 0 |
| CPU fresh start, C16 gate mesh (same radial grid and IC fill), 8 ranks, 2 cycles, ghost zones | rc 0, no FATAL, no events, dt 11.09 s | rc 0, no FATAL, no events, dt 11.18 s |
| T vs IC table (`gate/gate_ic_*.txt`): max |dT/T| 1e-6..1e-3 / 1e-3..1 / 1..300 bar | 1.3e-2 / 6.9e-4 / 6.8e-4 | 1.0e-2 / 2.1e-3 / 8.1e-4 |
| deepest active cells | 225.0 bar 4334 K (pole), 273.6 bar 4458 K, on the table | 225.4 bar 4333 K, 272.8 bar 4455 K |
| r(p_ref) vs tool pole / equator | +1.1e-5 / +1.2e-5 | +1.1e-5 / -6.8e-6 |
| hydrostatic residual (median, centred difference) | 0.45-0.89 % | 0.43-0.95 % |
| min density / top cells | 1.37e-18 / 2.3-3.1e-18 (top on the floor, as W121) | 1.53e-18 / 2.6-3.4e-18 |
| coverage (`gate/coverage_*.txt`) | IC 2140-6158 K, day (x1.3) to 5448 K, night (x0.6) 1325 K, hot bound 4508 K: inside ck/CE T 100-10100 K; p > 1e-8 bar inside (49 top rows clamped, as W121) | IC 2091-6154 K, day 5468 K, night 1369 K: inside |
| GPU smoke (apudev job 12528431) | rc 0, no FATAL, 10900 cycles, eos_fail/c2p_it/fofc/vceil 0 | not run (CPU-gated only) |

## 10. Flags / what needs a decision

1. **Seidel's Teq 2398 K is not reproducible** from any published parameter set; code Teq 2452.3 K
   (Lothringer star). Host Teff is the largest uncertainty (KELT SED 8640 K -> Teq 2303 K).
2. **Hardest UV host:** 16.5 % of the stellar flux below 0.26 um is deposited by band 10 at
   1e-4..5e-2 bar (50 % at 3.5e-3 bar); in reality much of it is absorbed above 1e-5 bar. A UV pseudo-band
   test (as kelt20b_work/rce_f5_uv.py) is advisable before trusting the 1e-6..1e-3 bar temperatures.
3. Metallicity weakly constrained and element-dependent (O-rich, C-poor); 3x nominal / 1x arm (user).
4. Bond albedo 0.15 (CHEOPS A_g 0.1-0.35, A_g < 0.23).
5. Cost: 4.4x W121 per rotation (measured), 120 rot ~33 h on 1 viper node.
6. The RCEs ran at keys 0.2-0.4 % from the final ones (no extra pass for the nominal). The first face-0
   launch was killed by another worker's pkill (05:22) and restarted from its checkpoint; no killed
   state was used as an IC.
7. Dynamo scaling ~37 G vs Seidel's 1.6-2.7 G.
8. The tide is negligible (<= 4.5 % of g at the top): the Roche port is optional here.

## 11. Files and commands

- `deck/wasp178b_hfresh.athinput` (3x NOMINAL), `deck/wasp178b_hfresh_1x_arm.athinput` (1x arm),
  `deck/ovr_3x.json`, `deck/ovr_1x_arm.json`.
- `ic/ic_wasp178b_3x_f5.txt`, `ic/ic_wasp178b_1x_f5_arm.txt`.
- `ck/CK3X/`, `ck/CK1X/`, `ck/CK10X/` (unused): symlinks + `sw_flux/sw_band_flux_W178_11.txt`.
- `sed/` (make_sed.py, ATLAS9 files, sw file, make_sed.out; `old_damasceno/` superseded).
- `rce/r3x_A`, `rce/r3x_B` (+ `f5/`), `rce/uvdep.txt`; `arms/r1x`, `arms/r1x_C`, `arms/r10x` (stopped).
- `grid/nominal/` (nx1 120 keys, faces, cost_120.txt, system.txt, roche_anchor.txt), `grid/arm1x/`.
- `gate/` (gate.sh, gate_ic.py, parse/run/IC/coverage outputs for nominal and arm1x; `superseded/`),
  `gpu/` (ana_out.txt; job script and analysis in ../wasp76b/gpu/, run dir `gpu/run_wasp178b/`).
- Scripts as in ../wasp76b (chain.sh, finish_f5.sh, run_arm.sh, arms/stop_at40.sh, post.sh + planet.env,
  scripts/).
- Command lines actually used (in order):
  ```
  chain.sh rce/r3x_A 1510 1.18e10 rce/ck_r3x_A1.npz 40 nof5                 (restart of the killed launch)
  chain.sh rce/r3x_B 1507.195 1.181230e10 rce/r3x_A/ckrce_nq2_B.npz 400       (stopped at it 146; then 2 it)
  finish_f5.sh rce/r3x_B 1507.195 1.181230e10                                (face-5 to it 40, stopped)
  FIN_RESTOL=1 finish_f5.sh rce/r3x_B 1507.195 1.181230e10 rce/r3x_B/f5/ckpt_it40.npz
  post.sh rce/r3x_B 3x ic_wasp178b_3x_f5.txt grid/nominal "112 116 120 124"   (nx1 120 chosen; cost_120 rerun)
  run_arm.sh arms/r1x 1507.195 1.181230e10 rce/r3x_B_ckpt400.npz scripts/rce_env_1x.sh "<label>"
  FIN_RESTOL=2e-3 ENVF=scripts/rce_env_1x.sh finish_f5.sh arms/r1x 1507.195 1.181230e10 arms/r1x/f5/ckpt_it40.npz
  FIN_RESTOL=1 ENVF=scripts/rce_env_1x.sh finish_f5.sh arms/r1x_C 1512.8786 1.179009e10 arms/r1x_C/f5/ckpt_run1.npz
  post.sh arms/r1x_C 1x ic_wasp178b_1x_f5_arm.txt grid/arm1x "116 120 124"
  scripts/make_ovr.py ... deck/ovr_3x.json 3x ; scripts/make_ovr.py ... deck/ovr_1x_arm.json 1x ; scripts/make_deck.py ...
  gate/gate.sh deck/wasp178b_hfresh.athinput s3x   (renamed nominal) ; gate/gate.sh deck/wasp178b_hfresh_1x_arm.athinput arm1x
  sbatch ../wasp76b/gpu/job.sh                                                 (job 12528431)
  ```
- `MD5SUMS`: md5 of every deck, IC, sw-flux file and script here. Files > 50 MB: only the gate and smoke restarts (`gate/run_*/rst/`, `gpu/run_*/rst/`), not deliverables.
