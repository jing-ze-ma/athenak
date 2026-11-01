# toi1518b: hydro-first deep-hot-Jupiter setup for TOI-1518b (dhj_planets_1101, viper 2026-11-01)

Built with the WASP-121b / KELT-20b pipeline (RCE -> transit p_ref -> planet_setup.py -> IC -> grid -> deck -> gates)
on the WASP-121b hydro fresh deck `w121hfresh_1101/deck/hfresh.athinput` (fork/fin51-r-fix 182477553).
Caltech runs it; viper did CPU gates and one apudev smoke only. Nothing was pushed or committed.

## 0. What to run
- **Deck:** `deck/toi1518b_hfresh.athinput` (3x solar, the adopted metallicity; section 3).
  - Every changed key is marked `# toi1518b:` with the WASP-121b value.
  - Everything else is the hfresh deck verbatim: hydro, lhllc, `<hydro>/am_vector = true`,
    `problem/am_transport = false` (required), cs_seam keys, top sponge 1e-7..1e-6 bar (pressure based), no bottom sponge,
    ck T4/C2/E4 (`ck_impl_every 4`, `xstep 8`, `once`, `ck_nquad 2`), floors dfloor 1e-18 / pfloor 1e-7 / tfloor 200 K,
    eos_logd_min -19, 3x EOS composition (eos_xh 0.7188, eos_yhe 0.2420), met / rad_met 0.4771.
- **1x sensitivity deck:** `deck/toi1518b_hfresh_1x.athinput` (section 3d; CPU-gated, not GPU-smoked; needs `ck/CK1X`).
  Caltech has the 1x tables in its repo `data/exo_fms_ck`; add `sw_band_flux_T1518_11.txt` to that set's `sw_flux/`.
- **Run directory:** must contain `ck` and `ic` (copies of, or symlinks to, `ck/` and `ic/` here). Table paths are relative,
  as in the s10fresh package; do not use `-d`.
  - `ck/CK3X/` here = symlinks to the WASP-121b package tables (`w121_s10fresh_1027/pkg/ck/CK3X/{ck,CE_tables,cia,ray}`),
    plus a real `sw_flux/` with **the new file `sw_band_flux_T1518_11.txt`**. Caltech already has CK3X from the s10fresh
    package: copy `ck/CK3X/sw_flux/sw_band_flux_T1518_11.txt` into its `sw_flux/` and use it.
  - `ic/ic_toi1518b_3x_f5.txt` (7 kB).
- **Commands (as the s10fresh README):**
  ```bash
  cd <rundir> || exit 1      # holds ck -> <pkg>/ck, ic -> <pkg>/ic
  srun -n 4 <athena> -i <pkg>/deck/toi1518b_hfresh.athinput time/nlim=200 > smoke.log 2>&1     # smoke
  srun -n 4 <athena> -i <pkg>/deck/toi1518b_hfresh.athinput -t <hh:mm:ss> > run.log 2>&1         # production
  ```
  Outputs: hst every 0.01 rot (1643.86 s), log every 0.1 rot, hydro_w bins every 0.5 rot, rst every rot (164385.6 s).
  `tlim` = 200 rot (section 9).
- **Binary:** any fin51-r-fix build at or after 182477553 with PROBLEM=deep_hot_jupiter_rt. No code change is needed.

## 1. System parameters

| quantity | adopted | source | other values |
|---|---|---|---|
| P | 1.90261131 +- 0.00000043 d | Simonnin+2025 A&A 698, A314, Table 7 (EXOFASTv2: TESS + FIES + SOPHIE) | 1.902603 (Cabot+2021 AJ 162, 218) |
| a/R_* | 4.109 +0.021/-0.020 | Simonnin+2025 Table 7 | 4.291 +0.057/-0.061 (Cabot+2021) |
| a | 0.03712 +- 0.00082 AU (a/R_* x R_* = 0.03711 AU; Kepler with M_* + M_p: 0.03710 AU) | Simonnin+2025 | 0.0389 +- 0.0011 (Cabot) |
| R_p/R_* (TESS) | 0.09939 +0.00050/-0.00051 | Simonnin+2025 | 0.0988 (Cabot) |
| R_* | 1.942 +- 0.043 R_sun | Simonnin+2025 | 1.950 +- 0.048 (Cabot) |
| M_* | 1.88 +0.13/-0.12 M_sun | Simonnin+2025 | 1.79 +- 0.26 (Cabot) |
| T_eff | 7299 +98/-100 K | Simonnin+2025 | 7300 +- 100 (Cabot) |
| log g_* | 4.136 +- 0.011 | Simonnin+2025 | 4.1 +- 0.2 (Cabot) |
| [Fe/H]_* | -0.10 +- 0.12 | Simonnin+2025 (= Cabot+2021) | |
| v sin i_* | 85.1 +- 6.3 km/s (F0/A host, fast rotator; misaligned orbit, b 0.88-0.90, precessing, Watanabe+2024) | Cabot+2021 | |
| **M_p** | **1.83 +- 0.47 M_J** (first RV mass, SOPHIE) | Simonnin+2025 | < 2.3 M_J (2 sigma, Cabot) |
| **R_p** | **R_p/R_* x R_* = 1.8783 R_J = 1.342834e10 cm** (TESS band) | Simonnin+2025 (1.878 +- 0.042) | 1.875 +- 0.053 (Cabot) |
| e | circular; e < 0.01 (2 sigma) | Cabot+2021 | |
| Teq (zero albedo, full redistribution) | **2546.1 K** = 7299 / sqrt(2 x 4.109) | derived; Simonnin Table 7 gives 2546 +35/-36 | **check value 2491 +- 38 K** (Seidel+2026; = Cabot's 7300 K with a/R_* 4.291: 2491.9 K) |
| omega | 3.822223e-5 s^-1 (P = 164385.6 s), tidally locked | derived | |
| g(R_p) | 1285.7 cm s^-2 (log g 3.109) | derived | Simonnin 3.110 +0.097/-0.13 |
| T_int | 522.0 K (Thorngren+2019 fit in the code, get_Tint(Teq)) | code | |

- **Teq vs the Seidel 2026 value:** our 2546.1 K is 2.2 % (55 K) above the paper's 2491 K. The whole difference is a/R_*:
  Seidel uses Cabot's 4.291, Simonnin's newer joint fit gives 4.109 (same T_eff). We adopt the newer fit because P, a/R_*,
  R_p/R_*, M_p and the star come from one consistent solution. A run at 2491 K would need its own RCE (`scripts/rce_env.sh`,
  RCE_TEQ) and would move T by ~2 %.
- The planet has 14 species detected at high resolution (Simonnin+2025, MAROON-X), Fe I / Fe II (PEPSI, Basinger+2025
  MNRAS 543, 4136), Fe I (Cabot+2021). Simonnin+2025 infer strong drag (tau_drag ~ 1e3-1e4 s) from the Fe / Fe+ trails;
  Ca+ reaches 1.97 +- 0.04 R_p, i.e. the Roche lobe (section 8). Seidel+2026 Fe I wind: 3.82 +- 0.50 km/s.
- Literature PDFs/text: `lit/` (cabot2021.pdf/.txt, simonnin2025.txt, basinger2025.txt).
- `scripts/system.py` -> `system.txt`: all derived numbers above and the tide table.

## 2. Planet mass and the gravity range

- **Nominal M_p = 1.83 M_J** (Simonnin+2025). The error is +-0.47 M_J = **+-26 %**.
- Same IC and p_ref, planet_setup at the 1 sigma ends (`sens/setup_mp*.txt`), grid112 at nx1 120 and cost.py
  (`sens/mp*/grid.log`, `sens/mp*/cost.txt`):

| M_p [M_J] | grav (at ap) | g(R_p) | ap = x1min [cm] | x1max [cm] (W121 rule) | x1max / R_p | radial dt est. -> dt | cost / rot vs W121 | min/rot (1 viper node) |
|---|---|---|---|---|---|---|---|---|
| 1.354 (-26 %) | 1161.83 | 951 | 1.215076e10 | 2.353278e10 | 1.75 | 12.13 s -> 14.6 s | 1.51 | 5.6 |
| **1.83** | **1493.49** | **1285.7** | **1.245919e10** | **1.968103e10** | **1.47** | **9.07 s -> 10.9 s** | **2.02** | **7.5** |
| 2.306 (+26 %) | 1826.29 | 1620 | 1.264766e10 | 1.795536e10 | 1.34 | 7.33 s -> 8.9 s | 2.49 | 9.2 |

- H ~ 1/g: the domain height (250 bar to the W121 top pressure) runs from 1.14e10 to 0.53e10 cm, dt about as 1/g, so the
  cost per rotation spans 1.5-2.5x W121 (+-24 %). nx1 stays 120 (it counts cells per H). The IC in pressure also depends
  on g (tau = kappa p/g); these rows reuse the nominal IC, a consistent off-nominal run needs its own RCE.
- Recommendation: run the nominal mass. The mass is now measured (unlike KELT-20b), so no bracket run is needed for
  hydro first.

## 3. Metallicity

### 3a. Evidence (ranked as asked: planet first, then mass-metallicity, then host)
1. **The planet's own atmospheric metallicity:**
   - Simonnin+2025 (MAROON-X, high-res retrieval): log10 Fe = -4.88 +0.63/-0.76, i.e. **0.07-1.62x solar in Fe (1 sigma)**,
     median ~0.3x; "consistent with solar (and stellar) metallicity". The other species are biased (different Kp/Vsys
     from Fe; Fe+ and Ca+ far above solar, Mg low) and are not used by the authors.
   - Basinger+2025 (PEPSI, PIRANGA III): detections and winds only, solar H/He assumed, **no metallicity**.
   - Cabot+2021: Fe detection only. **No JWST or HST spectrum and no volatile (H2O / CO) abundance exists** (searched
     2026-11-01). So the species that carry most of the model's IR opacity (H2O, CO, H-) are unconstrained.
2. **Mass-metallicity expectation (no volatile retrieval exists):**
   - Thorngren+2016 bulk metallicity Z_p/Z_* = 9.7 M^-0.45: **~7x stellar in bulk** at 1.83 M_J (scatter about a factor 2).
     The envelope is at most the bulk value.
   - Welbanks+2019 (H2O, Na, K): giant planets of ~1-2 M_J are substellar to stellar in H2O; the trend is weak and a 2024
     reanalysis (arXiv:2402.08292) cannot distinguish it from no trend.
   - Together: the atmosphere is expected between ~0.3x and a few x solar.
3. **Host:** [Fe/H]_* = -0.10 +- 0.12 (0.79x solar). Context only.

### 3b. Table sets that exist on viper (all complete: ck premixed g8 11-band hiT2 + FastChem CE hiT2 + EOS dump)

| set | ck / CE tables | EOS composition (eos_xh, eos_yhe; met) | EOS dump for the RCE / setup |
|---|---|---|---|
| 1x | repo `data/exo_fms_ck/{ck,CE_tables}` (Premixed_1x_g8_11_hiT2, FastChem_ck_1x_int_hiT2); also `csfloor_1008/data/ck/CK1X`, `kelt20_0925/ckdata1` | 0.7381, 0.2485; met 0.0 | `deepconv_0925/dump` |
| **3x** | `ck3x_0930/ckdata3` = `w121_s10fresh_1027/pkg/ck/CK3X` (Premixed_3x_g8_11_hiT2, FastChem_ck_3x_int_hiT2) | 0.7188, 0.2420; met 0.4771 | `w121prod_0929/met3/dump` |
| 10x | `wasp121_0925/ckdata10` (Premixed_10x_g8_11_hiT2, FastChem_ck_10x_int_hiT2; gen10.sh) | 0.6588, 0.2218; met 1.0 | `wasp121_0925/met10/dump` |

- Missing: anything between these, and 30x / 100x (premixed upstream tables exist only at 1x/10x/100x/1000x; 3x was
  generated). Nothing is missing for a 1x / 3x / 10x bracket.
- The ck tables are clamped at P 1e-8..1e3 bar and T 100..10100 K (hiT2): TOI-1518b's column (2300-4600 K below 1e-6 bar,
  4100 K at the top) is inside.

### 3c. Choice
- **Adopted (deck): 3x solar.** It is 0.58 dex above solar and **0.68 dex (3.8x) above the host** ([Fe/H]_* -0.10).
  - The planet's only measurement is one refractory (Fe) at the limb, in the upper atmosphere, by a retrieval whose other
    species are biased; it puts 3x about 1.5 sigma above its median and 0.27 dex above its 1 sigma upper end.
    Volatiles are unconstrained, and the mass-metallicity relation allows up to ~7x in bulk.
  - 3x is the WASP-121b production metallicity, so the two planets differ only in planet, star, albedo and orbit, which is
    what the Seidel 2026 wind-vs-Teq comparison needs.
- **Recommended sensitivity pair: 1x and 3x** (1x = +0.10 dex above the host, at the top of the Fe 1 sigma range).
  10x is disfavoured by the Fe retrieval (> 2 sigma) and not recommended.
- **User decision:** if the Fe retrieval is to be followed strictly (planet evidence first), 1x is the nominal and 3x the
  sensitivity arm instead. Both decks are in `deck/`.

### 3d. 1x arm
- **Deck `deck/toi1518b_hfresh_1x.athinput` (ready; CPU-gated, not GPU-smoked).** Same as the 3x deck except the planet
  keys from the 1x chain, the 1x grid and the metallicity keys of the WASP-121b 1x input: eos_xh 0.7381, eos_yhe 0.2485,
  met 0.0, rad_met 0.0, `ck_data_dir = ck/CK1X`, `ck_table = ck/CK1X/ck/Premixed_1x_g8_11_hiT2.txt`,
  `ck_ce_table = CE_tables/FastChem_ck_1x_int_hiT2.txt`, `ic_profile = ic/ic_toi1518b_1x_f5.txt` (`deck/ovr_1x.json`).
  `ck/CK1X/` = symlinks to the repo `data/exo_fms_ck/{ck,CE_tables,cia,ray}` + `sw_flux/` with the T1518 file.
- Chain (`chain_1x.sh`, `finish_f5_1x.sh`, env `scripts/rce_env_1x.sh`, EOS dump `deepconv_0925/dump`):
  rce1xA face 0 from the 3x column (stopped at it 100, res 1.8e-4: limit cycle) -> rce1xB face 5 at the 3x keys
  (limit cycle at res ~1.4e-3, accepted with RCE_RESTOL 2e-3 = the level the W121 3x f5 run stopped at) -> planet_setup
  -> **rce1xC** face 5 at ap 1.248208e10 / g 1488.0156 (keys then move by 1.2e-5: final).
  Energy closure 0.967 (3x: 0.9997): **the 1x column is less converged** (limit cycle at i 104, 1.8e-3 bar).
- 1x results: **p_ref (TESS) = 1.824e-3 bar** (3x: 1.072e-3; limb -/+300 K: 1.49e-3 / 2.27e-3); ap = x1min = 1.248193e10,
  grav 1488.0510; x1max (W121 rule) 2.007760e10; nx1 120 (max dr/target 1.016); radial dt estimate 9.08 s, cost estimate
  2.01x W121 per rotation (the same as 3x within 1 %; the measured 3x factor 2.41 applies to both).
- Global-mean T [K] at 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 1, 10, 100, 250 bar: 4089, 4071, 3550, 3317, 2404, 2288, 2245,
  2871, 4036, 4587 (3x: 4091, 4051, 3525, 3058, 2309, 2293, 2319, 2884, 4038, 4588). 1x is hotter at 1e-3..1e-2 bar
  (inversion deeper, as p_ref) and 70 K cooler at 1 bar.
- **Caveat: the 1x IC has a 420 K step at 1.45e-3 bar** (3221 -> 2798 K between adjacent rows; the face-0 1x column has the
  same step, so it is the 1x inversion onset, not a face-5 artefact). The WASP-121b 1x IC `ic_w121_1x_f5.txt` has
  steps of the same kind (335 K at 2.7e-4 bar, 277 K at 6.5e-4 bar). t_rad there is ~1e2 s, so the run relaxes it within
  tens of cycles. The IC gate therefore shows max |dT/T| 4.9 % in 1e-3..1 bar (the step sampled by the cells).
- Gates (`gate/parse_1x.log`, `gate/run_fresh_1x`, `gate/gate_ic_1x.txt`): parse rc 0; CPU fresh start 2 cycles rc 0,
  no FATAL, no events, dt 11.80 s; deepest active cells 223.9 bar / 4516 K and 408.0 bar / 4919 K on the table;
  r(p_ref) over 9600 columns +4e-6 / -8e-6 against the tool's pole / equator; min density 1.16e-18.

## 4. Bond albedo
- No bolometric albedo exists. TESS (Cabot+2021): eclipse 364 +- 28 ppm, nightside 43 +- 27 ppm, T_day 3237 +- 59 K
  (A_g = 0), T_night 1700 +700/-1200 K. With T_0 = T_eff (R_*/a)^1/2 = 3601 K and Cowan & Agol (2011) these give
  (1 - A_B) > 1, i.e. A_B ~ 0 (formally negative): a TESS brightness temperature of an inverted dayside is above the
  bolometric one, so this is a lower bound.
- The only bolometric UHJ albedo, WASP-121b's 0.277 +- 0.016 (JWST, Splinter+2025), and KELT-20b's adopted 0.32 are higher;
  the code's Parmentier fit gives ~0.01 at this Teq.
- **Adopted A_B = 0.10, range 0-0.3.** `problem/albedo = 0.10`, the RCE uses the same absorbed flux (RCE_ALB). Teq stays the
  zero-albedo 2546.1 K.

## 5. Host SED and UV coverage (`sed/`)
- `sed/make_sed.py` = kelt20_0925/sed/make_sed.py generalised to any T_eff and [M/H] in [-0.5, 0] (ATLAS9 ck04 `ckp00`
  and `ckm05` from STScI CDBS, trilinear in T_eff, log g, [M/H]); regression: it reproduces sw_band_flux_K20_11.txt
  exactly (max rel 0.0).
- `sw_band_flux_T1518_11.txt`: T_eff 7299, log g 4.136, [M/H] -0.10; same band edges and order as sw_band_flux_W121_11.txt;
  integral 0.9993 sigma T_eff^4 (`sed/make_sed.out`). `ck_star_teff = 0`, `ck_swflux = sw_band_flux_T1518_11.txt`.
- **4.31 % of the stellar flux lies below 0.26 um** (blackbody 5.2 %); it is folded into band 10 (0.26-0.42 um), which
  carries 24.5 % in total (blackbody 28.8 %; Balmer jump). KELT-20 had 14.4 % below 0.26 um, WASP-121 ~3 %.
  So the misplaced (too deep) UV heating is ~1/3 of KELT-20's; the same caveat applies at a smaller level.
  Flux beyond 324.68 um is folded into the last band (negligible).

## 6. RCE, transit radius and IC (`rce/`)
- Chain (`chain.sh`, `finish_f5.sh`): rce_w121.py route B, ck_nquad 2, 20 cells/dex (env `scripts/rce_env.sh`: Teq 2546.1,
  A_B 0.10, 3x tables, the T1518 sw file) -> rce_f5.py (face-5 + vacuum top exact RCE, F5_INLINE=1, as the WASP-121b
  `ic_w121_3x_f5.txt`) -> transit.py (TESS 0.6-1.0 um, the R_p band) -> make_ic.py -> planet_setup.py.
  - rceA (guess ap 1.22e10, g 1558; stopped at it 80, res 9e-6) -> rceA2 (3-it finisher) -> planet_setup ->
  - rceB at ap 1.249739e10 / g 1484.37: face 0 converged in 46 it (res 8.5e-8); face 5 stopped at res ~1e-5 (plateau in
    the top cells; restarted from its it-40 checkpoint with RCE_RESTOL 2e-5; the W121 3x f5 run stopped at 2e-3) ->
  - **rceC** (adopted) at ap 1.246085e10 / g 1493.091: face 5, energy closure OLR_exact / (absorbed + L_int) = 0.99971
    (W121 3x: 0.98777). The next planet_setup moves ap by -1.3e-4 and grav by +2.7e-4 (1.6e5 cm, 0.003 H): final.
  - Python: python-waterboa/2025.06 (numpy 2; the RCE scripts need np.trapezoid). The face-0 RCE (rceA, rceA2) ran under
    system python 3.9 (numpy 1.26).
- **p_ref (TESS) = 1.072e-3 bar** (R_eff definition, as W121/KELT); tau_chord = 0.56 level 1.04e-3 bar.
  Limb T -300 / +300 K above 0.1 bar: 7.8e-4 / 1.30e-3 bar (0.2 H). H(p_ref) = 1.47e8 cm.
  Other bands: NRS1 1.42e-3, WFC3 G141 1.76e-3, G395H 1.36e-3 bar (`rce/rceC/transit.txt`).
- p_ref isobar with rot_potential: pole 1.337697e10, equator 1.348021e10 cm, terminator area-equivalent = R_p.
- Global-mean T [K] (`rce/rceC/f5/ckrce_nq2_B.txt`):

| profile | 1e-8 | 1e-6 | 1e-4 | 1e-3 | 1e-2 | 0.1 | 1 | 10 | 100 | 250 bar |
|---|---|---|---|---|---|---|---|---|---|---|
| TOI-1518b 3x, A_B 0.10 | 4091 | 4051 | 3525 | 3058 | 2309 | 2293 | 2319 | 2884 | 4038 | 4588 |
| WASP-121b 3x, A_B 0.277 (ic_w121_3x_f5) | 3712 | 3459 | 3170 | 2278 | 2091 | 2071 | 2045 | 2628 | 3832 | 4375 |

  Radiative-convective boundary ~5.4 bar in both (grid face).
- **IC** `ic/ic_toi1518b_3x_f5.txt`: the rceC face-5 column + isothermal rows to 9.5e-11 bar + the exact 3x-EOS adiabat
  below the RCE bottom (283 bar) to 2241 bar (269 uniform rows, 0.0499 dex).
- **Keys (`rce/rceC/setup.txt`):** ap = x1min = 1.245919e10 cm, grav = 1493.4888, Teq 2546.1, omega 3.822223e-5.
  planet_setup's own x1max rule (3 % inside r(1e-9 bar)) would give 1.658411e10; the deck uses the W121 raised-top rule
  (section 7).

## 7. Radial grid (`scripts/grid112.py`, `grid/final/`)
- Rule: carry the W121 hfresh grid (C32 x nx1 112, x1max 2.1e10, 8 coefficients + 2 bumps) over at **equal pressure**:
  1. cells per H of every W121 cell on its own fresh-start column (ic_w121_3x_f5, point mass);
  2. W121 x1max = 2.1e10 sits at p = 6.1e-18 bar on that column (isothermal continuation above the IC top);
  3. TOI-1518b x1max = r(6.1e-18 bar) on its IC column = **1.968103e10** (the same number of top scale heights above the
     1e-6 bar level, which the W121 top was sized for, 3 H_p above the highest dayside 1e-6 bar level in 3-D);
  4. minimax of max(dr / dr_target) over the 8 coefficients + 2 bumps, dr_target = H / (W121 cells per H at that p).
  - Regression: on WASP-121b itself it returns x1max 2.100000e10 and the hfresh keys to 1e-3 (max ratio 1.004 from the
    6-decimal rounding; `gate/grid112_regress_w121.txt`).
- The exact-target cell count is 112.1 for both planets. nx1 112 reaches max ratio 1.078 (some cells 8 % coarser than
  W121); **nx1 = 120 reaches 1.0006: adopted** (+7 % per cycle). MeshBlocks 120 x 16 x 16, 24 MBs.
- Cells per H on the IC column, nx1 120 (W121 hfresh in brackets): 1e-6..1e-4 bar min 3.10 (2.98), 1e-4..1e-2 4.26 (4.04),
  1e-2..1 4.23 (4.10), 1..10 3.92 (3.80), 10..300 4.11 (4.01); 1e-9..1e-6 1.81 (1.67); 32 cells above 1e-6 bar (31).
  dr 1.25e7 .. 1.73e8 cm.
- Faces: `grid/final/faces_120.txt`; keys: `grid/final/gridkeys_120.txt`.

## 8. Tide (Roche) and domain top (`system.txt`)
- a/R_* = 4.109; q = M_p/M_* = 9.29e-4. Roche sizes: L1 (Hill) a (q/3)^1/3 = 3.756e10 cm = 2.80 R_p; Eggleton volume
  radius 2.621e10 cm = 1.95 R_p (Simonnin+2025 quote R_Roche ~ 1.96 R_p; Ca+ observed at 1.97 R_p).
- Tide / gravity (point-mass g; 2(r/a)^3/q = stellar tide at the substellar point; 3(r/a)^3/q adds the orbital centrifugal
  term; omega^2 r/g = the spin centrifugal term already in rot_potential):

| r | r/R_p | r/r_L1 | 2(r/a)^3/q | 3(r/a)^3/q | omega^2 r/g |
|---|---|---|---|---|---|
| x1min 1.2459e10 | 0.93 | 0.33 | 0.024 | 0.037 | 0.012 |
| R_p 1.3428e10 | 1.00 | 0.36 | 0.030 | 0.046 | 0.015 |
| r(1e-6 bar) 1.4771e10 | 1.10 | 0.39 | 0.041 | 0.061 | 0.020 |
| **x1max 1.9681e10** | **1.47** | **0.52** | **0.096** | **0.144** | 0.048 |

- The tide is 3-5 % of g at the photosphere and 10-14 % at the domain top: **it matters at the top**, but less than for
  WASP-121b (same formula: 8 % at R_p, 40 % (59 % with the orbital term) at its x1max 2.1e10 = 0.84 r_L1, a/R_* 3.8).
  The Roche port at Caltech is more urgent for WASP-121b than for TOI-1518b. At M_p 1.354 M_J the top is 1.75 R_p, 0.69 r_L1 and the tide ~22 % of
  g there. The domain stays well inside the Roche lobe (0.75 of the Eggleton radius at the nominal mass).

### 8b. Roche anchor, 1x (`gate/roche_anchor_1x.txt`; Caltech NOTE 272395a5b: round-1 decks run with the tide on)
- Tool: `inputs/production/tools/roche_anchor.py` from fork/fin51-r-fix c8aac483c (git show). Inputs are the section 1
  values, the 1x IC `ic/ic_toi1518b_1x_f5.txt` (md5 9ec9fa67...), the 1x EOS dump `deepconv_0925/dump/eos_table.txt`
  (md5 697ca675...), the 1x p_ref (TESS) 1.824e-3 bar, and x1max 2.007760e10 for the L1 check:
  `roche_anchor.py --mp 1.83 --rp 1.8783 --pref 1.824e-3 --a_au 0.0371092 --porb 1.90261131 --teff 7299 --rstar 1.942
   --mstar 1.88 --profile ic_toi1518b_1x_f5.txt --eos eos_table.txt --x1max 2.007760e10`
  (a_au = ARS R_* = 5.551462e11 cm / AU, rounded to 7 digits, so the tool's a is 5.551457e11).
- **Result (tide on): ap = x1min = 1.257097e10, grav = 1467.0475, roche_a = 5.551457e11.** The tide-off anchor of 3d
  (1.248193e10 / 1488.0510) is reproduced to all digits by the same tool with --no_roche, which checks the tool against
  planet_setup.py. With the tide, ap moves +0.71 % and grav -1.41 %.
- The p_ref surface (R_J) is pole 1.87122, terminator equator 1.88545, substellar 1.93325, antistellar 1.93157. Inner-wall
  pressure is 128 bar at the pole, 250 bar at the terminator equator, 1287 bar substellar and 1234 bar antistellar. The
  deepest is 5855 K, on the adiabat extension below the IC table. L1 = 3.668e10 cm, ap/r_L1 0.343, x1max/r_L1 0.547.
  Phi_tot on x1max is below Phi(L1), so the whole top is inside the lobe.
- **Kepler M_*:** the tool (like the pgen) derives G M_* = omega^2 a^3 - G M_p. That implies M_* = 1.8816 M_sun vs the
  literature 1.88 (+0.09 %). The mismatch is negligible here; WASP-121b has +2.63 % and WASP-189b 14 %.
- Grid at the Roche anchor (`grid/roche_1x/`, grid112.py, nx1 120, PLANET_EOS 1x): x1max 2.030897e10, max dr/target
  1.0098, min cells/H 3.29 at p >= 1e-6 bar, radial dt estimate 9.27 s. Keys in `grid/roche_1x/gridkeys_120.txt`.
  It is **not** run or gated, and the deck is unchanged (tide off).

## 9. Cost and rotations (`scripts/cost.py`, `grid/final/cost_120.txt`, GPU smoke section 10)
- Estimate: dt = 13.79 s (W121 measured) x radial-CFL estimate ratio (9.07 / 11.43) = **10.9 s**; P = 164386 s ->
  15032 cycles/rot (W121 7988); x 120/112 per cycle -> **2.02x WASP-121b per rotation, ~7.5 min/rot on 1 viper node**
  (W121 3.7). Per unit simulated time 1.35x.
- **Measured (apudev job 12527052, vipa1001, same node, same binary athena_dhj_gpu72_18247755_hfresh, 2 MI300A,
  HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1; `gpu/ana.py` -> `gpu/ana_out.txt`):**

  | arm | cycles | rot | s/cycle | mean dt | min/rot (1 node) |
  |---|---|---|---|---|---|
  | TOI-1518b deck | 11000 | 0.595 | 0.0290 | 8.85 s (8.6-9.0) | 9.0 |
  | WASP-121b hfresh deck | 11600 | 1.433 | 0.0276 | 13.61 s (last 13.80) | 3.7 |

  **TOI-1518b costs 2.41x WASP-121b per rotation** (1.05x per cycle, dt 0.65x, rotation 1.49x), **9.0 min/rot on one viper
  node**; per unit simulated time 1.6x. The measured dt is 19 % below the radial-CFL estimate (10.9 s), so the mass-range
  costs of section 2 scale by the same 1.19 (about 1.8x / 2.4x / 3.0x W121 for 1.354 / 1.83 / 2.306 M_J).
  Only the first 0.6 rot were measured; dt may drift as flows develop.
- **Rotations:** thermal relaxation times t_rad = c_p p/(4 sigma g T^3) on the IC columns, TOI / W121: 0.67 (1 bar),
  0.77 (10 bar), 0.79 (100 bar) in seconds; in rotations 0.45-0.53 (the rotation is 1.49x longer). 300 W121 rotations =
  13.7 t_rad(100 bar); TOI-1518b reaches that in ~160 rot.
  **Recommended: 200 rot** (= the W121 300-rot physical time 3.29e7 s, 17 t_rad(100 bar)); check the deep drift at 100 and
  150 rot and stop earlier if the 10-100 bar layer has settled. tlim in the deck = 200 rot.
  Cost: 200 rot ~ 30 h on 1 viper node at the measured 9.0 min/rot (W121 300 rot: 18.5 h), i.e. ~1.6x the WASP-121b
  campaign; on Caltech scale by the H200/MI300A ratio measured there.

## 10. Gates
- **Parse** (`gate/parse.log`, CPU binary `builds/bin/athena_dhj_cpu_8bcd50c8` = fin51-r-fix 182477553 + the dhj-remap-mhd
  merge, default-off MHD remap path): rc 0; nx1 120, x1max 1.968103e10, am_vector, am_transport false, ck_swflux, albedo
  as intended. The 9 unread-key warnings are the same 9 as the W121 hfresh smoke (rad_* route-B leftovers, ic_profile_mu).
- **CPU fresh start** (`gate/gate.sh`, production grid, 12 MPI ranks, ghost zones written, 2 cycles): rc 0, no FATAL,
  **no event of any kind in cycles 0-2** (dhj.log has no rows; the W121 hfresh smoke had 1744 dfloor at cycle 0),
  dt 11.75 s at cycles 0-2. Only warning: the expected am_transport = false note.
- **IC gate** (`gate/gate_ic.py` -> `gate/gate_ic.txt`, t = 0 bin, T and p from (rho, e) through the 3x EOS dump, 84 % of
  active cells inside the dump's rho >= 1e-14):
  - T vs ic_profile: max |dT/T| 1.4e-3 (1e-6..1e-3 bar), 5.6e-3 (1e-3..1 bar), 6.7e-4 (1..300 bar); 2.1e-2 only above
    1e-6 bar (thin top, as KELT-20b's 2 %);
  - deepest active cells 225.0 bar / 4520 K (pole-like) and 414.2 bar / 4931 K (equator-like), on the table adiabat;
    deepest inner ghost 596 bar, rel -2.5e-6; table reaches 2241 bar;
  - **r(p_ref) over 9600 columns 1.337705e10 .. 1.348024e10 vs the tool's pole 1.337697e10 / equator 1.348021e10:
    +6e-6 / +2e-6**;
  - hydrostatic residual |dp/dr + rho g|/(rho g) of the centred difference: median 0.5 % (pole g) to 1.2 % (equator
    g_eff), the size of the O((dr/H)^2) differencing error at ~4 cells/H, not an imbalance (the code is well balanced);
  - min density 1.03e-18 (dfloor 1e-18); top active cells 1.8e-18 .. 5.0e-18, i.e. the top starts at the floor as in
    W121 hfresh (same top pressure by construction).
- **EOS range:** the code's table is log rho -19..0, log T 1.5..6; the IC spans 4.5e3 K (250 bar) to 4100 K at the top
  and rho down to 1e-18: inside. The 3x EOS dump used by the offline tools starts at log rho -14 (~1e-9 bar); above that
  the tools use the exact isothermal continuation.
- **GPU smoke** (apudev job 12527052, the deck verbatim, no overrides, 5.5 min wall): rc 0, no FATAL, 11000 cycles,
  0.595 rot; eos_fail, c2p_it, fofc and vceil all 0. Floor events (dfloor 5.6e8, efloor 9.0e8, tfloor 3.0e8 summed over the
  7 log rows) are about 2x per cycle those of the W121 arm in the same job (2.6e8 / 5.8e8 / 1.8e8 over 11600 cycles),
  i.e. the raised top starts and stays on the floor like W121's.

## 11. What needs a user decision / caveats
- **Metallicity** 3x (adopted, = W121) vs 1x (strict reading of the Fe retrieval): section 3c.
- **Albedo** 0.10 is a choice (TESS implies ~0, the JWST UHJ value is 0.28): range 0-0.3.
- **Teq** 2546 K (Simonnin a/R_*) vs the Seidel check value 2491 K (Cabot a/R_*): +2.2 %.
- **Tide** reaches 10-14 % of g at the domain top: the Roche potential port matters for this planet.
- **Drag:** Simonnin+2025 find tau_drag ~ 1e3-1e4 s; a drag-free hydro run is expected to overshoot the observed Fe winds
  (as for WASP-121b). MHD later.
- **UV below 0.26 um** (4.3 %) is deposited with band-10 opacities.
- No clouds (TOI-1518b is an UHJ; nightside may condense, as for WASP-121b).

## Files
README.md, system.txt, chain.sh, finish_f5.sh (+ _1x variants), MD5SUMS;
deck/ (toi1518b_hfresh.athinput, ovr_3x.json; 1x deck), ic/, ck/CK3X (+ ck/CK1X), sed/ (make_sed.py, atlas9/*.fits,
sw_band_flux_T1518_11.txt, make_sed.out), scripts/ (planet_setup.py = fork/fin51-r-fix 3e0bd80a copy, identical to
182477553; dhj_remap.py; grid112.py; cost.py; system.py; make_ic.py; make_deck.py; make_ovr.py; rce_env*.sh; run_rce*.sh),
rce/ (rceA, rceA2, rceB, rceC = adopted; rce1x*), grid/ (final, prelim_rceA2), sens/ (mass range), gate/, gpu/, lit/.
- The restarts written by the gate and smoke runs (gate/run_fresh*/rst, gpu/run_*/rst; ~2 GB, > 50 MB each) are not part of
  the package and can be deleted. MD5SUMS covers deck/, ic/, ck/*/sw_flux/, scripts/, sed/, the shell scripts, gate and grid
  outputs, rce/rceC, rce/rce1xC, sens/ and gpu/ana*.
