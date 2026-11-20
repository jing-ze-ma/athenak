# fresh_prod_1120: viper inputs for the fresh v20.1 production runs (reply to NOTE-2026-11-20h)

USER decisions (11-20, PST):
- **Ideal gas for all stars, following the papers.** This is an exception to the general-EOS rule.
- **Orion runs these on Caltech.**

Every μ, GM and L below was checked against its IC; the check script is `ic_consistency_check.py`. The μ check uses
total-pressure hydrostatic balance only, so the opacity plays no part in it. The L check is the radiative-diffusion
luminosity of the IC with the deck's Rosseland table.

Data files are gzipped; run `gunzip -r .` first. In the decks, replace the `FRESH_DIR` placeholder with the absolute
path of this folder. The BSG deck instead points at the Caltech originals in `/resnick/.../bsgli_1014/files/`; the
copies here have the same content. Checksums are in `MD5SUMS`; they are of the gzipped files.

| star | deck | IC / profile | EOS | GM | L | opacity | grid | layout run here |
|---|---|---|---|---|---|---|---|---|
| BSG (Ma+2026 reproduction) | `bsg/bsg.athinput` | Linhao profile `linhao_profile_full.txt` (he_ic_linhao) + `ic_ma2026_f0.txt` | ideal, γ 5/3, **μ 0.62 (his)** | 2.62502e27 (19.78 Msun, his deck; USER 11-20 keep) | 5.581461e38 | Ma+2026 tables (padded); dfloor 6.15245e-14 (his, keep) | sp 20-80 Rsun, 256³, π/3 x π/3, SMR merge_lat_rc | 16 blocks, 4 Raven nodes x 4 A100 |
| M35MMS (Will, Schultz+22 T35L5.2) | `m35mms/m35mms.athinput` (paper deck) | `ic_m35mms.txt` (radiative-only HSE, Jiang+18 construction) + `prof_m35mms.txt` (inner ghosts); MESA `35M_save_profile_MS.data` | ideal, γ 5/3, μ 0.6215 | 35 Msun point mass | 1.063065e39 (MESA) | MESA Rosseland (`tables/rosseland_mesakap...`) + TOPS Planck (`tables/planck_ext2...`) | sp 9.697 Rsun-1.114 R_ph, 336x256x256, 20°x20° | 4 blocks, 1 node |
| M35ZAMS (Will, Schultz+22 T42L5.0), NEW | **no deck yet**: grid from `GRID_REFERENCE_m35zams_ours_mesaR_1blk.athinput` (a general-EOS deck: take the mesh/meshblock only); EOS/IC/table/BC keys as in m35mms.athinput | `ic_m35zams_ideal.txt` + `prof_m35zams_ideal.txt`; MESA `35M_save_profile_ZAMS.data`; builder `build_ic_ideal.py` | ideal, γ 5/3, μ 0.6220 | 35 Msun | 6.61828e38 | same as M35MMS | sp 6.798-8.200 Rsun, 384x128x128, θ 1.53875-1.60284, φ 0-0.064 | - |
| T19L6.4 (Yan-Fei, Jiang+18) | `t19l64/t19l64.athinput` | `ic_t19l64_LIC.txt` = `init_80.txt` as is (+ ghost extension below 14.83 Rsun, `build_ic.py`) | ideal, γ 5/3, μ 0.674 | 1.0268751734e28 (77.356 Msun) | 4.02933e39 (log 6.022) | same as M35MMS | sp 14.83-336.5 Rsun, 512x256x512, θ π/4-3π/4, φ 0-π | 8 blocks, 4 viper nodes |

## Consistency checks (IC vs deck)

- **T19L6.4:**
  - Total-pressure HSE fit of init_80: μ 0.6742, GM 77.35 Msun.
  - Diffusion L: log 6.020-6.023 at every radius from 14.9 to 28 Rsun, across the Fe bump.
  - Table 1 says log L 6.40, but that is impossible for this IC: even pure electron scattering gives at most 6.16 in
    the deep layers. Following the IC file is deliberate.
- **M35MMS / M35ZAMS:** exact by construction.
  - HSE 5e-8 / 4e-8, and the diffusion flux equals L.
  - ZAMS: Γ_Edd peaks at 0.88 (Fe bump), there is no density inversion, and ρ is within 1-4 % of MESA.
  - MMS: Γ_Edd reaches 1.04 at the Fe bump (radiative-only construction).
- **BSG (USER 11-20: use HIS GM, μ and floor):**
  - μ 0.62 (Ma+2026 text / his deck), GM 19.78 Msun (his deck), dfloor 6.15245e-14 (his).
  - For information only: the Linhao profile by itself balances at μ 0.614 with its own M(r) (19.82-19.99 Msun).
    With his values, the HSE mismatch is 0.7 % (μ) plus 0.2-1.1 % (GM), the same as in his own run.
- **All decks:** `temp_unit_kelvin` = μ m_u/k_B to within 6e-8.

## T19L6.4 corona test (`t19l64/corona_test/`): please run this on Caltech before any long T19/T9 run

The viper T19 run of 11-19 hit trouble from 30 h on. It used the same deck as `t19l64.athinput`, the binary
athena_he_gpu72_x_72562926_t9fix and 4 MI300A.
- **implicit_wb_vkick_reject:** the threshold is 3e8 cm/s. From cycle 3305 there were 4462 backward-Euler solves with
  kicks (49k cells), 379 stage rejects, and two dt collapses to 0.4-0.5 s at 38 h and 42 h.
- **Where it happens** (from the 3-h dumps; `kick_locations.py`):
  - The bulk flow at 1000-3000 km/s is dense (ρ ~1e-10), outward and above v_esc, at 20-30 Rsun. It looks like a
    physical super-Eddington Fe-bump outflow (the IC's Γ reaches 1.13).
  - The cells at the 3e8 cap are near-empty (ρ 5e-19-1e-14) and hot (2-7e8 K). They sit at 37-45 Rsun, just above
    the photosphere (31 Rsun), and by 44.6 h also in evacuated pockets inside the envelope (23-25 Rsun).
  - It is not the T9L6.2-type far interface, where the envelope met infalling floor gas at 108-116 Rsun
    (viper DTDROP_1115).

The test deck `t19l64_corona_test.athinput` changes three things, tagged [CORONA 11-20], plus basename t19cor:
1. **IC `ic_t19l64_corona.txt`** (`build_ic_corona.py`):
   - Above r_j = 31.635 Rsun the floor gas is replaced by a static hydrostatic corona: ρ = 2e-16 (r/x1max)^-2, so
     2.3e-14 → 2e-16.
   - P comes from HSE with gravity GM (1 - κ_R L/(4π c GM)), with Γ_R 0.35. T runs from 8.2e6 to 7.7e5 K.
   - Residuals: HSE 5e-7. Pressure is continuous at r_j, with a static density jump of 122.
   - Mass 1.5e-8 Msun; sound speed ≤ 4e7 cm/s, which allows about 330 s against the run's 38 s.
   - The F_r and E columns are unchanged, so only the gas changes.
2. **dfloor 5e-19 → 1e-16.**
3. **he_sponge_dmax 1e-18 → 1e-15:** v_r is damped in every near-empty cell.

Suggested run: a fresh start to t = 1.8e5 s (50 h). Compare with the baseline above: kick counts and capped cells per
dump. Caveats:
- With dmax 1e-15 the sponge also damps the outer corona beyond about 150 Rsun, along with any outflow thinner than
  1e-15 there.
- The radiation energy above the photosphere is the IC's constant E, which is 4-400 x F/c. If kicks remain, that
  reservoir is the next suspect.

## T9L6.2 first (USER 11-20): `t9l62/` and `t9l62/corona_test/`

USER: try T9L6.2 first. T9L6.2 is the cool (T_eff 9000 K) outbursting LBV model of Jiang+18; its IC is
`init_outburst.txt`.

**Deck:** `t9l62/t9l62.athinput` is the viper our-wall deck (production closed inner wall). Grid: sp 24.23-809.8 Rsun,
512x256x512, θ π/4-3π/4, φ 0-π. 8 blocks on 4 viper nodes.

**IC:** `ic_t9l62_LIC.txt` = `init_outburst.txt` as is, plus the ghost extension (`build_ic.py`).
- ideal gas, γ 5/3, μ 0.674, GM 7.4033817458e27 (55.77 Msun); L_IC 2.47875e39 (log 5.811).
- Total-pressure HSE fit: μ 0.6738-0.6742 and GM 55.76 Msun.
- Diffusion L: log 5.809-5.814 at 24-88 Rsun. Table 1 says 6.20, which is impossible for this IC (at most 5.97 even
  with pure electron scattering).

**Viper history** (setup and diagnosis in `T9L62_SETUP_1113.md` and `DTDROP_1115.md`):
- The paper-BC run (deck t9l62_final) died at ~189 h.
  - The envelope expanded at +1.9e7 cm/s into floor gas that had been free-falling since t = 0 (-2.3e7 at
    107-118 Rsun).
  - At that interface (108-116 Rsun) the implicit M1 solve kicked up to 2.9e4 cells beyond c, producing
    T 1e11-1e13 K and a dt collapse to 0.2 s.
- The our-wall run (this deck) reached 7.54e5 s (~210 h) at dt 12 s. It did not crash; the USER stopped it on 11-19.
- The dt fell from 99 s to 14 s at 175-189 h. That part looks physical: radiation-dominated voids at 34-38 Rsun once
  convection breaks the density inversion, as Jiang+18 describe, about 3x later than in the paper.

**Corona test** (`t9l62/corona_test/t9l62_corona_test.athinput`; changes tagged [CORONA 11-20], basename t9cor). It
applies the same recipe as the T19 test above.
1. **IC `ic_t9l62_corona.txt`:**
   - Above r_j = 95.49 Rsun (just above the IC's τ = 1 at 93.9 Rsun) the floor gas is replaced by a static hydrostatic
     corona: ρ 1.44e-14 → 2e-16 at x1max, T 2.1e6 → 2.5e5 K, Γ_R 0.30-0.33.
   - Residuals: HSE 8e-7. Pressure is continuous, with a density jump of 50.
   - Mass 2e-7 Msun; sound speed ≤ 2.1e7 cm/s, so no dt cost.
2. **dfloor 3.6e-19 → 1e-16.**
3. **New floor sponge, as in T19:** he_sponge_rate 1e-3, mode radial, r0 0.9 x1max, dmax 1e-15. The T9 deck had no
   sponge before.

Suggested: a fresh start to at least t = 8e5 s (~220 h), past the 189 h crash window, on the same Caltech
layout as T19. Compare dt and the vkick/NOT ADMISSIBLE counts with the baseline deck (or with the viper our-wall
run, whose run.log and dumps are in /viper/ptmp2/jinma/t9l62_1113/prod_viper, on request).

## Not included

- **Will's M13 MMS, 24.75 Msun EMS and 65 Msun MMS:** their grids are unknown (waiting for Will).
