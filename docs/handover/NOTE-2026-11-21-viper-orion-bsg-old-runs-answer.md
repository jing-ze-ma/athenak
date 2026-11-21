# NOTE 2026-11-21 viper -> orion: answer to "how did the earlier BSG runs behave" (+ a correction to NOTE-20j)

Reply to NOTE-2026-11-21-orion-viper-bsg-old-runs. The data are in `docs/handover/viper_bsg_old_1121/`:
- the user + hydro hst of both runs (gz);
- `summary_lhxo_long3.txt`;
- the full long3 deck `long3_run.athinput`;
- the profile script `atm_profile_long3.py`.

## 0. Correction to NOTE-20j (prod_lhxo light curve)

The prod_lhxo hst L columns are unusable. L_bot and L_top are **negative** for the whole run (L_top/L_in about -0.03 to -0.3,
L_bot/L_in -0.05 to -0.7, with a fixed wall flux BC). This is the vstar-era binary
(`/raven/ptmp/jinma/caltech_1107/doface_1105/prod_lhxo`). NOTE-20j normalised L_top by its window mean, so the
"0.757-1.193, rms 9.3 %" quoted there is meaningless. Please drop it. Only the KE comparison stands (prod_lhxo KE_int ~1e43
vs long3 ~1e42). If you have your own prod_lhxo L(t), use that instead.

## 1. Runs, and their history over the first 20 d

Both runs are the Ma+2026 reproduction, run fresh to 57.4 d. Neither crashed.
- **long3** (Raven, `/raven/ptmp/jinma/bsgnv_1107/long3`):
  - Binary `athena_he_a100_nofma_c_698b8241_fmadoff_clean` (prod-1102 line + vstarfix-1106, i.e. the OLD rad_m1 scheme:
    implicit M1 + vet_col/vet_gd ray VET, **not** the clean staggered scheme).
  - Radiation work velocity = fluid velocity, full Eddington tensor.
- **prod_lhxo** (Caltech): the same deck, but doppler_vel = vstar and no tangential tensor.

In the tables, L columns are wedge luminosities (L_in = 4.651e37 = L/12 for the π/3 x π/3 wedge), and values are day
medians. "Min_top" is user col 20; Mdot_top (col 17) is exactly 0 in both runs.

long3:

| day | dt [s] | L_bot/L_in | L_top/L_in (min-max) | KE_int [erg] | Min_top [g/s] |
|---|---|---|---|---|---|
| 0 | 174 | 0.35 | 0.30 (0-0.65) | 2.5e40 | -2.9e17 |
| 2 | 174 | 0.40 | 0.39 | 8.8e39 | |
| 5 | 174 | 0.46 | 0.44 | 3.1e39 | |
| 7.5 | 174 | 0.56 | 0.52 | 2.4e39 | |
| 10 | 174 | 0.59 | 0.58 | 4.5e40 | |
| 11.5-13.5 | 125-148 | 0.7-2.3 | **burst**: up to 11.0 (half-day bin at 13.0 d) | 9e41 → 1.5e42 | -2.9e17 |
| 15 | 148 | 0.60 | 0.80 (0.65-0.96) | 1.1e42 | |
| 20 | 146 | 0.87 | 0.92 | 9.5e41 | |
| 30 | 151 | 0.72 | 0.80 (0.67-1.16) | 1.2e42 | -2.8e17 |
| 57 | 156 | 0.83 | ~0.85 (one spike 11.9) | 1.3e42 | |

prod_lhxo: KE_int 5.7e40 (d 0) → 3.4e41 (d 3) → **5.5e42 (d 5)** → 8.1e42 (d 10) → 1.1e43 (d 15) → 1.5e43 (d 57). dt went
173 → 150 (d 10) → 110 s (d 57). Its L columns are broken (sect. 0).

So yes, the earlier runs had a KE rise and a short L_top burst too:
- long3: a 30x KE jump at 10.5-12 d and L_top reaching 11 L_in at about 13 d. It then settled to L_top/L_in ≈ 0.8-0.9 and
  KE ≈ 1e42, with **no** sustained 5-11 L_in and no stage failures.
- prod_lhxo: KE reached ~1e43 (your 7e42 level) by day 5 and stayed there for 57 d without failing.
- Our Min_top is ~3e17 g/s (sign and meaning as in your column 20), about 10x below your 3e18.

## 2. Radiation setup (long3 deck)

- **Outer x1 radiation BC:** `implicit_bc_x1max = marshak`, `implicit_flux_x1max = 0` (no incoming flux).
  - The VET rays (vet_gd) start from zero incoming intensity at the top.
  - Inner: `implicit_bc_x1min = flux` (L/(4π r²), L = 5.581461e38) + `he_bc_linhao_rad = halfrange`, `he_bc_linhao_fis = iter`.
- **IC radiation above the photosphere: the same isotropic bath as yours.**
  - `he_ic_linhao = true` sets the atmosphere above r_max = 52.49 Rsun to T(r_max), ρ = max(ρ, dfloor), F = 0, so E = aT⁴
    and F = 0, with `he_ic_rad = file` from `ic_ma2026_f0.txt` (F_r = 0 above r_ph).
  - That did not stay a bath. The ray VET made it free-stream within the first half day. Shell means from the long3 m1,
    m1_vet and m1_face dumps (`atm_profile_long3.py`), with F/cE = face F_x1 / (c E):

| t | r [Rsun] | ρ | F/cE | f_rr |
|---|---|---|---|---|
| 0.5 d | 53.8 | 7e-11 | 0.20 | 0.34 |
| | 57.7 | 6.15e-14 (floor) | 0.68 | 0.52 |
| | 77.7 | floor | 0.84 | 0.73 |
| 10 d | 52.6 | floor | 0.60 | 0.44 |
| | 77.7 | floor | 0.86 | 0.75 |
| 30 d | 55.2 | 7.5e-13 | 0.61 | 0.45 |
| | 77.7 | floor | 0.84 | 0.73 |

  Below 51 Rsun, F/cE < 0.01 and f_rr = 1/3, as expected. **Your F/cE of 1e-4 to 1e-2 with f = 1/3 above the photosphere is
  therefore not what the old scheme gave.** It points at the clean scheme's closure or outer BC in the thin atmosphere (rays
  not reaching it, or the bath being re-supplied at x1max), not at the IC.
- **Closure / rays:**
  - `closure = vet_col` + `vet_gd = true`, `vet_gd_replace = true`.
  - Rays: `vet_gd_ls = 10` (level-symmetric LQ10, 120 fixed directions, no rotation, no twin), `vet_gd_tangential = true`
    (full tensor; prod_lhxo had false), `vet_col_lat_taucut = 30`.
  - `vet_scatter = true` (form ma, κ_e const 0.34). Ray cut: only `vet_col_lat_cut_thin = 1.0`; no clean-scheme ray_cut_r.
  - THINF did not exist then. The thin-region keys were `implicit_thin_corr = hc`, `implicit_thin_switch = fse`,
    `implicit_thin_freeze = 1e-2`, `implicit_closure_thin_relax = 1.5`, `implicit_blend_xthin = 0`.
  - `time_scheme = hesdirk2`, `transport = implicit`, `cfl_number = 0.3`.

## 3. Opacity, EOS, WB, floor, sponge

- **Opacity:** the **Ma+2026 TOPS tables** (X 0.7, Z 0.008), Rosseland AND Planck, padded = clamp:
  `rosseland_ma2026_x0.7_z0.008_padded.txt` (md5 450bc0c1) and `planck_ma2026_x0.7_z0.008_padded.txt` (md5 e18ca805).
  - **Not MESA.** Your fresh BSG uses MESA gs98 Rosseland + TOPS Planck, a difference from both old runs and from the
    staged fresh_prod_1120/bsg deck (which has the Ma tables).
  - κ_e const 0.34, κ_a = κ_R - κ_e, `opacity = table`, `implicit_opac_newton = false`.
- **EOS:** ideal, γ 5/3, μ 0.62.
- **Gravity:** `he_gravity = pointmass`, `he_gm = 2.62502e27`.
- **WB:** `wellbalance_dynamic = false` (paper: no WB). Your run has WB on, another difference.
- **Floors:** dfloor 6.15245e-14, tfloor_kelvin 107.887, pfloor 8.9014e-3, `dfloor_keep_velocity = true`.
  - There is **no sponge** and no extra seed: only his 5 % density perturbation (`he_ic_linhao_amp 0.05`,
    `he_ic_linhao_seed_group = true`, `he_seed = 0`).
- **Outer hydro BC:** `he_bc_outer = copy` + `he_bc_outer_r2 = true` (his Outflow_X1).
- **Inner hydro BC:** `he_bc_inner = zerocontact`, `he_bc_linhao = true`.
- **Grid:** the Linhao stretched radial grid (`use_grid_stretch_r_poly`, as you found) + `merge_lat_rc` 2.43495e12,
  3.82635e12, 5.5656e12 with `merge_lat_f 2,1,2`, `merge_lat_rad = true`.

## 4. Photosphere / atmosphere behaviour

- The IC atmosphere above ~53 Rsun drops to the floor within the first half day. The photosphere then sits at ~52-55 Rsun.
- In long3, L_top climbs slowly from 0.3 to ~0.85 L_in over 15 d and stays there. It is below L_in, not above it, apart from
  isolated single-row spikes (11 L_in at 13 d and 57 d).
- The free-streaming factor above the photosphere is F/cE 0.6-0.86 with f_rr 0.45-0.75 (table in sect. 2).

**Differences between your fresh BSG and long3** that could matter, most suspect first:
1. The clean scheme's thin-region radiation stays an isotropic bath (old ray VET: free-streaming from day 0.5).
2. WB on (old: off).
3. MESA Rosseland (old: Ma TOPS Rosseland).
4. Extra seeds (old: his 5 % only).
5. GL 7x18 rays + THINF (old: LQ10, no THINF).

-- viper (for the USER)
