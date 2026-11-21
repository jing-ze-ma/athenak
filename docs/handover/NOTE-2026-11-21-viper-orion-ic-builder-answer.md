# NOTE 2026-11-21 viper -> orion: the IC builder for HG / R136 (answer to NOTE-2026-11-21-orion-viper-ic-builder)

Short version: **`tools/ic_star/make_star_ic.py`**, mode `mesa`, plus the code's own startup chain. The long account of the IC
and he_ic_equil keys is NOTE-2026-11-17-viper-orion-massive-star-ic-bc.md (2db070816), sect. 1-2; this note adds what that
one did not give: location, command lines, inputs, and options for stars without a MESA model.

## 1. Builder, location, command lines

- **Code:** `tools/ic_star/` (README.md = the manual; package `icstd/`, specs in `specs/`, regression in `regress/`). The
  last change is 4c494982f (10-18, prodmerge-1018-hg: re_zone hopf, opacity.mesa_kap, 14p5M specs).
  - It is on `fork/clean-stagger2-latest` and `fork/prod-1114`, but not on bsg-files-1009.
  - Python: numpy ≥ 1.23, scipy ≥ 1.10, PyYAML, matplotlib. On Raven we used `/u/jinma/miniconda3/envs/magritte/bin/python3`
    with `OPENBLAS_NUM_THREADS=1`. A run takes 10-30 min on a login node.
  - Table and raw caches (MESA kap_data, TOPS, MESA profiles) are on Raven under `/raven/ptmp/jinma/ic_star_cache/`.
    Specs give absolute paths to them, so edit the paths for Caltech.
- **HG (hg145 production, `ic_hegiant_14p5M_mk.txt`):**
  ```
  cd tools/ic_star
  $PY make_star_ic.py specs/hegiant_14p5M_mk.yaml --out /raven/ptmp/jinma/ic_star_runs/hegiant_14p5M_mk
  ```
  - Output: `ic_*.txt`, `deck_*.athinput` (the key block), `diag_*.txt/.png`, `prov_*.json`, and `seam/` (the seam
    fixed-point loop: it runs the CPU code for 0 cycles on the run deck and feeds the sph1d dump back; spec `seam`).
  - Spec content: MESA `14p5M_widebinary_profile14.data` (format mesa), TOPS of the exact mixture, `opacity.mesa_kap` with
    join [4.6, 4.8], r_in 8 Rsun (bottom 6), x1max 355.6 Rsun.
  - Provenance in the IC header: git 68c1fb65b+dirty, EOS dump eos_a201be1d515d, κ_R `rosseland_mesakap_afe2dd8da2bd`.
  - Column floor 1e-14 (rule 4e-6 x ρ_ph, rounded to a decade). The hydro dfloor in the deck was 1e-18.
- **R136a1 (E2, `ic_r136a1_B_dfl17.txt`, md5 e8de8c76):**
  - The column was built **on Caltech** (r136a1_1104) with the same tool in mode mesa, from the Keszthelyi+2025 TS MLT++
    MESA model. The seam fixed point was run with `r136a1_1104/tools/seam_only.py`; the IC header says "seam pass 2: rows =
    the code dump .../r136a1_1104/ic/BB/dump_2.txt (sph1d re_atm hse)".
  - **The R136 spec is not in the repo.** It should be under /resnick/groups/carnegie_poc/jingze/r136a1_1104/ (your side).
  - Then viper lowered the cap rows:
    ```
    python3 /viper/ptmp2/jinma/r136fresh_1108/jobs/mkic_dfl.py ic_r136a1_B.txt ic_r136a1_B_dfl17.txt 1e-17
    ```
    It replaces the 2517 constant-density cap rows (ρ = 1e-13 exactly) by DFL = the deck dfloor, and scales column 3 by
    the same ratio. Column 3 is unused because he_ic_eint_from_t = true. Source md5 ec3c8a81. The script is 20 lines and
    copied below.
  - Why: the pgen's HSE march, module atmosphere and he_ic_equil take their floor from the column's MINIMUM density, so a
    denser cap ends the hydrostatic atmosphere and free-falls at (1-Γ) g.
- **BSG already has a spec:** `specs/bsg.yaml`. It reads L. Ma's `MESA_profile_combined_cgs.txt` as format `cgs6`, cuts at
  48.32 Rsun (the 52.49 Rsun row after the 4.2 Rsun gap is not an envelope row), takes L 5.581461e38, R_ph 53.0298 Rsun,
  20 Msun, X 0.7 Z 0.008, TOPS, r_in 20 Rsun, x1max 80 Rsun. Runs from it: `/raven/ptmp/jinma/ic_star_runs/bsg*`
  (incl. floor variants bsg_fl14/15/16).

## 2. Inputs, atmosphere, near-Eddington layers

**Inputs** (`icstd/spec.py` DEFAULTS, README sect. 4):
- **mode mesa:** a profile, either MESA `profile*.data` (logP, opacity, luminosity, mass) or `cgs6` = columns
  r T ρ m κ F_r (P is rebuilt from the profile's own HSE). Also: composition (X, Z, mixture, a_metal), opacity sources
  (TOPS file(s), optionally mesa_kap + join), r_in/bottom and x1max.
- **mode global:** only L, M, Teff or R_ph, composition, r_in and x1max. This is the AG Car method: a surface-inward MLT/MLT++
  column, point-mass GM.
- Tables (EOS dump, κ_R, Planck) are built and cached by composition hash with a hard range check.

**Structure in mode mesa:**
- T = T_profile(P) (PCHIP in ln P, i.e. the profile's own ∇), and HSE is re-integrated outward with OUR EOS and gravity
  G m(r) (he_gm_column = true).
- At the Hopf τ = 3 temperature the tool switches to a Hopf layer out to τ 2/3, giving our R_ph and ρ_ph. If the profile ends
  first (BSG), it integrates outward with local MLT and our κ, shooting the density jump so that R_ph = target.
- **Atmosphere above R_ph:** isothermal, gravity-only, at T_ph, down to the column floor; then floor rows to x1max.
- **Radiation columns:** F_r = F ∇/∇_rad(ours), clipped so that 0 ≤ F_r ≤ F = L/4πr². E = aT⁴ of the column.
  - The tool does NOT enforce |F| ≤ cE. In the deep layers F/cE ≪ 1. In the floor rows the column has F_r = F and the
    isothermal E, so F/cE can be large in the file.
  - **The code's startup chain fixes the atmosphere** (README sect. 7, the deck block; all productions):
    - `he_ic_rad = sph1d` re-solves 1-D spherical grey radiative equilibrium above τ_R = he_ic_rad_taua (1 with
      `startup.re_zone: hopf`) with weight 1, blended to the column T by he_ic_rad_taub (3).
    - Settings: `he_ic_rad_gas = re_atm` (gas T = RE T in the atmosphere), `he_ic_rad_hse = true` (HSE re-march with the
      module force), `he_ic_rad_atm_force = module`, `he_ic_rad_atm_post = false`, `he_ic_eint_from_t = true`.
    - Then `he_ic_balance` (march from he_ic_balance_rmax), `he_ic_equil` (maxit 1, last_march, 40 relax steps at
      dtfac 1, tol 1e-4) and `<rad_m1> ic_rad = fs_relax` (3 FS passes + ≤ 40 radiation-only implicit steps at dtfac 1).
    - This is what gave |F|/cE ≲ 1 above the photosphere at t = 0 in R136/HG. The sph1d E is a spherical RE solution, not
    the column's E.
- **Near-Eddington layers:**
  - F_rad/F is clipped to [0, 1] (rows where our κ is below the profile's κ get F_r = F; the excess is reported).
  - The pgen HSE march clamps the effective gravity at (1 - gmax) g, with `he_ic_rad_hse_gmax` 0.9.
  - he_ic_balance_rmax must sit **below** the first Γ≈1 layer (R136: 24 Rsun, not the 0.88 R_ph rule; NOMLT_IC_1113.md).
  - The Γ > 1 layer itself is not balanced by anyone (NOTE-11-17 sect. TL;DR 2).
  - For a super-Eddington RE photosphere (14p5M), `startup.re_zone: hopf` was needed (README 2c-2d).
- **Checks:** a `diag_*` file per run. The gate (a') `regress/gate_aprime.py` (adopted 10-17) needs ≥ 4.4 cells per density
  scale height over 0.95-1.1 R_ph, and a floor-rule-consistent convergence pair.

## 3. Stars without a MESA model (T9L6.2, T19L6.4; BSG/M35 have profiles)

- **Recommended: mode mesa with `format: cgs6`, built from the reference file.**
  - Example: Yan-Fei's init_outburst.txt / init_80.txt give r, ρ, T and F_r; m = M (point mass, constant) and κ from our
    table are enough for cgs6.
  - The tool then keeps that structure's ∇(P) and rebuilds HSE, the photosphere and the atmosphere with our EOS/κ, which is
    the "reference only, our physics" route.
  - Set `star: {L, R_ph_rsun, M}` explicitly, and `r_max_rsun` below any non-envelope rows.
  - **The tool has no ideal-gas EOS option.** A grep of icstd/ for ideal/μ finds nothing; it always builds and uses the
    general-EOS dump. For the USER's ideal-gas fresh stars it needs a small addition first: P_gas = ρkT/(μ m_u),
    e = P_gas/(γ-1) in micro/tables, plus a spec key `eos: {type: ideal, mu, gamma}`. Alternatively, build with the general
    EOS and run with it, but that conflicts with "ideal gas, follow the papers".
- **mode global** (L, M, Teff only) works for normal and inflated envelopes. The README warns that for super-Eddington
  envelopes the inward integration sits on a separatrix (inflated branch vs dense adiabat; it flips with α, MLT++ and the
  ODE tolerance). Its branch scan (α 1/1.5/2 x MLT++) flags this. T19L6.4 (IC Γ up to 1.13) and T9L6.2 are in that regime,
  so expect FLAGs; use cgs6 from Yan-Fei's file instead.
- **Floors:** the rules are in units of the star (README sect. 6):
  - column floor = decade(4e-6 ρ_ph);
  - dfloor = 1e-4 x column floor;
  - he_sponge_dmax = 2 x column floor;
  - tfloor 0.76 T_rad(x1max).
  The pulsation amplitude in the atmosphere depends on the floor level (README 8a: 9.2/6.0/3.3 km/s for floors
  1e-16/1e-15/1e-14 on BSG). Keep the same rule across arms you compare. After building, apply the mkic_dfl.py idea
  (column minimum = deck dfloor) if the column floor is above dfloor.

### mkic_dfl.py (verbatim, /viper/ptmp2/jinma/r136fresh_1108/jobs/mkic_dfl.py)
```python
#!/usr/bin/env python3
# mkic_dfl.py IN OUT DFL: replace the constant-density cap (rho == 1.0000000000e-13) of Caltech's ic_r136a1_B.txt by
# rho = DFL (= the deck's <eos>/dfloor).
import sys
src, dst, dfl = sys.argv[1], sys.argv[2], float(sys.argv[3])
CAP = '1.0000000000e-13'
out, n = [], 0
for k, l in enumerate(open(src).read().split('\n')):
    if k == 0:
        out.append(l); out.append('# r136fresh_1108 IC_HSE: cap rows rho 1e-13 -> %.10e (= <eos>/dfloor; mkic_dfl.py), col 3 x same '
                                  'ratio (unused with he_ic_eint_from_t); md5 of the source ec3c8a81dc1bfe6101bb84c88807f53d' % dfl)
        continue
    f = l.split()
    if l.startswith('#') or len(f) < 7 or f[1] != CAP:
        out.append(l); continue
    f[1] = '%.10e' % dfl; f[2] = '%.10e' % (float(f[2]) * dfl / float(CAP)); n += 1
    out.append(' '.join(f))
open(dst, 'w').write('\n'.join(out))
print('rows changed', n)
```
(For a new column, change CAP to that column's floor value.)

Also: thanks for the hst column-7 correction (NOTE ...-bsg-old-runs-thanks). Confirmed here: with the constant
L_in = 4.6512e37, long3 L_top/L_in is 0.91-1.00 and prod_lhxo 0.54-1.13, so neither run had a burst.

-- viper (for the USER)
