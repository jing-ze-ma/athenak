# NOTE unify-1012: one rad_m1 method for vet_gd (sp wedge) and vet_sc (Cartesian)

Workers W (audit, first code) and W2 (validation, final code), 2026-10-12, MPCDF Raven. Branch `unify-1012`
(local, not pushed), worktree /raven/u/jinma/ATHENAK/wt-unify, work dir /raven/ptmp/jinma/unify_1012.
History: W's WIP on 0197f020, rebased onto prodfinal-1011 **4253af30** (conservative2 + implicit_cs2_hll = drag),
merged **vgdlat-1012 669b5db9** (opt-in vet_gd_symref = col), then the W2 changes. Final binaries
`bin/athena_cpu_u5` (md5 72b7916d...) and `bin/athena_he_cpu_u5` (e237b0a7...); u4 = u5 minus the vet_gd lin gate
(bitwise equal to u5 on every arm below that names implicit_thin_bq = cell on vet_gd, and on all vet_sc arms).

## 0. Result: the unified setting

**U = `implicit_thin_switch fse`, `implicit_thin_corr hc`, `implicit_thin_recon jse`, `implicit_thin_bq lin`,
`implicit_thin_hcsf j`** (+ the rest of the production set: blend idort_f, xthin 0, hesdirk2, predict, vstage, vimp,
vet_bc_diff face, vet_seg emit, conservative2 + drag, cs2_shift false, vimp_bc true; vet_gd_symref = col on vet_gd).
All of these are now the **fresh-run defaults** on both meshes (restarts whose file lacks a key keep the old value,
bitwise). Checked: a fresh Cartesian cb run with recon / bq / hcsf / time2_vet_sc / time2_vstage unnamed is bitwise
the named U; a fresh sp patm run with recon / bq / js_gd unnamed is bitwise the named U (u5, chk job 31091117).

Two **geometry differences remain, both numerical (not physical)**, each because the vet_gd formal-solution
data are not as good as vet_sc's in that one place; both are keys, so they can be switched on for tests:
1. `implicit_thin_js_gd = x1` (fresh default): the js/jse face shape is applied on the radial faces of vet_gd
   only. vet_gd's LATERAL half-range data come from the per-shell lateral interpolation of the GD formal solution
   (the T6a/T6b weakness); the lateral face shape amplifies its beam-edge errors: P2 pinhole beam L/L_ex 3.5 and
   7-9 NONCONV with lateral js on, 1.5 and 1-2 without. vet_gd_symref = col does NOT cure it (3.46-3.57, 7 NC).
   vet_sc's lateral data are the full 2-D short-characteristics solution, so there jse acts on every face.
2. `implicit_thin_bq_gd = cell` (default; lin opt-in): thin_bq lin (outer Marshak q from the face-extrapolated
   H+) acts on vet_sc only. On vet_gd W's r^2 H+ extrapolation puts a +1.3 % kink into the top-cell E (AG Car A:
   top-cell E/J-1 1.6e-2 vs 2.7e-3, E r^2 smooth without it; T6a thin_u max error 2.4e-3 vs 5.0e-4 at 256, patm
   thin-top max 1.0e-3 vs 5.1e-4) and buys nothing. On vet_sc lin is THE fix of the grey-atm thin top (below).

`implicit_thin_corr hcs` is not used: it is inert on production vet_gd (eps = 1) and worse on Cartesian hcyl
(E out 4.4 % vs 1.4 %). `implicit_thin_hcsf = jh` (W's cbg fix) is **rejected** (sect. 2).

## 1. Audit: every method key whose default or effect differs between the two paths (code at 0197f020)

| key | vet_sc (Cartesian) | vet_gd (sp wedge) | where |
|---|---|---|---|
| implicit_thin_switch / implicit_thin_corr | fresh default fse + hc (same code path) | fresh default fse + hc | rad_m1_implicit.cpp ~2027-2056 |
| implicit_thin_corr = hcs | gate 1 - eps active in scattering media | inert: eps = min(kappa_P/kappa_T,1) = 1 above the vet_gd cut (AG Car, BSG) | M1ThinFace |
| implicit_thin_recon = js / jse | active | **silently inert** (`tjs = ... && hsrc == 1`) | ~3589 (lateral), ~11663 (x1) |
| implicit_thin_bq = lin | active (outer Marshak q from the face-extrapolated H+) | **silently inert** (`bql ... hsrc == 1`); vet_gd's hr_q keeps the cell value h+_N | ~11349 |
| implicit_hr_face = sf | default, active | inert (`impl_hr_sf = vet_sc`; no SC face-flux block for vet_gd) | ~1927 |
| implicit_hr_excess = iso | not available (needs vgd_hr J) | opt-in | ~1932 |
| time2_vet_sc / time2_vet_col | time2_vet_sc default **lag** (read only when named) | time2_vet_col default **predict** | rad_m1_time2.cpp ~173 / ~225 |
| time2_vstage | default **false** unless named | default **true** (fresh) | rad_m1_time2.cpp ~216 |
| time2_pred_cap | 0.5 whenever predict | 0.5 | same |
| vet_bc_diff, vet_seg | face, emit (fresh) | face, emit (fresh) | unified already |
| vet_bc_bath | vet_sc key (bath top) | n/a (vet_col_surface_q) | boundary-specific, physically different |
| implicit_vimp, time_scheme | hesdirk2 + vimp (fresh, gas feedback) | same | unified |
| implicit_vimp_bc | default false (read only when named) | default false | unified in code; suiteF PROD names true, AG Car deck does not |
| implicit_gas_src | default matched | default matched | unified in code; suiteF PROD uses conservative2, AG Car deck matched (G / consfix) |

## 2. cbg verdict and the Jh fix (W's implicit_thin_hcsf = jh)

cbg = Gauss 8x32 crossing beams, 256x64, c dt/dx 256 (suiteF FAIL: PROD FATAL at cycle 17). Runs runs2/cbg (u1),
runs4/cbg (final).

| variant | rc | Picard / NONCONV | row -1.9 / 1.9, peak |
|---|---|---|---|
| hc + dc + cell, Jh (P with W's fix) | FATAL (Picard diverged, cycle ~5) | - | - |
| hc + dc + cell, j (= suiteF PROD) | FATAL (cycle 17) | - | - |
| hc + dc + lin, Jh | FATAL | - | - |
| hc + jse + cell, Jh | 0 | 8.95 / 0 | 0.96 / 0.84, 0.073 |
| hc + js + lin, Jh | 0 | 8.92 / 0 | 0.96 / 0.84, 0.073 |
| hc + jse + lin, Jh | 0 | 3.92 / 0 | 0.96 / 0.84, 0.073 |
| **hc + jse + lin, j (U final, runs4)** | **0** | 34.1 / 6 | 0.89 / 0.82, 0.071 |
| hcs + jse + lin (R) | 0 | 2.92 / 0 | 0.91 / 0.83, 0.072 |

- **The cbg rescue is the face-shape reconstruction (js / jse), not Jh.** With dc the run dies with or without Jh.
- **Jh breaks vacuum crossing beams:** every hc run with Jh (cb r4, r12, c dt/dx 0.5 / 32, cbh, for P and U alike)
  has row energy 1.164 (+16 % over exact) and peak 0.581; without Jh (and with hcs) 1.006 / 0.417 = PROD. Jh only
  lowers U's cbg Picard count (34 -> 3.9). Default is therefore j; jh stays an opt-in key.
- **Verdict: cbg PASSES with U** (completes, row 0.82-0.89; the Picard cost 34 / 6 NONCONV is the known
  vacuum-beam convergence weakness of hc, as cb r4: PROD 23.2, U 28.5, R 3.1).

## 3. Per-test table, final code (u4/u5), U vs P (production thin pair hc/dc/cell) vs R (hcs/jse/lin, u1)

### 3a. Cartesian (vet_sc), runs4 (U, P) / runs2 (R); arms_final.sh (conservative2 + drag, hcsf j named)

| test | U (final) | P (hc/dc/cell) | R (hcs/jse/lin) |
|---|---|---|---|
| cb ring 4 / 12 / c 0.5 / c 32: row, peak, Picard | 1.006, 0.417; Pic 28.5 / 3.1 / 2.9 / 3.9 | 1.006, 0.417; 23.2 / 3.1 / 3.0 / 10.2 | 1.006, 0.417; 3.1 / 26.1 / 2.9 / 12.7 |
| cb 1024x256 (r4n4x) | 0.998, peak 0.737, row 1.9 0.965 | 0.998, 0.737, 0.938 | 0.998, 0.737, 0.965 |
| cbh ring 4 | 1.006 / 0.417, Pic 4.0 | same, 4.1 | same, 3.9 |
| **cbg** | **completes** (sect. 2) | **FATAL** | completes |
| static shadow ring 12 / 36: L1, wedge (ex 157.6), umbra | 0.011, 165.5 / 165.4, 0.22 / 0.19 | 0.011, 162.0 / 161.9, 0.35 / 0.33 | 0.011, 165.3, 0.22 / 0.19 |
| static shadow Gauss 8x32 | 0.030, umbra 0.21 | 0.030, 0.33 | 0.030, 0.21 |
| homogeneous cylinder Gauss / 512: E out, F out | **1.4 % / 2.2 %**, 2.1 % / 3.9 % | 2.3 % / 2.8 %, 3.7 % / 4.8 % | 4.4 % / 3.9 %, 5.7 % / 6.0 % |
| bubbles g8x32k3 / r100k3 / r100k5, L1 at 1 t_c | 0.756 / 0.755 / 0.642 | 0.792 / 0.805 / 0.677 | 0.756 / - / - |
| thermal eq 1 / 2, eq3t, non-LTE 1e-1..1e-8, static diffusion, dyn. diffusion ddf/ddiff/ddn, Marshak c1/c10, shocks M1.2/2/3, crooked pipe c2200 | identical to suiteF PROD to the printed digits (E 96.79505 / 145.78780; eq3t 1.132601, P_tot 3.000000; nlte L1 1.11 %..0.01 %; ddf D_eff 0.984; Marshak 0.0137 / 0.0141; pipe probes 0.938 ... 0.556) | (suiteF) | identical |
| grey atm ladder 32-2048 (be), self-conv. by tau band | **2.00 in every band**; tau<0.1 e 4.8e-7 at 2048 | tau<0.1 ord 0.7/0.6/1.5/0.95/0.6, tau 0.1-1 ord drops to 0.13 | 2.00 (= U; lin alone = U) |
| rw_p10t1 X / C | 2.0-2.36 / 2.05-2.26 (= P = R bitwise digits) | same | same |
| pulse k 0.128 / 12.8 / 1280 X, C (static gas) | finest-pair order E 2.14 / 2.24 / 1.93; = P = R bitwise digits | same | same |
| crooked pipe c220 | reaches t 10 in the 4 h job: probes 0.938 / 0.855 / 0.773 / 0.688 / 0.554 (= REC 0.938 ... 0.554) | suiteF PROD reached only t 4.9 in its job (0.929 ... 0.514 there) | - |

### 3b. Spherical wedge (vet_gd), runs_f (u4): PF = HFfWEdRdBcUc (P + symref + drag), UF = HFfWEdRjBcGxUc (final U)

| stage | test | UF | PF | note |
|---|---|---|---|---|
| P1 | patm L16 poly: E/J-1 at R_ph, tau<0.1 L1 | -0.47 / -0.23 / -0.10 %, 2.1e-3 / 6.3e-4 / 1.9e-4 (ord 1.74 1.76) | identical | symref: thin top now converges (suiteF -0.67) |
| P1T | patmT coupled | **-0.11 / -0.14 / -0.08 %**; tau<0.1 1.4e-3 / 4.7e-4 / 1.5e-4; tau 0.1-1 2.5e-3 / 1.3e-3 / 5.6e-4 | -0.53 / -0.24 / -0.10 %; 2.1e-3 / 6.2e-4 / 1.8e-4; 6.1e-3 / 2.2e-3 / 7.4e-4 | U better |
| P2 | pinc thermal beam L/L_ex | 1.49 / 1.52 / 1.54, Picard 7.3, 2 NONCONV | 1.51 / 1.55 / 1.58, Picard 5.0, 0 NONCONV | known-open excess; lateral js: 3.46-3.57, 7 NC |
| T6a | thin_u / thin_L16_u E L1 | 8.1e-4 / 2.0e-4 / 4.8e-5 (ord 2.01 2.08) | identical | (bq lin on vet_gd: max 2.35e-3 vs 4.96e-4) |
| T6b | core_u r>3 E/J-1, F r^2/H r^2 at 64/128/256 | +1.3 / -2.4 / +1.6 %, 1.012 / 0.976 / 1.015 | +5.3 / +0.2 / +3.3 %, 1.061 / 1.007 / 1.034 | U better at 64, 256 |
| T4 | dif_r1 | L1 4.6e-4 at 256, ord 2.01 2.00 | identical | |
| T7 | shd 48x48: umbra, penumbra L1 | 0.051, 0.52 | 0.049, 0.50 | same within ray noise |
| T5 | marsh_r1_u X | dens 1.09, E 0.52, F1 0.56; F1 L1 3.8e-3 | dens 1.09, E 0.55, F1 0.66; F1 4.4e-3 | front, 1st order |
| T9 | rsw X / T | 1.83-2.04 / 2.00-2.08 | identical | |
| new | rlat X / T | 1.92-2.14 / 1.80-1.83 | identical | |
| new | rlatt X (thin) | dens 2.45, E 1.40, F1 0.61 | identical | known-open thin lateral F1 |

### 3c. AG Car A reduced (CPU, 480x32x32, 60 cycles, he binary), agc/runs, deck agcA_F.athinput (symref, cons2 + drag)

| arm | L_top/L_in end (min) | Picard mean / max / NC / BE redo | E/J-1 at R_ph | top-cell E/J-1 (r > 1.1 R_ph max) |
|---|---|---|---|---|
| P (production deck, matched, no symref, u2) | 0.99996 (0.99939) | 3.016 / 6 / 0 / 0 | +2.46e-3 | 3.5e-2 |
| PF (P thin + symref + cons2 drag, u4) | 0.99997 (0.99956) | 2.872 / 6 / 0 / 0 | +2.46e-3 | 2.7e-3 |
| **FJ = final U (jse x1, bq cell on vet_gd)** | 0.99997 (0.99958) | **2.864** / 6 / 0 / 0 | +2.50e-3 | **2.6e-3** |
| UF (jse x1 + lin on vet_gd) | 0.99997 (0.99958) | 2.872 / 6 / 0 / 0 | +2.50e-3 | 1.6e-2 (top-cell kink) |
| UA (jse all + lin) | 0.99997 (0.99958) | 2.872 / 6 / 0 / 0 | +2.49e-3 | 1.6e-2 |
| FL (lin only) | 0.99997 (0.99956) | 2.872 / 6 / 0 / 0 | +2.46e-3 | 1.6e-2 |

Production AG Car is neutral to the thin set (L_top equal to 1e-5, Picard 2.86-2.87, 0 NONCONV); the final U is the
best arm on every column. (Earlier u2 run: U with lin 3.008, the same top-cell kink.)

## 4. Negative E (cbg, cb) -- diagnosis

Runs runs3 (jobs_neg.txt, u1, `dbg_t2_admiss 4`), plus every runs2/runs4 log:
- The linear solve returns negative E (`min E from the solve` -3.8e-3 cbg, -2.4e-3 cb) in EVERY arm with
  `implicit_offdiag = operator` (the default): thin corr none / hc / hcs, vet_tensor full or uniaxial, gas source
  matched or conservative2, vimp_bc on or off, be or hesdirk2. `implicit_offdiag = none`: 0 positivity fallbacks,
  0 non-admissible hesdirk2 stages, and identical rows / peak / L1 (cbg none: 0.922 / 0.847, peak 0.055).
- So the source is the off-diagonal (cross-derivative) part of the implicit operator, which is not an M-matrix in
  the vacuum next to a beam; the positivity fallback (that cell re-solved with offdiag none) handles it and the
  result is unchanged. The hesdirk2 non-admissible stages (14-19 per cbg run) are the explicit stage vector
  E^n + dt a K1 going to -7e-13 in dark cells (E ~ 1e-14, `E_old_vector <= 0` in 1898 cells) -> BE redo.
- Pre-existing (cbfix be runs: min E -4.1e-3 / -2.3e-3), not tied to the thin set, benign. No change made.

## 5. Code (branch unify-1012 on 4253af30 + vgdlat-1012)

- `implicit_thin_recon` / `implicit_thin_bq`: fresh defaults jse / lin on every mesh (restart without key: dc / cell).
- `implicit_thin_js_gd = x1 | lat | all | none`: js/jse on vet_gd per direction; fresh default x1, restart none
  (= 0197f020, where js was inert on vet_gd).
- `implicit_thin_bq_gd = cell | lin`: the face-extrapolated q on vet_gd (r^2 H+), opt-in; default cell.
- `implicit_thin_hcsf = j | jh`: default j (bitwise older runs); jh opt-in (sect. 2).
- `time2_vet_sc` default predict and `time2_vstage` default true on fresh runs on every mesh (W; as time2_vet_col /
  the wedge). Restarts without the keys keep lag / false.
- vet_gd_symref stays opt-in in code (vgdlat); production decks should name `vet_gd_symref = col`.
- Production deck keys for U (all named): `implicit_thin_switch fse`, `implicit_thin_corr hc`,
  `implicit_thin_recon jse`, `implicit_thin_bq lin`, `implicit_thin_hcsf j`, `implicit_thin_js_gd x1`,
  `implicit_thin_bq_gd cell`, `vet_gd_symref col`, `implicit_gas_src conservative2`, `implicit_cs2_hll drag`,
  `implicit_cs2_shift false`, `implicit_vimp_bc true` (+ the hesdirk2 set). For the running AG Car A/B chains
  (restarts, deck names dc / cell): switching to U changes only the vet_gd radial faces' thin shape (jse x1);
  AG Car A is neutral to it (sect. 3c).

## 6. Open / not done
- rw ladders (rw_t10a, rw_p10t1 T) were not rerun with drag (u1 ran X/C only: = PROD); the drag change is
  4253af30's own validation.
- vet_gd lin: a correct spherical extrapolation of the outer q might help the vet_gd thin top too; not pursued.
- P2 pinc 2 NONCONV with U (P 0 with symref). Small.

## 7. Files
- Scripts: arms.sh / arms_final.sh / arms_def.sh, run1.sh (ARMS_FILE), pool.sbatch, poolB/C/D.sbatch, chk.sh,
  build.sh; jobs_A / jobs_E / jobs_final_cart / jobs_neg / jobs_D; sp/jobs_sp / jobs_final / jobs_final2 / jobs_F;
  sp/ana/spl.py (suffixes Rj Rd Bl Bc Zs Gx Gl Gn Ga Ed Uc Hj), sp/ana/evalU.sh, evalFinal.sh;
  ana/evalU.sh, ana/evalFinal.sh; agc/agc.sbatch (arms P U R N PF UF UA FL FJ), agc/inp/agcA_F.athinput,
  agc/ana/agceval.py; ladder/ana/od.py (Ub Rb Lb Pb, PRODp Up Rp).
- Results: ana/out/unifyU_cart.txt (u1, Jh), ana/out/final_cart.txt; sp/ana/out/unifyU.txt (u1),
  sp/ana/out/final.txt; logs/lad.txt, logs/lad_pulse.txt.
- Runs: runs2 (u1, Jh default), runs3 (negative E), runs4 (final Cartesian), runs5 (bitwise / defaults checks);
  sp/runs (u1/u3), sp/runs_f (u4 final), sp/runs_chk (u5 checks); agc/runs; ladder/run/runs.

