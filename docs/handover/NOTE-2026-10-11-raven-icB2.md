# NOTE icB2 (worker K, 2026-10-11): AG Car B IC with MLT++ continued to tau_R = 2/3

Work dir /raven/ptmp/jinma/icB2_1011. Generator committed on branch `icB2-1011` (from prodfinal-1011 df483770,
worktree /raven/u/jinma/ATHENAK/wt-icB2), commit 7955f73e: docs/handover/scripts/make_ic_agcar_mltpp.py. Not pushed.

## 1. Generator reconstructed on Raven (step 1)

make_ic_agcar_ge.py exists only on Viper. make_ic_agcar_mltpp.py rebuilds it on top of make_ic_mlt_star.py,
using the same inputs:
- EOS dump agcar_1008/shakeB_ana/eos_dump.txt (X 0.36, Y 0.62, Saha H/He). eint is recovered to 8e-11.
- TOPS Rosseland (the file header's table; md5 49aea217). It equals ext2 to 3 digits at T > 20 kK.
- L, Teff, M: log L 6.17, Teff 20 kK, M 55 Msun, alpha_MLT 1.5.

What had to be recovered from the file itself:
- **MLT++.** In every convective row, (nabla - nabla_ad)/(nabla_MLT - nabla_ad) = 1.000e-3 +- 1e-5, which is
  f2, wherever the local-MLT excess nabla_MLT - nabla_ad > f1 = 1e-4. Rows with a smaller excess are untouched.
  The switch is sharp, between excesses 9.68e-5 and 1.004e-4, at the CZ base, rows 1926/1927.
- **Grid.**
  - Interior: dln r = 2e-4, but not anchored at r_join, and base is inserted. The script takes the radii from the
    old file (`--grid_from`).
  - Hopf layer: r_join + k DRN, with DRN = 2e-4 R_ph, plus the R_ph row.
  - Atmosphere: R_ph + k DRN, with the floor 1e-16 every 50 DRN out to 350 Rsun.
- **Iso atmosphere.** Gravity only (atm_rad 0).
- **E above R_ph.** E = (F/c)(1 + (c a T_ph^4/F_ph - 1)(R_ph/r)^2). The sph1d IC replaces this.

Reproduction: `--mode inward --grid_from ic_agcar_B_ge.txt` gives ic_agcar_B_ge_repro.txt (md5 9dc151a6). Same
5822 rows, r to 1.4e-11. Max relative difference per column (cmp.py):

| region | rho | eint | F_r (col 4) | E | T | F_MLT/F (col 7, abs) |
|---|---|---|---|---|---|---|
| r < R_ph | 3.6e-7 | 4.2e-7 | 1.9e-7 | 3.2e-7 | 8.0e-8 | 1.2e-7 |
| r >= R_ph | 6.1e-9 | 5.6e-9 | 0 | 0 | 0 | 0 |

The join point is the same: r_j 6.9935333883e12 (tau 3), T 25930.42 K, rho 8.0789232e-11, against the file's
8.0789232070e-11.

## 2. New column (step 2): `--mode outward`

- Rows r < r_j are copied **verbatim** from ic_agcar_B_ge.txt (5419 rows; `cmp` shows them identical).
- From (r_j, T_j, rho_j), the same MLT++ equations are integrated **outward** (dPtot/dr = -rho g,
  dlnT/dlnP = nabla_MLT++, log grid continued) until the old column's tau scale reaches 2/3. The tau scale is
  the Hopf layer's, with tau(r_j) = 3. MLT++ therefore replaces the grey Hopf layer over tau 3 -> 2/3, as asked.
- Above R_ph': the same isothermal gravity-only atmosphere at T_ph' = T(R_ph'), starting from rho(R_ph').

**File: /raven/ptmp/jinma/agcar_1008/files/ic_agcar_B_ge_mltpp23.txt, md5 c6b8b932134f8d119f7235d3d4ec7d0e**
(5845 lines, also in MD5SUMS.icB2). The original ic_agcar_B_ge.txt (d39c9afa) is untouched.

| quantity | old (Hopf to tau 3) | new (MLT++ to tau 2/3) |
|---|---|---|
| R_ph (generator) | 7.04718e12 | 7.08089e12 (+0.48 %) |
| R_ph (tau 2/3 from the top, ext2, file) | 7.04721e12 | 7.05738e12 (+0.14 %) |
| R_ph after the sph1d HSE (run log) | 7.12966e12 | 7.10872e12 |
| T_ph | 20113 K | 20788 K (Teff' 19952 K) |
| Gamma_rad = kappa F_r/(c g), 6e12 < r < R_ph | max 1.337 (14 cells > 1) | **max 0.8655** |
| Gamma_rad, tau 3..2/3 | 0.97-1.34 | <= 0.8655 |
| density inversions below R_ph | 39 rows | **0** (and 0 above) |
| F_MLT/F at r_j-, r_j, r_j+ | 0.1139, 0 (jump) | 0.1139, 0.1119, 0.1093 (continuous) |
| F_MLT/F, tau 3 -> R_ph' | 0 | 0.13 -> 0.029 (tau 1.9) -> 0.047 at R_ph', then 0 |
| L: F_r/(1 - col 7) / F - 1 | 2e-7 | 2e-7 (col-7 print precision) |

- Gamma_tot (all flux radiative) in the new layer is 0.89-0.99. The density no longer rises outward, so kappa
  stays lower.
- MLT++ holds Gamma_rad at about 4 nabla_ad (1 - beta), about 0.865.

**Caveat (tau bookkeeping).** A gravity-only cap on (rho, T)(R_ph') has its own tau of only 0.19. Even on r_j
it has only 0.48.
- The file's true tau_R at R_ph' is therefore 0.19. In the file, MLT++ reaches true tau 0.19-0.24, and the true
  tau at r_j is 2.52.
- The run rebuilds everything above tau 3000 anyway (sph1d re_atm + HSE). The code's R_ph after the HSE is
  7.1087e12, outside R_ph'.

Variant (`--stop atm --cap_rad 1`): **ic_agcar_B_ge_mltpp23r.txt**, md5 d0e4e257.
- R_ph' is placed where a cap with g - kappa F/c has tau 2/3, which is the true tau_R 2/3. MLT++ extends to
  7.11749e12 (+1.0 %), T_ph' 18652 K, and the true tau at r_j is 3.59.
- Kept only as a comparison.

**How the code uses col 7.** mlt_flux_frozen sets F_MLT = F_r fm/(1 - fm) on the faces (col 4 and col 7,
interpolated). esrc = -div(A F_MLT) w(t) is deposited into the radiation (he_star_m1.cpp:1769-1868), with the
same time-only ramp w(t) (mlt_ramp_start/time) on every face. The new nonzero F_MLT above tau 3 is therefore
ramped off with the rest of the scaffold.
- The step at the top of the MLT layer moves from tau 3 (0.112 L) to R_ph' (0.047 L).
- The run log line is unchanged: max F_MLT/F on the faces is 0.7589 at 5.116e12, and F_MLT at the top face is 0.
- mlt_closure_taper is used only with mlt_closure = adaptive.

## 3. Validation, reduced CPU (step 4; job 31085267, binary athena_he_cpu_df483770 md5 27e5cf67)

Setup:
- 480x32x32, 4 ranks x 3 threads, 60 cycles, production dt.
- Decks inp/agcar_B_vh_{new,newr,old}.athinput: copies of prod/agcar_B_fresh with he_ic_file set, hse true,
  gmax/damp/maxit 0.9/0.5/400 named, outputs every 10 cycles.
- Analysis: ana/dyn.py (E's script); window 0.97-1.03 x 7.04718e12; ana/dyn_31085267.txt, ana/prof.py.

| arm | sph1d HSE | clamped cells | R_ph IC | Picard mean/max | NONCONV | max\|v\| c10/c30/c60 [km/s] | rho drift c10/c30/c60 | L_top c60 |
|---|---|---|---|---|---|---|---|---|
| Bold_h (old column = E) | 146 outer, 49 s | 538 | 7.04718 -> 7.12966e12 | 2.73 / 8 | 0 | 0.22 / 0.53 / 0.82 | 0.003 / 0.013 / 0.026 | 1.0005 |
| **Bnew_h (mltpp23)** | **327 outer / 2274 VET, 111 s, converged (dln rho 1.0e-6)** | **0** | 7.05736 -> 7.10872e12 | **2.78 / 8** | **0** | **0.11 / 0.31 / 0.59** | **0.0004 / 0.004 / 0.016** | 1.0002 |
| Bnewr_h (mltpp23r) | 285 outer, 89 s | 0 | 7.11751 -> 7.10318e12 | 2.76 / 7 | 0 | 0.11 / 0.32 / 0.61 | 0.0005 / 0.004 / 0.017 | 1.0053 |
| Bnew_n (hse false) | - | - | 7.05736 -> 7.06094e12 | 2.72 / 8 | 0 | 3.2 / 7.9 / 13.6 | 0.07 / 0.35 / 1.62 | 0.9990 |
| Bnew_f (he_ic_rad file) | - | - | file | 2.68 / 7 | 0 | 3.3 / 7.9 / 13.6 | 0.07 / 0.35 / 1.64 | 0.9990 |

- Bold_h on df483770 reproduces E's eef1f8f6 numbers digit for digit.
- **Where the motion sits (c60, shell means, ana/prof.py):**
  - The old column's outflow from its clamped layer (<v_r> 0.64 km/s at 0.994 R) is **gone**.
  - Inside 0.97-1.02 R: <v_r> <= 0.06 km/s, max|v| <= 0.35 km/s (the seed level), |drho| <= 0.13 %.
  - The remaining window maximum (0.59 km/s, 1.6 %) is at 1.025-1.05 R in the near-floor atmosphere (rho
    1e-12..1e-16). It is the same there for the old column (0.55-0.63 km/s): the floor feature, as in A.
- **sph1d vs column at the photosphere.** T/T_col is 1.014 at R_ph and rho/rho_col 0.96. With the old column
  these were 1.196 and 0.50. The new column is close to the code's own RE + HSE solution.
- Without hse (n, f), B still expands, at 13.6 km/s (old column: 19). hse = true stays required for B.

## 4. A100 full size (gpudev job 31085275, athena_he_a100_df483770 md5 41db0c68)

Full size, 4 A100, 40 cycles, deck inp/agcar_B_vh_new.athinput; runs/j.31085275/Bnew_h:
- rc 0, FATAL 0, **NON-CONVERGED 0, Picard 2.65 / 8**.
- The sph1d solve is identical to the CPU one: 327 outers, 2274 VET passes, **0 clamped cells**, R_ph 7.05736 ->
  7.10872e12. It takes 101 s on the GPU node host.
- Photosphere window (0.97-1.03 x 7.04718e12), c10 / c20 / c30 / c40:
  - max|v| 0.106 / 0.208 / 0.308 / 0.405 km/s;
  - rho drift 0.0004 / 0.0017 / 0.0039 / 0.0069. This matches the CPU run (0.107 / 0.311 at c10 / c30).
- L_top c1 / c40: 1.0000 / 1.0002. 40 cycles take 155 s wall including the IC.
- ana/dyn_gpu_31085275.txt.
- The DeltaAI copy (task_deltaai/) is WITHDRAWN, because Raven finished first.

## 5. Production package (step 5; NOT submitted)

- prod/agcar_B_fresh_v2.athinput (md5 023c8fe79ecd242037f7b8d34cb1af3b).
  - Made from agcar_B_fresh.athinput by mkdeck.py, with he_ic_file -> ic_agcar_B_ge_mltpp23.txt and
    he_ic_rad = sph1d, he_ic_rad_gas = re_atm, **he_ic_rad_hse = true** named.
  - The three hse keys are named: he_ic_rad_hse_gmax 0.9, he_ic_rad_hse_damp 0.5, he_ic_rad_hse_maxit 400.
  - implicit_det_reduce, implicit_vimp_bc and time2_start are not named, so they take their defaults (vimp_bc
    and slope start off).
- **Marked line, pending the coordinator's confirmation:** `he_sponge_dmax = 2.0e-16` (was 1.0e-17).
  - Reason: worker J found that the IC floor atmosphere (rho 1e-16) falls in while dmax is below the floor.
  - The user approved this value on the condition that J's CPU job 31085684 passes.
  - If J's test does not pass, revert the line.
- `implicit_gas_src = matched` is left in the deck; the coordinator switches it after the consfix-1011 merge.
- prod/run.cfg.fresh_B_v2: BIN athena_he_a100_df483770, BIN_MD5 41db0c6875c16244ae9808efb32a8365.

## Open
- MLT++ now reaches the generator's photosphere. The code's photosphere after the HSE (7.109e12) lies above it,
  in a layer with F_MLT = 0. The 0.047 L step at R_ph' is deposited at true tau of about 0.2-0.3 until the ramp
  removes it.
- The deck change to v2 is the user's decision (old B deck: hse false).
