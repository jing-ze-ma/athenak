# NOTE coreL-1012: the M1 luminosity of an opaque emitting surface (T6b), `vet_seg = jump` and `implicit_thin_hcface = tau`

Worker Y, 2026-10-12. Local branches, not pushed. Commit identity Jingze Ma.

| branch | worktree | based on | commits |
|---|---|---|---|
| `coreL-1012` | /raven/u/jinma/ATHENAK/wt-coreL | vgdlat-1012 669b5db9 | b5cb6c8c, d487c142 |
| `coreLu-1012` | /raven/u/jinma/ATHENAK/wt-coreLu | unify-1012 7986b1d8 (W2 WIP, which includes vgdlat) | the same two commits cherry-picked: 0f101630, 369efa9c. The hcsf / js_gd conflicts in rad_m1_implicit.cpp are resolved by keeping both changes. |

- Work dir: /raven/ptmp/jinma/coreL_1012.
- Binaries in bin/:
  - j4 = b5cb6c8c, plus lateral tau on the sp wedge (the behaviour before d487c142);
  - j5 = d487c142;
  - u1 / he_u1 = coreLu 0f101630;
  - u2 / he_u2 = 369efa9c.
- In both worktrees `kokkos` is a symlink to wt-prodfinal/kokkos. It is not committed.

## 0. The reference was off

- s3_eval's L_ex is H r^2 from a GL mu quadrature at r about 7. That quadrature is noisy by +-0.2 % (the discontinuous limb).
- The code's core is the cell set with centre r < 1. Its edge is the face 0.996875 at N 128 / 256 / 512 / 1024, and 0.94 or 1.05 at N 64.
- All numbers below are therefore against L_ex(R_code) = 0.5 ∫_0^R (1 - e^{-2 kappa sqrt(R^2 - p^2)}) p dp, which is exact. At N 128-1024 this is 0.9937 L_ex(1).
- Script: ana/lsum.py.

## 1. Diagnosis (4x4 lateral-trivial core ladder; ana/edge.py, ana/fsconv*.py, ana/symref2.py)

**Base result** (HFfMY1 = production keys + symref; vet_seg emit; fse + hc). L_M1/L_ex(R) at N 64 / 128 / 256 / 512 / 1024 is -3.7 / +1.4 / +4.2 / +4.9 / +2.9 %. It does not converge.

**The steady budget is exact.** The integrated emission ∫ chi c (B - E) dV equals the exterior F r^2 to 1e-4. So the luminosity excess is an E deficit of 1-4 % in the 5 cells below the surface, at N 512, in the thermalisation layer where 1 - E is O(0.1-0.5).

**Those cells have a correct closure.**
- They use the central flux (tw = 0), and D_rr matches the exact f_K to 1e-3.
- E is pulled down from the switch face: w is about 1 at the face between the edge cell and the one below it. There hc imposes the SC's face flux, c (H_L + H_R)/2.
- The SC's H in the edge cells is +5 % and its J -9 % (N 512).

**Cause 1: the SC segment rule at a material discontinuity (`vet_seg = emit`, rad_m1_vetseg.hpp).**
1. Inward rays. Entering the emitting body, the face S on the opaque side starts from the thin cell's S (about 0). S then ramps linearly to S_0 over the own half cell, so I- inside is too low: H is too high and J too low.
2. Outward rays.
   - The emit limiter's weight th carries the factor (1 - e^{-hu}). Once the edge cell is optically thin (tau_cell 0.13 at N 512), the face S falls back to the mean ramp S_u -> S_0, about 0.06 at N 512.
   - Half of the outer half-cell's emission is lost: the FS L is -6 % at N 512.
3. The 1-D symref reproduces the code's FS: edge J -8.7 %, H +5-6 %, FS L -5.7 % (N 512).
4. The FS errors fall at order 1 or less, but the hc face mean turns them into an M1 luminosity error that is non-monotone in N.

**Cause 2: the hc face flux.** hc's face flux Hf = (H_L + H_R)/2 is the arithmetic mean across the kink of H at the surface: H is steep inside (it changes over 1/chi) and constant outside (r^2 H = const). This is first order low: with an exact FS the M1 L is still -3 % at N 512.

## 2. Fixes (opt-in, defaults bitwise)

**`<rad_m1>/vet_seg = jump`.** It is emit plus the following (rad_m1_vetseg.hpp, mode 2):
- **Outward.** S_f = S_u + hu*(gf + th*((1 - e^{-hu}) mm - gf)), with th = (hu - h0)/(hu + h0), the chi contrast alone.
  - The opaque-cell slope estimate gb = S_u - I_u keeps the factor 1 - e^{-hu}, because it is only valid in an opaque cell.
  - A thin edge cell of an emitting body therefore keeps its own S up to the face.
- **Inward.** Where the own cell is the thicker and the hotter (h0 > hu and S_0 > S_u), its part of the segment starts from S_f0 = S_f + (h0 - hu)/(h0 + hu)(S_0 - S_f). S may jump at the face.
- **The emit restriction is kept.** A cold opaque cell heated from its thin side keeps the mean face value: the heated walls are unchanged.
- **Turned ray (vet_gd sweep typ 3, and the C++ symref).**
  - When the ray dips below the own inner face, the chord there takes chi and S of the cell below. These are the log-mean of the own column and the lateral upwind point at i-1. The rule is piecewise constant per cell.
  - This is the limb of the body inside the edge cell, which the shell-centred sweep otherwise omits.
- **Scope.** It applies to vet_sc, vet_gd, vet_col_lat and symref through the shared VetSegSplit. Smooth fields see O(dx) contrasts, so second order is kept.
- **Effect on the 1-D FS.** The interior FS becomes exact (<= 2e-3). The FS L is -1.8 / -0.9 / -0.4 / -0.16 % at N 64-512 (first order); emit gives -19 / -14 / -10 / -6 %.

**`<rad_m1>/implicit_thin_hcface = tau`** (rad_m1_implicit.cpp, M1ThinFace and its x1 / transverse call sites):
- The formal-solution face flux of thin_corr (hc / hcs / lag) is interpolated linearly in optical depth: Hf = (tau_R H_L + tau_L H_R)/(tau_L + tau_R).
  - tau is the half-cell depth, kappa_T times the distance from the centre to the face.
  - On the sp x1 faces it is applied to r^2 H, which is exact for a free-streaming exterior.
- Next to a thin cell, Hf is the flux the surface emits. With equal depths it is the arithmetic mean, so smooth fields stay second order.
- **`implicit_thin_hcface_gd = x1 (default) | all`** decides whether the sp wedge's transverse faces take it too.
  - Cartesian meshes always take every face.
  - On P2, `all` gives the same channel L, with Picard 5.5 against 4.7 for x1.

## 3. T6b convergence (4x4, `core_u`, L_M1/L_ex(R_code) - 1)

| arm | N 64 | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|---|
| base HFfMY1 (vgdlat binary = j-binaries default, bitwise) | -3.7 % | +1.4 % | +4.2 % | +4.9 % | +2.9 % |
| jump (Sj) | -19 % | -12 % | -6.6 % | -3.0 % | -0.6 % |
| tau (Ht) | | +3.9 % | +4.9 % | +5.0 % | |
| no hc (C0) + jump | | -0.9 % | +0.1 % | +1.1 % | +2.0 % |
| **jump + tau (SjHt)** | **-3.5 %** | **-1.6 %** | **-0.58 %** | **-0.14 %** | **+0.02 %** |
| jump + tau, no symref (HFfMSjHt) | | -0.7 % | -0.2 % | +0.07 % | |
| unify branch (u1) base | -8.1 % | -1.8 % | +2.3 % | +3.9 % | +2.6 % |
| unify branch (u1) jump + tau | -3.6 % | -1.6 % | -0.58 % | -0.15 % | +0.01 % |

- The jump + tau order is 1.2 / 1.5 / 2.0 / about 3.
- Exterior E (r > 3) is -3.6 / -1.6 / -0.6 / -0.15 / +0.02 %, and the core interior is <= 5e-4.
- vet_bc_diff cell = face (bitwise). Before W2's branch, implicit_thin_recon jse / bq lin was inert on vet_gd (bitwise).
- 16x16 (core_L16_u): F r^2 / L_ex(R) is 0.983 / 0.993 at N 128 / 256, against 1.013 / 1.041 for base.

## 4. Regression

**vet_gd** (j4 / j5 against the vgdlat sr4 runs, i.e. V's runs4):

| test | base | jump + tau |
|---|---|---|
| P1 patm L16, E/J-1 at R_ph, N 64 / 128 / 256 / 512 | -4.7e-3 / -2.3e-3 / -1.0e-3 / -3.8e-4 | -2.0e-4 / -5.8e-4 / -3.0e-4 / -9.8e-5 |
| P1 tau < 0.1 band L1 (N 512) | 5.4e-5 | 1.1e-5 |
| P1 envelope (N 512) | -5.9e-5 | -1.2e-5 |
| P1 Picard | 2.0 | 2.0 |
| T6a thin L16, E L1, N 128 / 256 | 2.0e-4 / 4.8e-5 | 4.0e-6 / 3.9e-6 (floor) |
| T6a max \|F r^2 - 1\| | | same |
| T4 dif r0 3 | 1.8e-3 / 4.6e-4, order 2.0 | identical |
| T7 shd 48x48 | lit L1 0.51 | identical to 1e-3 |
| L pulse, r0/w 1000 vs Cartesian | 3.753e-4 | 3.735e-4 |
| L pulse, r0/w 1 | 8.86e-2 | 8.83e-2 |
| rsw (K2), X order | 2.09 | 2.05 (1e-10 level) |
| P2 pinc, L/L_ex r 4.3 / 4.6 / 4.9 | 1.51 / 1.55 / 1.58, Picard 5.0, 0 NC | 1.77 / 1.82 / 1.85, Picard 4.7, 1 NC |
| P2 core luminosity (r 1.6-1.95) / exact (R 1.5, chi 50, tau_cell 3.1) | 0.815 | 0.951 |
| P2 channel / core luminosity | 1.85 | 1.86 |

- **P2 is not a regression of the channel.** The base hid 18 % of a missing core luminosity, and the known moment-step channel excess (P2b) is unchanged.
- The P2 ablation: jump alone gives core 0.73 and channel 1.35-1.41; tau alone gives core 1.017 and channel 1.89-1.97.

**Cartesian vet_sc** (REC = tJe + hesdirk2; emit against jump + tau; cart/, binary j4):

| test | emit | jump + tau |
|---|---|---|
| slab, J_sc first thin cell, tau_cell 0.1 / 1 / 10 / 100 / 1e4 | 0.469 / 0.449 / 0.4996 / 0.5000 / 0.5000 | 0.496 / 0.5000 / 0.5000 / 0.5000 / 0.5000 |
| slab, M1 E/0.5 - 1, tau_cell 0.1 / 1 / 10 / 100 | +2.0 / +7.2 / +2.6 / +0.26 % | -0.9 / -2.2 / -2.4 / -0.24 % |
| hcyl rho kappa 10, L1(E) out | 0.044 | 0.016 |
| hcyl rho kappa 10, L1(F) out | 0.057 | 0.012 |
| hcyl rho kappa 10, rF/rF_ex | 1.052 | 0.992 |
| hcyl rho kappa 1000 (tau_cell 31), L1(E) | 0.218 | 0.230 |
| hcyl rho kappa 1000, rF/rF_ex | 0.703 | 0.691 |

- hcyl at rho kappa 1000: the residual is the moment-side defect at an unresolved surface (segfix section 2).
- **Identical:** heated wall wall1d (f1, f4 against the f64 reference), crossing beams (cb r12, cbh r4, cbg g8x32), static shadow (shd r12g2, shdg g8x32), eq1, sdifk, and non-LTE nltek1-8.
- Bubbles: L1 0.252 → 0.253.
- marsh: harness FATAL in both arms (vimp on a 1-D mesh, as in segfix).
- **Crooked pipe** (full, c dt/dx 220, t = 10):
  - T_end at the 5 probes is identical to 1e-3: 0.938 / 0.855 / 0.773 / 0.689 → 0.688 / 0.554.
  - t(T >= 0.5) is 0.10 / 0.10 / 0.30 / 0.50 / 3.90 → 0.10 / 0.10 / 0.30 / 0.70 / 4.00.
  - So probes 4 and 5 arrive one or two output intervals later, against the large emit/mean differences of segfix (1.9 vs 0.3).

**Reduced CPU AG Car A** (W2's agc setup: 480 x 32 x 32, 4 MPI x 18 OMP, 300 cycles; agcA_F deck = W2 final; arm J = + vet_seg jump + hcface tau; he_u1; agc/):

| quantity | F | J |
|---|---|---|
| L_top/L_in | 1.00022 | 1.00022 |
| min / max over the run | 0.99956 / 1.00022 | 0.99965 / 1.00022 |
| Picard | 3.215, 0 NC | 3.215, 0 NC |
| E/J-1 at R_ph | +2.4e-3 | -1.8e-4 |
| max \|E/J-1\| over 0.9-1.1 R_ph | 1.77 % | 1.62 % |
| max \|E/J-1\| above 1.1 R_ph | 2.7e-3 | 3.1e-3 |
| CPU | | +0.3 % |

- Rerun with he_u2 (hcface_gd = x1 default): identical numbers.
- The F arm is bitwise between u1 and u2, and the CPU cost is the same.

## 5. Status, recommendation

- Both keys are opt-in. Defaults are bitwise to sr4: the thin_L16, patm_L16, core_u and shd outputs are identical.
- The T6b target is met: L_core -> exact at order up to 2-3, and the M1 core luminosity no longer exceeds the black-body limit.
- **Candidate for default on fresh runs** (vet_sc and vet_gd alike): `vet_seg = jump` together with `implicit_thin_hcface = tau`. Neither helps alone.
- Before a production default:
  - a GPU A100 cost / determinism check;
  - a longer AG Car A segment and a BSG gate.
- **Open:**
  - The P2 channel excess (about 1.85x the source) is the moment step, as in P2b.
  - The unresolved opaque surface (tau_cell >> 1: hcyl rho kappa 1000, P2 core at -5 %) is still first order / moment-side.
