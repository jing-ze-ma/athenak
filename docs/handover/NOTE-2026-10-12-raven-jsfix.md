# NOTE jsfix-1012: the E pile-up in front of a sharp thin -> opaque step (implicit_thin_recon = jse)

Worker JS, 2026-10-12, MPCDF Raven. Branch **jsfix-1012** (local, not pushed), from prodfinal-1011 ed32f0fa.
Worktree /raven/u/jinma/ATHENAK/wt-jsfix (kokkos symlinked, not committed). Commit **8bfeae8c**.
Binaries: bin/athena_cpu_j1 (md5 728038c3...), bin/athena_he_cpu_j1 (74a17df8...); build job 31106449.
Main merged 8bfeae8c into bsgcand-1012 (d3bf42d7).

## 1. Cause

M1ThinJs (rad_m1_implicit.cpp) builds the js/jse face shape from the face value of the formal-solution half-range
fluxes, H+-_f. That value was the arithmetic mean of the two cell-centre values:

    h+_f = h+_L * ph,  ph = clamp(0.5*(1 + H+_R/H+_L), 1/jslim, jslim)

At a sharp thin (L) -> opaque cold absorber (R) step, H+_R at the absorber's centre is attenuated by e^{-tau_R/2},
so it is about 0. Then ph = 0.5, the donor's upwind ratio h+ is halved and the thin cell must hold about 2x E
to push its flux into the absorber. The true face value of H+ (the intensity is continuous along the ray) is the
thin side's value.

- **vet_sc has the same defect.** The Cartesian static shadow hides it because its edge is a smooth tanh, so
  neighbouring tau differ only by O(dx). A sharp 1-D slab shows it (sect. 3). So the cause is not anything
  specific to vet_gd: not the lateral interpolation, the twin/symref, eps or r^2.
- **The eps weighting is not the cause.** jse uses wt = eps, and eps = 1 for a cold absorber, so the full
  shape acts there.
- **The limiter only caps the damage.** jslim 2 bounds the halving at 0.5.

## 2. Fix (one rule on every mesh and face): `implicit_thin_js_face = hcface | mean`

- **hcface (fresh-run default).** The js/jse face value of H+- takes the SAME face interpolation as the hc
  correction's face flux, `implicit_thin_hcface`:
  - **tau:** optical-depth weights H_f = (tau_R H_L + tau_L H_R)/(tau_L + tau_R), times (r/r_f)^2 on the sp x1
    faces. This is coreL's rule.
  - **mean:** the old arithmetic mean.
- **So js and hc now see one face value.** At a material step it is the thin side's value. With equal depths it
  is exactly the arithmetic mean, so homogeneous Cartesian runs are bitwise and smooth fields stay 2nd order.
- **mean** is the old behaviour. A restart whose file lacks the key keeps mean, which is bitwise to ed32f0fa
  (checked: T7 with `js_face=mean` on j1 is identical to ed32f0fa on every metric).
- **Code.**
  - M1ThinJs takes hwl/hwr. The hcface weight blocks are moved before the js call at both call sites: the
    x1 kernel and the lateral-face kernel (where js_gd lets js act).
  - The key is read in the thin-switch setup next to implicit_thin_jslim.
  - `thsw_jsmean` is declared in rad_m1.hpp.
- **Scope.** vet_gd keeps js on x1 only (js_gd x1), and its lateral faces keep the arithmetic mean
  (hcface_gd x1), as before. vet_sc: every face.

## 3. Reproduction and ladders (work dir runs: cart/runs, sp/runs; ana/lad.py, ana/t7j.py)

**T7 spherical shadow** (48x48x32, FIN keys, 120 directions; front = the thin cells radially in front of the clump,
exact E/E_free = 1):

| run | front E/E_free mean / max | umbra | lit L1 r 4.05 / 4.99 / 5.77 |
|---|---|---|---|
| ed32f0fa (B) | 1.56 / 1.84 | 0.144 | 0.329 / 0.482 / 0.507 |
| j1 js_face mean (M) | 1.56 / 1.84 (= B) | 0.144 | same |
| **j1 default hcface (J)** | **0.92 / 1.00** | 0.138 | 0.329 / 0.482 / 0.507 |
| recon dc (Q2 diag) | 0.88 / 0.96 | 0.134 | 0.329 / 0.479 / 0.502 |
| PRODd | 0.88 / 0.96 | 0.049 | 0.373 / 0.520 / 0.564 |

**1-D Cartesian slab (vet_sc, FIN keys).**
- Setup: hot opaque slab x < 0 (B = 1), a thin gap 0 < x < 0.5 where the exact E is 0.5, and a cold opaque
  absorber 0.5 < x < 1.
- Reported: E/0.5 of the FRONT cell (the last thin cell), and the gap L1.

| nx1 (absorber tau_cell 100 at 64) | B (= M) front, L1 | J front, L1 | D (recon dc) front |
|---|---|---|---|
| 32 | 1.976, 1.3e-1 | 0.988, 1.3e-2 | 0.988 |
| 64 | 1.995, 6.7e-2 | 0.998, 4.8e-3 | 0.998 |
| 128 | 2.005, 3.3e-2 | 1.003, 1.5e-3 | 1.003 |
| 256 | 2.000, 1.6e-2 | 1.000, 2.5e-4 | 1.000 |
| 512 | 1.999, 7.9e-3 | 1.000, 1.1e-4 | 1.000 |

- **Absorber tau_cell sweep at 64** (0.1 / 1 / 10 / 1000 / 1e5), front cell:
  - B: 0.996 / 1.041 / 1.986 / 1.991 / 1.990.
  - J: 0.995 / 0.995 / 0.995 / 0.995 / 0.995, equal to D.
- Picard: 2.05-2.13 in every arm.
- So B is a 0th-order defect (2x at every N) for tau_cell >~ 5. J converges and equals the dc answer.

## 4. Validation, arm FIN on S's harness (suiteF_1012: runs_js, sp/runs_js, sp/ana/out_js, ladder/run_js)

| group | result vs final_1012 FIN (78b2c961) |
|---|---|
| Cartesian: eq1/2, eq3t, ddf/ddn/ddiff, sdifk, rs_m1.2 | bitwise |
| nlte k1-8, shocks M2/M3 | ~1e-10 relative |
| cb r4/r12/c0.5/c32/n4x, cbh, cbg | rows, peak, Picard identical (cbg 0.891/0.829, Picard 34.1) |
| hcyl Gauss / 512 | E out 1.78 / 1.37 %, rF 0.988 / 0.993 (identical to 4 digits) |
| shd r12 / r36, shdg | wedge 161.6 -> 160.9 (exact 157.6), L1 0.011, umbra 0.22 / 0.20; shdg L1 0.030 |
| bubbles r100k5 / r100k3 / g8x32k3 at 1 t_c | L1 0.687 / 0.757 / 0.764 -> 0.714 / 0.762 / 0.769; E_tot 0.56 -> 0.60 (known FAIL, small shift) |
| Marshak c1/c10, cloud | unchanged |
| pipe c2200 | probes 0.938 ... 0.557 (0.556) |
| pipe c220 (t 10) | probes 0.938 / 0.856 / 0.773 / 0.689 / 0.555 (was 0.938 / 0.855 / 0.773 / 0.689 / 0.554); t(T >= 0.5) 0.1 / 0.2 / 0.4 / 0.6 / 3.81 (was 0.1 / 0.1 / 0.3 / 0.7 / 4.0) |
| grey atm (steady ladder) | identical errors and orders in every tau / Kn band (tau 0.1-1: 1.88 2.06 2.08) |
| pulse k 0.128 / 12.8 / 1280 X, C | identical (finest-pair E orders 2.14 / 2.24 / 1.93) |
| sp: P1 patm, patmT (+eps 1e-2/1e-4), T6a thin, T4 dif, T9 rsw/rswt/rswv/rswtv X/T, rlat/rlatt X/T, rswt_q, L pulse/nlte2/rsw | identical to the printed digits (patm E/J-1 at R_ph +2.1e-4 / -5.9e-4 / -2.9e-4) |
| T6b core_u, core_L16_u | r > 3: +6.40 / -2.22 / -1.23 % (was +6.34 / -2.23 / -1.23); F r^2/H_ex 1.062 / 0.977 / 0.987 |
| T5 marsh_r1 | X/T orders unchanged (lateral-component noise only) |
| P2 pinc | L/L_ex 0.61-0.63 -> 0.57-0.59 (PRODd keys at rotate 0: 0.54-0.57; known-open moment-step channel); NC 1 -> 0 |
| **AG Car A reduced CPU** (300 cyc) | L_top/L_in 1.00022 (min 0.99965), Picard 3.149 / 6 / 0 NC, E/J-1 at R_ph -1.61e-4 (-1.60e-4): identical |
| **BSG reduced** (bsgR, 30 cyc) | rc 0, Picard 2.797 / 9 / 0 NC (identical), L_top equal to 4e-9 |

- **Verdict:** PASS. The default changes only where js acts across an optical-depth jump. It fixes the T7 front
  and the slab, it is neutral on every gain W2 and Y documented, and AG Car and BSG are unchanged. Not run:
  the A100 binary (no GPU-specific code), and the Cartesian X/T/C ladders and rw (homogeneous media, bitwise by
  construction).

## 5. NEW finding (OPEN, not caused by js): the M1 cavity inside a complete cold opaque spherical shell

- **Setup.** sp/inp/sshell.athinput is a 1-D spherical shell (4x4 wedge): star flux 1 at r = 1, thin r < 2.5
  (core_r), and a cold opaque gas outside (rho kappa 1e3), with the FIN keys.
- **Result.** E/E_free in the cavity is far off and the same in EVERY arm: ed32f0fa, j1, recon dc,
  thin_switch none, thin_corr none, vet_seg emit, and vet_gd = false (vet_col alone).

  | N | first / mid / front cell E/E_free |
  |---|---|
  | 32 | 28 / 114 / 211 |
  | 64 | 13 / 52 / 109 |
  | 128 | 5.9 / 27 / 55 |
  | 256 | 2.7 / 13 / 28 |

  - The error falls as 1/N (E roughly proportional to dx) and the field is steady (nlim 200 gives the same).
  - With rho kappa 10 (tau_cell about 0.8) E is too LOW near the star: 0.23 at r = 1, rising to 1.9 at the front.
  - The no-shell control is exact to 1e-3.
- **Not seen elsewhere.** In T7 the clump subtends a small solid angle and the energy goes around it. The
  Cartesian slab with the same step is exact.
- **Hypothesis, not verified.**
  - The sp face thin | opaque is not reached by the upwind / thin-switch path. The switch keys change nothing,
    so w = 0 there, likely because the absorber cell has no usable formal solution (J_R = 0, so
    M1ThinFace's no-FS branch applies).
  - The face then takes the blend with the arithmetic-mean face tau, which is diffusive. It transmits about
    c E/(tau_cell), so the cavity fills up until E is about tau_cell x E_free.
- **Next step:** dump AL/HCL/DG at that face; compare with vet_sc's M1HrSf donor-scale path, which keeps such a
  face upwind.
- **Relevance:** it matters only for a fully enclosing cold opaque layer with an unresolved surface (tau_cell
  >> 1). It is not production-like: AG Car and BSG are unaffected. Runs: sp/runs/ss_*, ssd_*.

## 6. Open / pending
- pipe c220: done (sect. 4).
- Not pushed. Main has merged 8bfeae8c into bsgcand-1012.

## 7. Files
- Code: wt-jsfix 8bfeae8c (src/rad_m1/rad_m1_implicit.cpp M1ThinJs and its two call sites, the key read;
  src/rad_m1/rad_m1.hpp).
- Scripts:
  - build.sh, run1.sh, pool.sbatch, pool2.sbatch, cart.sbatch, spsuite.sbatch, lad.sbatch;
  - jobs_a/a2/b/b2/c/d/db/e.txt;
  - ana/lad.py (slab + shell), ana/t7j.py (T7 metrics, t7j.json), ana/evalJS.sh (Cartesian),
    ana/evalsp_js.sh (sp).
- Outputs: ana/pair_js.txt, ana/cart_js.txt; suiteF_1012/sp/ana/out_js/*; logs/.
