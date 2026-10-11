# NOTE sptests (sptests-1011): production-path (vet_col + vet_gd, spherical wedge) counterparts of the IDORT / Jiang tests

Branch `sptests-1011` (local, NOT pushed) in the worktree /raven/u/jinma/ATHENAK/wt-sptests, on deepatm-1011 3b50c192
(vetacc-1011 scheme + `vet_bc_diff = face`, default for vet_sc, opt-in for vet_gd). Work dir /raven/ptmp/jinma/sptests_1011.

## Verdict so far (Stage P, the priority subset; the rest of the suite is designed in Stage 0 and partly smoke-run)

- **Photosphere: fse+hc is right, production hr is wrong.** Against an exact spherical reference (vet_gd, poly r,
  production-like c dt/dx 1e3), hr is 2.2 % low in E at R_ph at N 256 (7.1 % at N 64), and it converges at first order.
  fse+hc is within 0.05-0.3 % and converges at order 1.5-2.4. The same holds for the thermal radiative-equilibrium
  version. AG Car's +4-6 % E at 1.000-1.006 R_ph under fse+hc is therefore the removal of hr's deficit.
- **Thin envelope:** the two arms are identical, about +1 % above J. The limit here is vet_gd's lateral angular
  resolution, not the switch.
- **Thermal channel beam on vet_gd:** both arms are stable at c dt/dx 30-3000 and over 1000 cycles. hr sends 41x the
  exact luminosity through the channel; fse+hc sends 2.3-3.4x. fse+hc costs Picard 30 at c dt/dx 3000 (hr: 11.5).
- New opt-in code: commit 82ef2e9c. With the new keys unnamed, output is bitwise identical to the base build.

---

## Stage 0. Design (written before any run)

### 0.1 Why a separate suite

- Accuracy evidence for the implicit VET so far is `closure = vet_sc` (Cartesian, uniform). Production (AG Car A, He giant, BSG)
  is `closure = vet_col` + `vet_gd = true` + `vet_gd_replace` on the 3-D sp wedge with a poly-stretched r grid.
- vet_sc refuses sp; vet_gd needs a 3-D wedge, ONE MeshBlock along x1, and a Marshak outer x1. So no problem runs on both.
- The two are tied here in two ways:
  1. **Large-radius control:** the same 1-D problem on a thin shell at r0 >> width (uniform r) has curvature terms O(dr/r0).
     It must converge to the Cartesian vet_sc answer as r0/width = 10, 100, 1000.
  2. **Genuinely spherical tests with exact (or high-accuracy 1-D ray) answers**, on uniform and poly r.

### 0.2 What the production path does at the inner boundary (and why code is needed)

- vet_gd sweeps the shells from `vlat_icut` (first shell where every column has tau_top >= `vet_col_lat_taucut` = 30) to the
  top. Below the cut, `vet_gd_replace` puts D = I/3. In a column that is thin everywhere, icut = 0: the sweep starts at the
  physical inner face.
- There the outward rays start from the **diffusion intensity** S + mu (-dS/dtau) (vet_bc_diff = cell: S of the first cell;
  face: 1.5 S_a - 0.5 S_b). In a thin first shell dS/dtau = dS/dr / chi is unbounded: this start is meaningless for an
  optically thin emitter, a point/sphere source or a beam injected at the inner boundary.
- **New opt-in key `vet_gd_inner = emit`** (default `diff` = bitwise): when icut = 0, the outward rays at the inner face
  carry the physical emission of the inner boundary:
  - `implicit_bc_x1min = flux`: a black-body sphere of flux F_in: isotropic outward I = 4 F_in / c (E units);
  - `implicit_bc_x1min = marshak`: the isotropic bath I = implicit_ebath_x1min (as vet_sc's vet_bc_bath);
  - `vet_gd_beam = true`: a collimated beam: the ONE set direction nearest (vet_gd_beam_th, vet_gd_beam_ph) carries
    I_b = e_b / w_d from the patch [th0, th1] x [ph0, ph1] of the inner sphere; nothing else enters. The M1 solve gets the
    matching lateral inner flux map F_in(theta, phi) = c e_b mu_d (patch area fraction of the cell), so moments and rays
    inject the same luminosity. The twin (ray-noise subtraction) gets the shell-mean inflow.
  - Below a cut (icut > 0) the start stays the diffusion intensity (the production thick-bottom case).
- pgen options (src/pgen/tests/rad_m1_tests2.cpp, `sph_shell` / `sph_atm`, read only when named):
  `core_r/core_rho/core_temp` (IDORT homogeneous sphere), `clump_*` (opaque clump), `gas_hold` (gas state re-imposed each
  stage: fixed emission), `atm_scale_h` (exponential atmosphere of IDORT 4.1.3 on a shell).

### 0.3 Arms (all with the explicit `vet_bc_diff` key)

| arm | keys on top of the test |
|---|---|
| base GD | closure vet_col, vet_gd, vet_gd_replace, taucut 30, GL nmu 10 x nphi 12 staggered, twin + twin_fuse + twin_det, vet_source_noesrc, vet_scatter ma/relaxed, vet_col_lat_offdiag operator, op_stencil, wall_interp, tangential false (agcar_A_xtf) |
| **P** (production) | GD + implicit_flux blend, faces all, beam halfrange, idort_f (tau0 1, alpha 1, flo 0.3, fhi 0.6), **xthin 30**, rotate_every 10, be |
| **F** (fse+hc) | P with implicit_blend_xthin 0, implicit_thin_switch fse, implicit_thin_corr hc |
| **c / f** suffix | vet_bc_diff = cell / face |
| **H** prefix | time_scheme hesdirk2, time2_vet_col predict, implicit_vimp, time2_vstage (vimp/vstage dropped where the gas is not coupled) |
| R0 suffix | rotate_every 0 (control for T ladders: vcpred found rotate 10 sets a dt-dependent ray-noise floor) |
| Cartesian reference | closure vet_sc, vet_tensor full, the same flux keys (hr <-> P, tEs-like fse+hc <-> F), same binary |

### 0.4 Design table: Cartesian reference test -> spherical-wedge counterpart

| # | Cartesian test (paper) | spherical-wedge counterpart (production path) | exact / reference solution | metric and pass criterion | map to the Cartesian result |
|---|---|---|---|---|---|
| T1 | thermal relaxation (JSD12 sec 4; IDORT 4.3.1 cases 1-2) | uniform static gas on the wedge r in [1,2], rho kappa_a = 100, c = 100, a = 1; (E,T) = (100,1) and (1,100)... scaled; reflect inner, Marshak top with the equilibrium bath; uniform and poly r | 0-D ODE dE/dt = c kappa (aT^4 - E), c_v dT/dt = -same (scipy, rtol 1e-12) | interior (r < 1.8) max rel. error of E, T at t = 10/(c kappa); T ladder dt/1..16: be order 1, H order 2 | same ODE as the IDORT/J21 figure; identical numbers expected |
| T2 | non-LTE scattering atmosphere (IDORT 4.1.3, J21 4.5) | thin shell r in [r0-10, r0+10], rho = 1e-3 exp(10 - (r - r0)), kappa_a = eps, kappa_s = 1 - eps, T = 1 held, Marshak inner bath 1, vacuum top; eps = 1e-2, 1e-4; r0/20 = 1, 10, 100, 1000 | IDORT analytic E = 1 - exp(-sqrt(3 eps) tau)/(1 + sqrt eps) (Eddington); high-accuracy 1-D reference: plane-parallel and spherical Lambda-matrix solve (ana/ref1d.py) | E/E_ref per tau band; max\|E_sp - E_cart\| decreasing ~10x per decade of r0 to the solver-difference floor | **large-radius control**: r0 -> inf must give the vet_sc (Cartesian, same N) answer |
| T3 | radiating pulse / pulse X-T ladders (thinsw/vcpred battery, rho kappa_s = 1) | sph_shell Gaussian pulse rho kappa_s = 1 between vacuum ends (Marshak bath 0 both ends; inner vet_gd_inner = emit with bath 0); r0/width = 1 (r in [1,3]), 10, 100, 1000 | self-convergence (X: N = 32..512; T: dt/1..16); r0 -> inf: the Cartesian vet_sc run of the same slab | X order 2, T order 1 (be) / 2 (H, rotate 0); \|E_sp - E_cart\| ~ 1/r0 | **large-radius control** + order |
| T4 | static diffusion (J21 4.4) / dynamic diffusion (IDORT 4.2.1) | Gaussian shell E = exp(-40 (r - r0)^2), rho kappa_s = 4e4, c = 10, r in [r0-1, r0+1], r0 = 2 (spherical) and r0 = 20, 200, 2000 | spherical: u = rE obeys 1-D diffusion: E = (1/r) int G(r,r') r' E0(r') dr', G = image-pair Gaussian (exact); plane limit (160Dt+1)^(-1/2) exp(...) | L1(E - E_exact)/L1(E) at sqrt(Dt) = 0.13, 0.18, 0.22; X order 2 (uniform, poly); T order be 1 / H 2 | Gaussian of J21 Fig. as r0 -> inf |
| T5 | Marshak / diffusion wave (JSD12 5.1; littests `marsh`) | spherical Marshak wave: cold gas (constant c_v), Marshak inner bath, r0/width = 1 and 1000 | r0 >> width: the planar reference of littests marsh (t6_marshak.py S_N); r0 ~ width: self-convergence | front position and E profile vs reference; X/T order | **large-radius control** |
| T6 | homogeneous sphere (IDORT 4.1.2, J21 4.6) / point source | (a) thin medium, black-body inner sphere r_in = 1, flux F_in (vet_gd_inner = emit), r in [1, 10]; (b) IDORT sphere: rho kappa_a = 10, T = 1 for r < 1, ~0 outside, r in [0.2, 6] | (a) J = (2 F_in/c)(1 - sqrt(1 - (1/r)^2)) (finite-sphere dilution), H = F_in/r^2, K/J exact; -> L/(4 pi r^2 c) at r >> 1. (b) I = B (1 - exp(-kappa l(r,mu))), quadrature (ana/sphere.py) | L1 and max rel. error of E, F, f_rr; X ladders uniform/poly; (a) emit vs diff | IDORT Fig 3 profiles (f = 1/3 inside, f -> 1 outside, F r^2 = const) |
| T7 | static shadow (IDORT 4.1.4, JSD12 6.5, littests shd) | black-body inner sphere r_in = 1 (emit) + opaque clump (tau ~ 1e3, cold, held) at r_c = 3, a = 0.5 on the wedge axis; 3-D wedge (32 x 32 lateral) | ray quadrature over the visible star disk minus the clump: umbra, penumbra (finite source) (ana/shadow_exact.py) | E/J_exact in umbra / penumbra / lit; penumbra width and umbra-edge position at r = 4, 5; Picard and stability | IDORT/J21 shadow geometry: sharp umbra + penumbra from the finite source |
| T8 | crossing beams (IDORT 4.1.1, J21 4.2) / DSJ12 beam | two beams injected at the inner sphere (vet_gd_beam, one set direction each, tilted +-4.4 deg in theta) crossing at r ~ 3; single radial beam | straight lines: E = e_b inside the tube {x : x - s n hits the patch}, F = c e_b n | beam E/e_b median, luminosity through shells (exact: constant), FWHM growth, centroid direction error, E beyond crossing | crossing beams pass through each other (moment closure with D from vet_gd) |
| T9 | linear radiation-acoustic waves (IDORT 4.3.2, JSD12 6.1) | rswm (vcpred): radiation-dominated moving shell, P_rad/P_gas ~ 10, rho kappa_s = 1000 and a thinner rho kappa_s = 10 variant; r in [1,3] | self-convergence X and T (dens, E, F1, v) | X order 2, T order H 2 / be 1 | the order criterion of the Cartesian wave test |
| T10 | extended grey atmosphere / Milne (thinfix/deepatm atm, Hopf) | sph_atm: rho = 100 (r/1)^-2 pure scattering, r in [1,5], flux bottom, Marshak top, steady; and r0 >> width (Hopf limit) | 1-D spherical Lambda-matrix solve S = J with r^2 H = const (ana/ref1d.py) | E/J_ref and self-conv. order per tau band (< 0.1, 0.1-1, 1-10, > 10): thin top >= 1st (target 2nd), thick 2nd | atm band orders of thinfix/deepatm (tJe+face: 2.0 everywhere) |

Grids: every test on uniform r and on a poly-stretched r (the `f_stretch_r_c*` poly of runs_5q, which refines toward r_in
by a factor ~2-3 like production). Radial tests: 4 x 4 lateral cells (theta [1.171, 1.971], phi [0, pi/3]); 3-D tests
(T7, T8): 24-32 cells per lateral direction.

---

## Stage P. PRIORITY SUBSET (coordinator 10-11): the photosphere question for AG Car, fse+hc vs production hr

Asked: which of `hr` (production, implicit_blend_xthin 30) and `fse+hc` is right at the photosphere / thin envelope, with a
resolution study on a stretched grid like production, and a thermal beam/emitter proof on vet_gd. Code 82ef2e9c, binary
bin/athena_cpu_82ef2e9c (md5 df208420...). Arms P = production hr (`all`/halfrange/idort_f, xthin 30, be, rotate 10,
GL 10x12, twin+fuse+det, scatter ma, offdiag operator, `vet_bc_diff = face` named), F = P with xthin 0, fse, hc.

### P1. Spherical grey atmosphere through the photosphere (case `patm`)

- rho kappa = 896 r^-8 on r in [1, 5] (4 x 4 / 8 x 8 / 16 x 16 / 32 x 32 / 64 x 64 lateral cells on the sp wedge),
  pure scattering S = J (the steady state equals grey radiative equilibrium), flux 1 in at r = 1, Marshak vacuum top.
  tau(base) = 128, tau = 1 at r = 2.0, **R_ph(tau = 2/3) = 2.1185**, thin envelope to 2.36 R_ph (AG Car runs to 3 R_ph).
  vet_col below tau = 30 (r < 1.23), vet_gd above (replace). Poly-stretched r (runs_5q poly: cells 2.4x finer at the top),
  N = 64 ... 1024.
- Steady state in two stages: A = 80 steps at implicit_cfl 1e6 (deep diffusion), restart, B = 300 steps at a
  production-like **c dt/dx_min = 1e3** (hr's xthin weight is dt-dependent; AG Car has c dt/dx ~ 750 at R_ph). Flux is
  conserved to max|F r^2 - 1| = 3.8e-4 (N 256), 3e-5 (N 1024). B vs a 600-step B: E changes <= 5e-3 (thin top).
- **Reference** (ana/ref1d.py `sphere_vef`, ref/patm_ref.npz): 1-D spherical p-ray formal solution (linear-in-tau source,
  exact exponentials, limb ray pair) inside a variable-Eddington-factor iteration with r^2 H = 1 exactly; inner diffusion
  intensity at tau = 20. 1160 vs 2320 nodes differ by 1.2e-4 (so ~3e-5 error). f(top) = 0.875, g = H/J(top) = 0.933.
  (A direct Lambda-matrix solve of S = J lost flux by 50 % in the thick part; not used.)

**Result (16 x 16 lateral cells, poly r; E/J_ref - 1 at R_ph, L1 per tau band; fig/patm_L16_PvsF.png):**

| N | P: E/J-1 at R_ph | P: dR_ph/R_ph | P: band 0.1-1 / 1-10 | F: E/J-1 at R_ph | F: dR_ph/R_ph | F: band 0.1-1 / 1-10 |
|---|---|---|---|---|---|---|
| 64 | -7.1 % | -1.1e-2 | 6.5e-2 / 2.2e-2 | -0.73 % | -2.7e-4 | 7.8e-3 / 7.4e-3 |
| 128 | -3.9 % | -6.6e-3 | 3.4e-2 / 1.4e-2 | -0.26 % | -2.3e-4 | 1.3e-3 / 1.8e-3 |
| 256 | -2.2 % | -3.8e-3 | 1.5e-2 / 1.1e-2 | -0.05 % | -4.9e-5 | 2.7e-3 / 8.6e-4 |
| 512 | -1.1 % | -2.0e-3 | 5.0e-3 / 9.4e-3 | +0.04 % | +8.5e-5 | 4.8e-3 / 7.0e-4 |
| order | **~0.9 (first order)** | | 0.9-1.1 / 0.4-0.6 | **1.5-2.4** | | 2.1-2.6 (then the floor) |

- **fse+hc is right at the photosphere; production hr is low there.** hr under-predicts E over 0.8-1.2 R_ph (a 2 % dip
  at N 256, 7 % at N 64, centred on tau ~ 0.5-1) and converges only at first order; fse+hc stays within 0.1-0.3 %. The
  photosphere radius (E = J_ref(2/3)) sits 0.4-1.1 % of R_ph too deep with hr, 0.005-0.03 % with fse+hc.
- This is the same sign as AG Car (fse+hc raises E by 4-6 % at 1.000-1.006 R_ph): there fse+hc removes the hr deficit.
  The AG Car magnitude (4-6 %) is that of hr's deficit at a coarse effective tau-resolution (here 4-7 % at N 64-128).
- Both arms are identical above ~1.2 R_ph: E is 1.0-1.4 % above J_ref in the thin envelope. This plateau is set by the
  vet_gd formal solution's LATERAL angular resolution, not by the switch: tau < 0.1 band L1 at N 256 = 2.1 % (8 x 8),
  0.9-1.2 % (16 x 16), 10 % (4 x 4; and it grew with N there, +29 % at N 1024 - fixed Cartesian ray directions
  interpolated bilinearly across 0.2-rad cells once per shell). AG Car has 0.008-rad cells; see P1b.
- Cost: Picard 2.8-2.9 per step for both arms (be, c dt/dx 1e3); stable at implicit_cfl 1e6 too.

- **Thermal version (patmT, absorption + gas in radiative equilibrium, rescaled, 16 x 16):** identical picture.
  P: E/J-1 at R_ph -7.1 % (N 64), -2.2 % (N 256); F: -0.80 %, +0.004 %. Picard 3.0 both. So the hr deficit is not a
  scattering-source artefact.
- Default bitwise: bin/athena_cpu_82ef2e9c vs base 3b50c192 on patm (P arm, vet_gd) and rsw (vet_gd + gas): every tab
  file identical (verify/result.txt).

### P2. Thermal beam / emitter on vet_gd (case `pinc`)

- Hot opaque core (r < 1.5, rho kappa_a 50, B = 1, held) under a cold opaque absorbing shell 2 < r < 4 (rho kappa_a 500,
  tau = 1000) pierced by a radial cone channel of half-angle 0.06 rad on the wedge axis; thin elsewhere; reflect +
  diffusion start at r = 1. 3-D wedge 64 x 48 x 48 (theta +-0.35 periodic, phi 0..pi/3), vet_gd covers everything
  (no column reaches tau 30 from the top through the channel). The analogue of the Cartesian channel beams xb20/ba0 on
  which hc diverged (thinsw battery).
- Exact J and H_r beyond the channel by ray quadrature (ana/pin_exact.py, ref/pinc_exact_64.npz).

| run (100 cycles unless noted) | Picard mean/max | NONCONV | L(r=4.6)/L_exact | beam E/J_ex median r 4.3 / 4.6 / 4.9 | spurious E at channel-wall cells |
|---|---|---|---|---|---|
| P, c dt/dx 300 | 5.0 / 7 | 0 | **41** | 14 / 9.9 / 8.0 | 0.36 (axis 0.05) |
| F, c dt/dx 300 | 3.6 / 5 | 0 | **2.3** | 0.90 / 0.72 / 0.47 | 0.064 (axis 0.048) |
| P, c dt/dx 3000 | 11.5 / 65 | 0 | 41 | 14 / 9.9 / 8.0 | |
| F, c dt/dx 3000 | 30.2 / 78 | 0 | 2.3 | 0.90 / 0.72 / 0.47 | |
| F, c dt/dx 30 | 3.5 / 5 | 0 | 2.3 | 0.91 / 0.72 / 0.49 | |
| F, 300 cycles | 3.7 / 5 | 0 | **3.5** | 2.3 / 2.0 / 1.3 | |

- **Both are stable** (no divergence, no NONCONVERGED, cfl 30-3000) - unlike the Cartesian vet_sc thermal beams with hc.
- **hr is badly wrong on a thermal channel beam**: 40x the exact luminosity leaves the channel, E piles up in the
  channel cells next to the opaque walls (7x the axis value). This is the Cartesian hr "beam E/J 2.9-7.1" in 3-D.
- fse+hc is ~20x closer, but the channel luminosity still grows between 100 and 300 cycles (2.3 -> 3.5); the 1000-cycle
  history shows it is bounded, not growing: L/L_ex = 2.31, 3.42, 2.60 at cycles 100, 200, 300 (hr: 41.1, 41.2, 40.9).
  The +-25 % wander is the rotate_every 10 ray pattern on a 4-cell-radius channel.
- Cost: at c dt/dx 3000 fse+hc needs Picard 30 (P 11.5): fse+hc is expensive on stiff thermal beams at large dt.

### P3. Thermal radiative-equilibrium photosphere (case `patmT`): failures explained

- All patmT runs of the first batches stopped at cycle 2 with "dt = -3.0e-01 ... non-finite or non-positive time step".
  The line before it is "non-finite hydro state ... top active cell, rho 2.6e-6, T = -nan". The negative dt is only the
  FATAL path after a NaN gas state; the cause is my input: a = 1e64 with T ~ 1e-16 (to keep the hydro dt out of the
  way) gave the top cells a gas heat capacity ~1e-20 of their radiation energy, and the gas-radiation exchange Newton
  returned T = NaN at c dt/dx >= 1e3 (earlier runs survived only because rad_signal_speed limited dt to 5e-4).
  Not a production regime. Fix: rho x 1e20, kappa x 1e-20 (rho kappa unchanged), rerun: all patmT runs now complete (Picard 3.0, no NaN); results in P1. Code note: in that
  extreme limit the exchange solve produces NaN instead of clamping / falling back.

## Status of the rest of the suite (at this hand-back)

- Done: Stage 0 design; code (82ef2e9c: vet_gd_inner = emit, vet_gd_beam, pgen core / clump / screen / gas_hold /
  atm_scale_h); harness (ana/spl.py, spt.py, ref1d.py, patm_eval.py, pin_exact.py, pin_eval.py, pin_hist.py,
  fig_patm.py); smoke runs of thin, core, nlte2, pulse, relax, shd and beam.
- Smoke findings:
  - thin sphere (T6a): `emit` is needed. With the old diffusion start, E/E_exact = 0.79 in the first cell. With emit it
    is 0.88 (P) and 0.99 (F).
  - The first held-gas inputs did not emit (gas heat capacity << exchange). They are now rescaled (kscale).
  - rad_signal_speed is now off for static gas: it had capped dt at 5e-4.
- Still running when this was written:
  - patm with production-like lateral cells (Lq: 64 x 64 at 0.008 rad), N 64-256;
  - L32 lateral at N 256;
  - L16 at N 1024.
- Not started: T1-T10 ladders, large-radius control, cell vs face, hesdirk2 arms, shadow / beam exact, ports.
