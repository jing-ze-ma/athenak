# NOTE sptests (sptests-1011): production-path (vet_col + vet_gd, spherical wedge) counterparts of the IDORT / Jiang tests

Branch `sptests-1011` (local, NOT pushed) in the worktree /raven/u/jinma/ATHENAK/wt-sptests, on deepatm-1011 3b50c192
(vetacc-1011 scheme + `vet_bc_diff = face`, default for vet_sc, opt-in for vet_gd). Work dir /raven/ptmp/jinma/sptests_1011.

## Verdict so far (Stage P, the priority subset; the rest of the suite is designed in Stage 0 and partly smoke-run)

- **Photosphere: fse+hc is right, production hr is wrong.** Against an exact spherical reference (vet_gd, poly r,
  production-like c dt/dx 1e3), hr is 2.2 % low in E at R_ph at N 256 (7.1 % at N 64), and it converges at first order.
  fse+hc is within 0.05-0.3 % and converges at order 1.5-2.4. The same holds for the thermal radiative-equilibrium
  version. AG Car's +4-6 % E at 1.000-1.006 R_ph under fse+hc is therefore the removal of hr's deficit.
- **Thin envelope:** the two arms are identical, about +1 % above J at 16 x 16 lateral cells. This is vet_gd's lateral
  angular resolution, not the switch. It falls about linearly with the lateral cell angle: +0.3 % at the production 0.008 rad
  (Stage P1b).
- **Large-radius control (Stage L):** vet_gd on a thin shell converges to the Cartesian vet_sc answer and to the exact
  answer as 1/r0. Pulse: 1e-3 floor. Non-LTE atmosphere: 0.25 %. This requires `vet_col_source = gas` in that
  held-temperature scattering test. The default `relaxed` source is -5 % (c dt/dx 1e3) to +13 % (cfl 1e6) there.
- **P2 diagnosis (Stage P2b):** the channel over-luminosity of fse+hc (2.2-2.5x) comes from the moment step, not from the
  formal solution: L_fs = 0.3 -> 0.8 of exact with 120 -> 480 directions. The cbfix sf mechanism is absent; no port.
- Final code: prodcand e66ee669 reproduces the Stage-P fse+hc results bitwise (with ic_rad = none).
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

## Binaries from here on (coordinator 10-11)

- **prodcand e66ee669** = bin/athena_cpu_e66ee669 (copy of /raven/ptmp/jinma/prodcand_1011/bin, md5 39ab667b48f4e15d29bf13e39f8179ce):
  prodcand-1011 with sptests 82ef2e9c + 9b8552ff merged, M1_VGD_NMAX kept at 432. Run tree `runs_pc/` (env SPT_RUNS).
  - prodcand turns on `ic_rad = fs_relax` by default on fresh vet_col + vet_gd starts with a fluid. It relaxes the
    initial radiation field over 40 steps before cycle 0. The time-dependent pulse then lost 62 % of its energy before
    cycle 0. All test inputs now set `ic_rad = none` (spt.rad_base). With it, prodcand reproduces 82ef2e9c bitwise
    (pulse, cart and sp).
- **sptests 9b8552ff** = bin/athena_cpu_gdhr2 (CPU only): adds the `m1_gdhr` output (the vet_gd half-range ratios) and
  M1_VGD_NMAX 2048. It is used only for the 2x/4x GL direction study of P2. Default output is bitwise identical to 82ef2e9c
  (nlte2 cmp).
- Scope change (user 10-11): production hr is no longer an arm. Only fse+hc arms are run; the hr numbers above are cited.

## Stage P1b. Lateral resolution of the thin-envelope error (patm, fse+hc)

Thin envelope = r > 1.2 R_ph. `env` is the volume-weighted SIGNED mean of E/J_ref - 1 there; the 'tau<0.1' column is L1.
Binary 82ef2e9c, poly r, c dt/dx_min 1e3.

| lateral cells (angular size) | N 64 env | N 128 env | N 256 env | N 512 env | N 1024 env | E/J-1 at R_ph, N 256 |
|---|---|---|---|---|---|---|
| 4 x 4 (0.2 rad) | +3.1 % | +5.5 % | +9.9 % | +16 % | +22 % | -0.19 % |
| 8 x 8 (0.1) | | | +2.3 % | | | -0.02 % |
| 16 x 16 (0.05) | +0.42 % | +0.87 % | +1.14 % | +1.34 % | +1.59 % | -0.05 % |
| 32 x 32 (0.025) | | | +0.64 % | | | -0.09 % |
| Lq 64 x 64 (0.008, production-like) | -0.03 % | +0.24 % | +0.34 % | | | -0.05 % |

- **The +1 % thin-envelope "floor" is not a floor.** It falls about linearly with the lateral cell angle: +9.9, 2.3, 1.14,
  0.64 % from 4x4 to 32x32 at N 256. At the production angular cell (0.008 rad) it is +0.24 % (N 128) and +0.34 % (N 256).
- At fixed lateral cells it grows slowly with N_r: 4x4 +3 % -> +22 %; 16x16 +0.4 % -> +1.6 % over N 64-1024.
  The likely cause is the per-shell bilinear lateral interpolation of the global-direction intensities. Each shell
  crossing adds an error of about (lateral cell angle) x (dI/dmu), so it accumulates with the number of shells.
- The photosphere (R_ph) values do not depend on the lateral grid (Lq = L16 to 1e-4). The P1 verdict on R_ph stands.
  fse+hc reaches +0.10 % at N 1024 (16x16), so its R_ph error has hit its floor. hr is at -0.33 %, still first order.
- Lq at N 64-128 has a smaller dt than L16, because its lateral cell sets dx_min. fse+hc is insensitive to this; the hr
  rows of the Lq runs are not dt-matched.

## Stage L. Large-radius control: vet_gd on a thin shell -> the Cartesian vet_sc answer

Same 1-D problem on a shell of width w at r0/w = 1, 10, 100, 1000 (uniform r, 4 x 4 lateral cells), against the same
binary's Cartesian vet_sc slab (CFf: vet_sc + the same fse+hc keys). Data: runs_pc (prodcand e66ee669) and runs (82ef2e9c,
identical where both ran). ana/lr_eval.py.

**T3 pulse** (Gaussian, rho kappa_s = 1, vacuum ends, t = 1, c dt/dx 0.5 at N 64; L1 of (E_sp - E_cart)/max E_cart):

| r0/w | N 64 | N 128 | N 256 | E_peak sp/cart (N 256) |
|---|---|---|---|---|
| 1 | 8.3e-2 | 7.9e-2 | 7.4e-2 | 1.166 (true spherical focusing) |
| 10 | 9.2e-3 | 8.4e-3 | 7.9e-3 | 1.011 |
| 100 | 1.8e-3 | 1.7e-3 | 1.7e-3 | 1.0005 |
| 1000 | 1.1e-3 | 1.1e-3 | 1.1e-3 | 0.9995 |

- The difference falls as 1/r0 down to a 1.1e-3 floor. The floor is independent of N, so it is the difference between
  the two formal solvers (vet_gd GL 10x12 product set vs vet_sc Gauss 8x32), not curvature or resolution.
  The hr arm (82ef2e9c) gives the same numbers to 2 digits.

**T2 non-LTE atmosphere** (IDORT 4.1.3: chi = 1e-3 e^(20-z), eps = 1e-2, held gas B = 1, bath 1 below, vacuum top).
Exact reference: plane-parallel Lambda-matrix solve (ana/ref_nlte.py -> ref/nlte2_ref.npz, 700 vs 1400 nodes 3e-3).
Two-stage steady state (`_tb`: 80 steps at cfl 1e6, then 300 at c dt/dx 1e3). Mean E/J_exact - 1 per tau band:

| run (N 320) | tau < 0.1 | 0.1-1 | 1-10 |
|---|---|---|---|
| Cartesian vet_sc CFf | -0.02 % | +0.05 % | +0.09 % |
| sp r0/w 1, FfSg | -40 % | -23 % | -8.4 % (true sphere: dilution) |
| sp r0/w 10, FfSg | -7.1 % | -3.8 % | -1.1 % |
| sp r0/w 100, FfSg | -1.05 % | -0.42 % | -0.04 % |
| sp r0/w 1000, FfSg | -0.26 % | -0.03 % | +0.07 % |
| sp r0/w 1000, FfSg, N 640 | -0.24 % | -0.02 % | +0.07 % |
| eps = 1e-4 (nlte4), r0/w 1000 vs cart, FfSg | -0.25 % | -0.10 % | -0.04 % |
| sp r0/w 1000, Ff (default vet_col_source = relaxed), c dt/dx 1e3 | **-5.4 %** | **-2.7 %** | +0.4 % |
| same at implicit_cfl 1e6 | **+13 %** | **+10 %** | -0.02 % |
| hr (82ef2e9c, relaxed), c dt/dx 1e3 / cfl 1e6 | -5.3 % / +8.5 % | -7.3 % / +23 % | -0.8 % / +7.5 % |

- **vet_gd converges to vet_sc and to the exact answer as r0 -> inf** (FfSg: 1/r0 down to 0.25 % at the thin top).
  This ties vet_gd to the Cartesian evidence.
- **But only with `vet_col_source = gas`** (arm suffix Sg). With the production default `relaxed`, vet_gd's formal-solution
  source in a scattering-dominated layer is wrong for held gas:
  - The relaxed source takes the scattering part from the 0-D end-of-step E* = (E^n + c dt kappa_a B)/(1 + c dt kappa_a).
    At large c dt kappa_a, E* -> B.
  - The dump (m1_fs) at implicit_cfl 1e6 shows S_fs = 0.999 where the exact S = eps B + (1 - eps) J = 0.36.
    J_fs = 5.4 J_exact at the top; the M1 E stays within 13 %, because only the closure shape is wrong.
  - This test holds T (B != J). In production the gas is coupled and near radiative equilibrium, so E* ~ J there. The
    effect is limited to non-equilibrium scattering layers at large c dt kappa_a. patmT (coupled gas, absorption) shows
    nothing. This is not proven for an electron-scattering envelope with coupled gas; see "open" below.
- Lateral cells do not matter in the plane limit: 16 x 16 equals 4 x 4 to 3e-3 at r0/w 1000.
- The steady state is reached: 240 vs 80 steps changes E by < 3e-3. rotate_every 0 vs 10 makes no difference.
- The deep base cells (first 2 of 320, tau ~ 5e5) carry a start-up artefact (E = 7e4, 5e-3). It is identical in every
  geometry and arm, and excluded from the bands.

## Stage P2b. Thermal channel beam diagnosis (why fse+hc is 2.3x the exact luminosity)

ana/pin_eval.py (now also reports the formal solution's own luminosity L_fs and the face luminosity) and ana/pin_hr.py (the
M1 budget of the channel and the vet_gd half-range ratios). 100 cycles, c dt/dx 300, N_r 64.

| run | lateral cells | vet_gd directions | rotate | L_M1/L_ex (r 4.3 / 4.9) | L_fs/L_ex | J_fs/J_ex in beam |
|---|---|---|---|---|---|---|
| F | 48 x 48 | 120 (GL 10x12) | 10 | 2.25 / 2.35 | 0.29 | 0.06 |
| F R0 | 48 x 48 | 120 | 0 | 0.88 / 0.92 | 0.53-0.58 | 0.04 |
| F G2 | 48 x 48 | 480 (20x24) | 10 | 2.22 / 2.30 | 0.79-0.82 | 0.59 |
| F G2 R0 | 48 x 48 | 480 | 0 | 1.94 / 2.01 | 0.80-0.84 | 0.57 |
| F | 96 x 96 | 120 | 10 | 2.44 / 2.56 | 0.13 | |
| F R0 | 96 x 96 | 120 | 0 | 1.23 / 1.28 | 0.38-0.41 | |
| hr (P), cited | 48 x 48 / G2 R0 / 96 x 96 | | | 41 / 41 / 20 | same FS | |

Exact references: smooth cone (ref/pinc_exact_64.npz) and the code's staircase hole (PIN_STAIR, ref/pincs_exact_64.npz).
The staircase changes L_ex by -11 % to +1 % (open area 1.01 x the cone). So the reference geometry is not the cause.

- **(a) The formal solution is not the cause.** Its own luminosity is LOW, 0.29 x exact with 120 directions, because the
  0.011-sr cone is undersampled. It rises to 0.8 with 480 directions. The M1 luminosity stays at 2-2.3x throughout.
  The moment step therefore adds the excess, about 2.5-8x the formal solution's own transmission.
- **(b) Angular resolution** fixes the formal solution (0.29 -> 0.8), not the M1 answer (2.25 -> 2.2; R0 1.9-2.0).
  The 120-direction R0 run at 0.9 is a coincidence: the fixed set happens to lie off the channel axis.
  4x GL (1920 directions, rotate 10 and 0): still running (job p2d) at this hand-back.
- **(c) Lateral resolution 2x** does not help fse+hc (2.25 -> 2.5). The formal solution gets worse (0.29 -> 0.13; more
  shell interpolations of a narrow beam). hr halves (41 -> 20).
- **(d) The walls are held:** gas_hold re-imposes rho, T, v of the whole gas every stage, with gas_feedback false. The
  shell is a fixed cold black absorber, as the exact solution assumes.
- **Budget of the M1 channel** (pin_hr, F):
  - 6.43e-3 enters at r = 2;
  - 94.2 % leaks laterally into the walls;
  - 5.8 % leaves at r = 4 (2.3x the exact outflow).
  - In the channel E/J_fs = 3-10, with the formal solution's ratios hp_r = 0.98 and |h_theta| <= 0.1 (beam-like).
  - The gap below the inlet has hp_r = 0.81.
  - So energy entering the channel at broad angles is given the narrow beam shape of the channel's own formal solution,
    and too little of it is sent into the walls.
- **cbfix mechanism (implicit_hr_face = sf) on vet_gd: not the cause here.**
  - sf triggers where a cell's net FS face flux exceeds c J (Jh/J > 1).
  - In the channel and on the axis Jh/J = 1.000 everywhere: the FS flux there is consistent with its own J.
  - It triggers only in the opaque wall cells (J_fs ~ 1e-11, 37 % of wall-ring cells). Their faces have dtau ~ 30 per
    cell, so the AP weight makes them central faces and the half-range ratios are not used.
  - The over-luminosity is the opposite mismatch: M1 E >> J_fs in cells whose FS shape is right. The sf donor scale does
    not see it.
  - **No sf port to vet_gd** (it would not change this case).

### Stage-P verdict on the final code (prodcand e66ee669, runs_pc, ic_rad = none)

- patm L16 fse+hc at N 64-512, patmT L16 fse+hc at N 64-256, and pinc fse+hc (rotate 10 and 0) are **bitwise identical** to
  82ef2e9c (tab/bin cmp). E/J-1 at R_ph = -0.73 / -0.26 / -0.05 % (patm N 64/128/256) and -0.80 / -0.24 / +0.004 % (patmT).
  pinc L/L_ex = 2.25 (rotate 10) and 0.88 (rotate 0). The Stage-P verdict holds on the final code.

## Stage C. Does production's vet_col_source = relaxed bias a coupled-gas atmosphere? (coordinator 10-11 item 1)

Case `patmT_L16_e<k>_s` is patmT (rho kappa = 896 r^-8, flux 1 in, vacuum top, coupled gas, poly r, 16 x 16 lateral
cells, two-stage steady state ending at c dt/dx 1e3). The extinction is split into absorption eps chi and scattering
(1 - eps) chi, with eps = 1e-2 and 1e-4. In the steady state the gas is in radiative equilibrium (B = J), so S = J. The
exact reference is the same as for patm. Binary xs1 (prodcand e66ee669 + the opt-in limiter, default off; data-bitwise to
prodcand), run tree runs_xs.

| arm | eps | N | E/J-1 at R_ph | dR_ph/R_ph | band 0.1-1 L1 | env | Picard |
|---|---|---|---|---|---|---|---|
| Ff (be, relaxed) | 1e-2 | 64 / 256 | -0.724 / -0.052 % | -2.7e-4 / -4.9e-5 | 7.8e-3 / 2.7e-3 | +1.14 % (256) | 2.98 / 2.87 |
| FfSg (be, gas) | 1e-2 | 64 / 256 | -0.724 / -0.052 % | -2.7e-4 / -4.9e-5 | 7.8e-3 / 2.7e-3 | +1.14 % | 2.97 / 2.95 |
| HFf (hesdirk2 predict, relaxed) | 1e-2 | 64 / 256 | -0.729 / -0.035 % | -2.7e-4 / -2.0e-5 | 7.8e-3 / 2.7e-3 | +1.14 % | 2.89 / 2.94 |
| HFfSt (hesdirk2, time2_vet_col_source state) | 1e-2 | 64 / 256 | -0.730 / -0.035 % | -2.8e-4 / -2.0e-5 | | +1.14 % | 2.98 / 3.00 |
| Ff / FfSg / HFf / HFfSt | 1e-4 | 256 | -0.051 / -0.051 / -0.035 / -0.035 % | | | +1.14 % | 2.80 / 2.81 / 2.99 / 3.00 |

- **No bias.** With coupled gas in radiative equilibrium, relaxed and gas agree to <= 1e-5 in E at R_ph, for both eps.
  So does hesdirk2's state source. Picard is unchanged (2.8-3.0).
- The relaxed source is biased only where the gas is held away from radiative equilibrium (B != J), as in the IDORT
  non-LTE test of Stage L. There E* -> B at large c dt kappa_a.
- In production such layers are dynamic: shocks, adiabatically cooled wind, a thermal time longer than dt. The size there
  is set by how far B is from J, not by the photosphere.
- **Recommendation: no change to production.** If a non-equilibrium scattering layer matters, `vet_col_source = gas` is
  exact there (Stage L) and costs no Picard (2.0 vs 2.0 in nlte2); `state` is not needed.
- vet_scatter_j = current (the cheap variant) is refused by the code with vet_scatter_form = ma (non-convergence on the BSG).

## Stage X. P2 consistency limiter `implicit_hr_excess = iso` (opt-in, vet_gd only; coordinator item 2)

Code (branch sptests-1011):
- 8af6b2ae: merge of prodcand e66ee669 into sptests-1011.
- f94ce268: the limiter.
- b6f91052: its steady-regime gate.

Binary bin/athena_cpu_xs2 (md5 4ac3fba2f2bf3e2915c907a95e08b1cd), run tree runs_xs2. With the key unnamed, the data are
bitwise those of prodcand (pinc bin arrays, pulse tab; the bin header differs only by the added default key in the
parameter dump).

- **Rule (M1HrExcess).** In a cell where the lagged E^n exceeds the vet_gd formal-solution J (vgd_hr HRJ), only the J
  part of E keeps the formal solution's half-range shape; the excess gets the isotropic +-1/4:
  `q = J/E^n, h+ -> q h+ + (1 - q)/4, h- -> q h- - (1 - q)/4`.
  - It is applied on the x1 and lateral faces, and only where c dt/dx1 >= implicit_hr_excess_cfl (default 10).
  - HCL >= 0 and HCR <= 0 are kept; the hc correction is untouched.
- **Why the gate.** Ungated, the time-resolved pulse (c dt/dx 0.5) changed by 9 %: in a light front E^n differs from the
  steady formal solution's J legitimately. Gated, the pulse is bitwise unchanged.

| pinc (fse+hc, 100 cycles, c dt/dx 300) | directions | rotate | L_M1/L_ex r 4.3 / 4.9 | L_fs/L_ex | beam E/J_ex median r 4.3 | Picard mean/max |
|---|---|---|---|---|---|---|
| Ff | 120 | 10 | 2.25 / 2.35 | 0.29 | 0.90 | 3.6 / 5 |
| **FfX** | 120 | 10 | **0.31 / 0.32** | 0.29 | 0.061 (J_fs 0.062) | 7.0 / **198** |
| Ff R0 | 120 | 0 | 0.88 / 0.92 | 0.53 | 0.19 | 2.2 / 4 |
| FfX R0 | 120 | 0 | 0.32 / 0.34 | 0.53 | 0.033 | 3.2 / 4 |
| FfQ | 432 (18x24) | 10 | 2.29 / 2.38 | 0.76 | 1.88 | 3.7 / 5 |
| **FfQX** | 432 | 10 | **0.69 / 0.71** | 0.76-0.80 | 0.55 (J_fs 0.54) | 4.3 / 8 |
| FfQ R0 | 432 | 0 | 2.02 / 2.10 | 0.78 | 1.55 | 2.2 / 4 |
| FfQX R0 | 432 | 0 | 0.63 / 0.65 | 0.78-0.82 | 0.53 | 3.9 / 5 |

- **The limiter makes the moment solution consistent with the formal solution.** The channel luminosity goes from a
  non-converging 2.0-2.4x to 0.31 (120 directions) and 0.69 (432 directions). That is within 10-40 % below the formal
  solution's own L_fs, which converges to exact with the direction count (0.29 -> 0.76; 0.8 at 480).
- The remaining error is an under-luminosity set by angular resolution, and it converges.
- **Atmosphere (P1) is barely touched** (16 x 16, poly r):
  - patm N 64: E/J-1 at R_ph -0.58 % (off -0.73 %), env +0.54 % (off +0.42 %).
  - patm N 256: R_ph -0.044 % (off -0.051 %), env +1.14 % (unchanged).
  - patmT N 256: R_ph +0.011 % (off +0.004 %).
  - Picard 2.8-3.0 in all of them.
- **Cartesian beams (vet_sc) are untouched by construction** (vet_gd src 2 only). CFfX equals CFf bitwise.
- **Open:**
  - The 120-direction FfX run has one step at Picard 198 (mean 7.0); the 432-direction runs do not (max 8).
  - Not tested on AG Car; its cost and effect there should be measured before any default change.

## Stage 3a. Exact-solution tests, fse+hc (binary xs2 = prodcand + opt-in, run tree runs_xs2, ic_rad = none)

ana/s3_eval.py (modes thin / dif / relax / core / shd).

- **T6a black-body sphere in a thin medium** (r_in = 1, flux 1, r in [1, 10], 4 x 4 lateral cells, vet_gd_inner = emit).
  Exact E = 2(1 - sqrt(1 - 1/r^2)), F r^2 = 1.

  | grid | N 64 | N 128 | N 256 |
  |---|---|---|---|
  | uniform r: L1 E/E_ex-1 | 1.2 % | 2.4 % | 3.2 % |
  | uniform r: max\|F r^2 - 1\| | 1.9e-2 | 5.0e-3 | 1.3e-3 |
  | poly r: L1 E/E_ex-1 | 1.2 % | 2.0 % | 2.8 % |
  | poly r: max\|F r^2 - 1\| | 2.7e-2 | 6.8e-3 | 1.7e-3 |

  - The flux converges at second order.
  - E grows with N: this is the 4 x 4 lateral-interpolation error of P1b (more shells, more interpolations).
    16 x 16: queued (s3b).
  - vet_bc_diff cell = face bitwise here, because the emit start does not use the diffusion intensity.
- **T4 spherical diffusion** (J21 static diffusion on the sphere, r0 = 3, rho kappa_s = 4e4, c = 10; exact = image-free
  Gaussian convolution of u = rE). L1 relative error at sqrt(Dt) = 0.22, dt = T/512:

  | arm | N 64 | N 128 | N 256 | X order |
  |---|---|---|---|---|
  | be, uniform | 4.3e-3 | 1.5e-3 | 7.8e-4 | 1.6, 0.9 (the dt error floor) |
  | be, poly r | 9.5e-3 | 2.8e-3 | 1.1e-3 | 1.8, 1.3 |
  | hesdirk2 predict, uniform | 3.8e-3 | 9.4e-4 | 2.4e-4 | **2.0, 2.0** |

  - T ladder at N 128 (self-convergence): be order 1.2-1.6; hesdirk2 order 2.3 (1e-6 level).
- **T1 thermal relaxation** (rho kappa_a = 100, c = 100: 1/(c kappa) = 1e-4).
  - Outputs at 2.5-10 relaxation times: both schemes are at 1e-9 by t = 5e-4. relax2 (T0 = 100) at dt0/1: be 4e-5 in E,
    hesdirk2 5e-4 (the stiff first step of hesdirk2).
  - The early-time ladder (0.5-2 relaxation times, `relax*_early_T`) is queued (s3b).
- **T9 radiation-acoustic shell** (rsw, P_rad/P_gas ~ 10, rho kappa_s = 1000):
  - Interior (r < 2.7) X self-convergence order 2.5-2.8 (be) and 1.8-2.0 (hesdirk2) in E, rho, v.
  - T order: be 1.0; hesdirk2 1.4-1.8.
  - The outermost cell pair at the Marshak top does NOT converge (rho 1.0001 vs 1.013 at N 128 / 256). This is a top
    boundary layer of the hydro with radiation force; it is excluded from the order.
- T6b IDORT sphere (exact by quadrature, s3_eval core), T7 shadow (exact visible-disk quadrature, s3_eval shd) and T8 beams:
  queued (s3b).
  - Smoke (hr, 32 radial cells): shadow lit-region L1 42 %. This is ray effects of 120 directions against a star of
    0.2 sr; the queued runs include the 432-direction set.

## Status of the rest of the suite (at this hand-back)

- Done: Stage 0 design; code (82ef2e9c: vet_gd_inner = emit, vet_gd_beam, pgen core / clump / screen / gas_hold /
  atm_scale_h); harness (ana/spl.py, spt.py, ref1d.py, patm_eval.py, pin_exact.py, pin_eval.py, pin_hist.py,
  fig_patm.py); smoke runs of thin, core, nlte2, pulse, relax, shd and beam.
- Smoke findings:
  - thin sphere (T6a): `emit` is needed. With the old diffusion start, E/E_exact = 0.79 in the first cell. With emit it
    is 0.88 (P) and 0.99 (F).
  - The first held-gas inputs did not emit (gas heat capacity << exchange). They are now rescaled (kscale).
  - rad_signal_speed is now off for static gas: it had capped dt at 5e-4.
- Done since: P1b lateral study, Stage L large-radius control (pulse, nlte2, nlte4), P2b diagnosis (a)-(d) and the cbfix
  check, prodcand confirmation.
- Running at this hand-back: pinc 4x GL directions (FfG4, FfG4R0; sptests build gdhr2), staircase exact at 96 x 96.
- Open:
  - vet_col_source = relaxed in a coupled-gas electron-scattering envelope (production): is E* ~ J there? Test: nlte2
    with coupled gas.
  - the fse+hc channel excess: a consistency limiter (M1 thin-cell flux vs the FS flux when E >> J_fs) would be the next
    code experiment;
  - Picard 30 at c dt/dx 3000.
- Not started: Stage-3 exact-solution tests (point/sphere source, shadow, inner-boundary beams, spherical Marshak,
  radial waves, 0-D relaxation), the cell/face and hesdirk2 arms (fse+hc only, per the scope change).
