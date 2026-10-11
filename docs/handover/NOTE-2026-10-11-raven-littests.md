# NOTE littests (littests-1011): the literature radiation tests, reproduced in AthenaK rad_m1 (implicit, vet_sc)

Work dir: /raven/ptmp/jinma/littests_1011. Branch `littests-1011` (local, NOT pushed) in the worktree
/raven/u/jinma/ATHENAK/wt-littests, on origin/vetacc-1011 (bb1ebff8).

**Bottom line.**
- **What was built.** The IDORT / J21 / JSD12 tests were reproduced in AthenaK with closure = vet_sc, which needed three small additions (Sect. 3):
  - a boundary-injected discrete-ordinate beam (`implicit_beam_n`);
  - the IDORT 2-D ring quadrature (`vet_ring`);
  - the IDORT shadow cloud (`shadow_dsq/peq`).
- **What passes.** Everything that is diffusive, local or a static shadow:
  - thermal equilibrium cases 1 and 2;
  - the non-LTE atmosphere, ε = 1e-1 to 1e-8;
  - static diffusion;
  - linear waves (2nd order for cen and tJe);
  - Marshak;
  - the static shadow, where E is within 1.1 % (L1) of exact for hr and tJe, against 8-12 % for cen.
- **What fails.** **The literature's own beam test, crossing beams in vacuum, fails for every arm.**
  - The formal solution J_sc is exact there.
  - The moment E does not follow it. hr and tJe lose the beam within about 5-10 cells of the boundary (row energy below 1e-3 by y = -1, 0 at the top); this holds at c dt/dx = 0.5 too. cen keeps 72 % of the energy but smears it over the whole box (beam peak 0.07 E_b on a 0.018 floor).
  - This is the face-flux defect the old beam tests could not isolate.
- **Two physics-mode problems** affect the old tests and the "static" ones:
  - (1) `gas_feedback = false` is not a static medium inside the implicit solve: the gas heats or cools within the step.
  - (2) with `gas_feedback = false` a moving medium loses the radiation work term (dynamic diffusion moves at 4/3 v).
- **Moving thermal equilibrium (v/c = 0.3) fails:** the gas is not decelerated (v 3.001 against 2.956).

Figures: /raven/ptmp/jinma/littests_1011/fig/*.png. Commits (local): 6c5a7260, e64b93a7, e2622521, 5586161e.

---

## 0. Papers

| short | paper | arXiv | role |
|---|---|---|---|
| IDORT | J.-Z. Ma, Pakmor, Justham & de Mink, AREPO-IDORT (implicit discrete ordinates RT), submitted to A&A ("Ma et al. 2025/2026") | 2503.16627 | **primary spec** (the user's paper) |
| J21 | Jiang 2021, ApJS 253, 49 (implicit discrete ordinates, Athena++) | 2102.02212 | tests IDORT copies, plus static diffusion, radiation shocks, dusty wind |
| JSD12 | Jiang, Stone & Davis 2012, ApJS 199, 14 (implicit VET, Athena) | 1201.2223 | the VET ancestor: thermal relaxation, Marshak, crooked pipe, ±14° shadow, shocks, photon bubble |
| DSJ12 | Davis, Stone & Jiang 2012, ApJS (short-characteristics VET; journal page not checked) | 1201.2222 | the formal-solution beam test (shallow beams; linear vs quadratic interpolation) |

The "Ma 2026 blue-supergiant" paper in the handover note (arXiv 2609.28656) is by **Linhao Ma**, Bildsten & Y.-F. Jiang,
not by the user. It uses the J21 scheme and has **no** transport tests.
The user's paper is IDORT (confirmed by the coordinator). I read every source from arXiv LaTeX.
The figure colour scales were read from the rendered figures.

## 1. Catalogue (Part 1)

Notation: c dt/dx is the radiation CFL. "(mine)" marks numbers I derived; they are not printed in the papers.
Cheapest first within each paper's own order.

### 1.1 IDORT (Ma et al., arXiv 2503.16627), in paper order

| § / Fig | test | setup | angles | c, dt | plotted / what is "correct" |
|---|---|---|---|---|---|
| 4.1.1 / Fig 2 | **crossing beams in vacuum** | Box [-0.5,0.5]x[-2,2] at 64x256 (and 256x1024), honeycomb mesh. kappa = 0. Two beams enter through one ghost cell of the bottom boundary: at x = -0.1 along (-1,1)/√2, and at x = +0.1 along (1,1)/√2, with I = 0.8c. Periodic in x, outflow at the top | 4 in-plane rays (±45°, ±135°) | a = 1, c = 1e3 ("do not matter"). Steady mode, 20 / 200 / 2000 iterations | E map, log scale ~3e-2 to 0.8. **Correct:** straight 1-cell beams that wrap periodically and pass through each other. Crossings (mine): x = 0 at y = -1.1, -0.1, 0.9, 1.9; x = ±0.5 at y = -1.6, -0.6, 0.4, 1.4. Spreading is numerical only; it is much smaller at 256x1024 |
| 4.1.2 / Fig 3 | homogeneous sphere | 1-D spherical, 1000 cells on r in [0.05, 7]. R0 = 1, ρκ = 10 inside, 0 outside, T0 = 1. Inner boundary fixed isotropic, outer outflow | 40 (J21 spherical set) | a = 1, c = 100. About 1000 iterations | E, F, P/E against r. Smit+97 analytic: f = 1/3 inside, sharp change at r = 1, F = cE ∝ r^-2 outside |
| 4.1.3 / Fig 4 | **non-LTE scattering atmosphere** | 1-D, 1280 cells on [-10,10]. ρ = 1e-3 exp(10 - x), T = 1, κ_a = ε, κ_s = 1 - ε, ε = 1e-1 ... 1e-8. Inner boundary isotropic B, outer outflow | 8 | a = c = 1; steady (1e3 to 1e6 iterations) | E against τ. Analytic: E = 1 - exp(-√(3ε) τ)/(1 + √ε), with τ = 1e-3 (e^(10-x) - 1) |
| 4.1.4 / Fig 5 | **static shadow** | [-0.5,0.5]x[-0.3,0.3] at 512x256. ρ = 1 + 9/{1 + exp[10((x/0.1)^2 + (y/0.06)^2 - 1)]}, P = ρT = 1, κ = ρ T^-3.5 (σ = ρ^5.5). Every left ghost cell carries I = 1031.3 c a along ±15°. Periodic y, outflow on the right | 24 (two rings) | a = 1e-7, c = 1.9e5; ≤ 2000 iterations | E linear 0-400, f_xx 0.2-0.9, f_yy 0.04-0.43. **Correct:** umbra tip at x ≈ 0.245 (mine). At x = 0.5: outer penumbra edge \|y\| ≈ 0.200, lit inner wedge \|y\| < 0.068. Wings carry half of the injected E. "Almost identical to J21" |
| 4.2.1 / Fig 6 | **dynamic diffusion** | 1280x16 on [-10,10]x[-.125,.125], periodic. ρ = 1, v = 1, κ_s = 4e4. E = exp(-40x^2), floor e^-10 | 8 | c = 1e3; 3-level local time stepping; α = 20 needed | E(x) at t = 4, 8, 16 against (160Dt+1)^-1/2 exp[-40(x-vt)^2/(160Dt+1)], D = c/(3ρκ) |
| 4.2.2 / Fig 7 | multi-source bubbles | [-0.5,0.5]^2 at 1024^2. 20 opaque discs; positions come from the Zenodo package | 200 | c = 10, dt_max = 1e-3 (≈ 25x CFL) | log E at t = 0, 0.05, 0.25, 1 t_c. Fronts at R0 + ct; bubbles pass through each other |
| 4.3.1 / Fig 8 | **thermal equilibrium** | 32^2, periodic, ρ = 1, P = ρT. Cases (T0, E0, κ_a, c): (1, 100, 100, 100); (100, 1, 1, 100); (1, 1, 1, 10) with v_x = 3 | 24 | - | E, T(t). Equilibria (mine): T = 3.137, E = 96.8; T = 3.475, E = 145.8. Moving case (J21): E → 1.13, v → 2.956 |
| 4.3.2 / Figs 9-10 | **linear waves** | [0,1]^2 at 16^2-256^2, ρ0 = T0 = 1, γ = 5/3, amplitude 1e-3. τ = ρκ_a in {0.01, 1, 100}; P = aT^4/P0 in {0.01, 10} | 8 | c = 10 | L1(ρ) ∝ N^-2 down to a 1e-6 floor; phase speed and damping at 256 match the analytic values |
| 4.3.3 / Figs 11-13 | radiation shocks | 1-D, 4096 cells, L = 0.115847, M = 1.2, 2, 3. κ = 577.4/ρ | 24 | a = 1e-4, c = 1.73e3 | profiles hold over 3 flow-crossing times against Ferguson+17 |
| 4.3.4 / Fig 14 | **irradiated cloud (RHD)** | the 4.1.4 shadow, relaxed, then evolved with hydro and AMR | 24 | dt_max = 3e-4 | ρ at t = 0.06, 0.1, 0.16 (t0 = sound crossing). Compare J21 Fig 5 |

### 1.2 J21 (arXiv 2102.02212): extra tests and differences from IDORT

- **4.1 Thermal equilibrium:** same three cases, but with ℂ = 100 printed for the moving case (IDORT says this is a typo for 10). dt ≈ 1e-3.
- **4.2 Crossing beams:** 4 angles at μ = ±0.577 (one per quadrant), I = 0.8, 64x256. dt = 4e-3, i.e. c dt/dx = 256 (mine). Log colour 0.05-1.9. Also run on SMR.
- **4.3 Shadow:** 40 angles; beam I = 103.13 (it is 1031.3 in IDORT; neither matches the plotted ~400, see F1 below). Dynamic cloud: ℙ = 1e-7, snapshots at 0.06, 0.1, 0.16 t0.
- **4.4 Static diffusion:** 1-D, 256 cells on (-1,1), outflow, ρκ_s = 4e4, ℂ = 10. Snapshots at √(Dt) = 0.13, 0.18, 0.22.
- **4.4 Dynamic diffusion:** as in IDORT; 1280 cells, ℂ = 1e3, t = 4, 8, 16.
- **4.5 Non-LTE atmosphere:** as in IDORT, with 1 angle per octant.
- **4.6 Homogeneous sphere:** ρ = 1e-7 outside.
- **4.7 Linear waves:** first order where RT dominates.
- **4.8 Radiation shocks:** M = 2, 3; ℙ = 7.72e-4, ℂ = 1.73e3.
- **4.9 Dusty wind:** 512x1024, 2e4 steps.

### 1.3 JSD12 (arXiv 1201.2223): extra tests

- **§4 Thermal relaxation:** 1-D, ℂ = 1e4, σ_a = 100, (E, T) = (100, 1) and (1, 100). Equivalent to IDORT cases 1 and 2.
- **§5.1 Marshak wave:** L = 0.1, 512 cells, c dt/dx = 1. Parameters are those of Sekora & Stone 2010 (not read). Compared with Su & Olson.
- **§5.2 Crooked pipe / tophat:** 512x256 on (0,7)x(-2,2). Opaque ρ = 10, σ_a = 2000; pipe ρ = 0.01. Isotropic 0.5 keV source on the left for \|y\| < 0.5. Probe histories against Gentile 2001.
- **§6.1 Linear waves:** 1-D, ℂ = 1e4, f = 1/3.
- **§6.2-6.3 Shocks:** Lowrie-Edwards; Sincell VET shock.
- **§6.4 Photon bubble.**
- **§6.5 Shadow at ±14°:** injected also through the top and bottom boundaries ("infinite plane"); y outflow; ℂ = 1.9e5. One step with dt = 1e-3 (≈ 500 light crossings). The number of VET angles is not stated. Colour E 0-200, f_xx 0.253-0.938. Geometry (mine): umbra tip at x = 0.261; at x = 0.5 the outer edge is at \|y\| = 0.190 and the lit inner wedge spans \|y\| < 0.060.
- **DSJ12 §5.2 beam test:** two shallow beams (n·y = 0.174) that wrap several times. Linear interpolation is much more diffusive than monotone quadratic interpolation. Resolution is not stated.

**Flags:**
- (F1) The beam intensity is printed differently in IDORT and J21 (1031.3ca vs 103.13), and neither matches the plotted level of about 400 aT0^4. I inject 200 per beam, so the full-face plateau is 400.
- (F2) The background optical depth (≈ 1) cannot be checked in the papers' saturated colour maps.
- (F3) The crossing-beam width is one ghost cell, inferred.

### 1.4 Status per test (Part 3 summary)

| test (paper §) | status | AthenaK set-up / reason |
|---|---|---|
| thermal equilibrium 1, 2 (IDORT 4.3.1, J21 4.1, JSD12 §4) | **reproduced, PASS** | `equil`, 32^2 periodic, 24 → Gauss 8x32 rays |
| moving thermal equilibrium (v = 3, c = 10) | **reproduced, FAIL** | radiation momentum F/c^2 = 0.041 is created but the gas is not decelerated |
| crossing beams (IDORT 4.1.1, J21 4.2) | **reproduced, FAIL** (all arms) | new `implicit_beam_n` + `vet_ring = 4`, rotated (paper y → x1) |
| static shadow (IDORT 4.1.4, J21 4.3; JSD12 6.5 at ±15°) | **reproduced, PASS** (hr, tJe); cen marginal | new `shadow_dsq/peq` + beams + `vet_ring = 12/36/60` x 2 groups. JSD12's ±14° set is not representable (a ring containing ±14° also contains a ray along x2); its top/bottom injection is equivalent to periodic y |
| non-LTE atmosphere (IDORT 4.1.3, J21 4.5) | **reproduced, PASS** (hr, tJe); cen surface +15 % | `atmosphere` pgen, as a 1280x4 slab |
| static diffusion (J21 4.4) | **reproduced, PASS** | `thick_pulse`, 256x4 |
| dynamic diffusion (IDORT 4.2.1, J21 4.4) | **reproduced; FAIL with gas_feedback = false (moves at 4/3 v); with feedback the position is right, peak 13 % low** | new `gas_vx` |
| linear waves (IDORT 4.3.2, J21 4.7) | **reproduced, PASS** (cen, tJe 2nd order); hr breaks at P = 10, τ = 1 | `radwave` with the exact eigenmode, 1-D on 2-D, N = 16-256 |
| Marshak (JSD12 5.1) | **PASS** (existing S_N-referenced T6 gate, not SS10's parameters) | `marshak`, 128 cells, c dt/dx = 1, 10 |
| irradiated cloud (IDORT 4.3.4, J21 Fig 5) | **reproduced qualitatively** | shadow + hydro; no quantitative reference |
| homogeneous sphere (IDORT 4.1.2, J21 4.6) | **not done** | needs a spherical-polar run; vet_sc is Cartesian only (the spherical closures are vet_col). A Cartesian disk analogue is possible, not done |
| multi-source bubbles (IDORT 4.2.2) | **not done** | source positions are only in the Zenodo package. The `cyl` pgen allows 8 objects (IDORT has 20); needs an outflow x2, but vet_sc needs periodic x2 |
| radiation shocks (IDORT 4.3.3, J21 4.8, JSD12 6.2-6.3) | **not done** (time) | `radshock` pgen exists but has explicit ghost BCs; it needs an implicit Dirichlet (efix) set-up plus a Ferguson/Lowrie-Edwards reference; about 1e4 steps at 4096 cells |
| crooked pipe (JSD12 5.2) | **not done** | needs an isotropic patch bath on x1min (a variant of `implicit_beam`), keV units and c_v = 1e15 |
| photon bubble, dusty wind, Sincell shock | **not feasible in short jobs** | 1e3-2e4 steps at 256^2-512x1024, and MHD |


---

## 2. Audit of our current beam tests (Part 2)

Tests audited:
- Inputs: /raven/ptmp/jinma/thinfix_1011/cpu/inp/{ba0,ba20,xb20,cyl,shd3b}.athinput.
- Runs: cpu/run/beams/*.
- References: xtf/ana collref.py, cylref.py, shdfs.py.
- Figures: /raven/ptmp/jinma/beamplots_1011.

All five tests are `m1_test = cyl` set-ups: opaque objects in a thin background (ρ = 1e-4, T = 1e-3), Cartesian 2-D, 128x256 (ba*, xb20) or 128x384 (shd3b). They use c = 1e3, c dt/dx ≈ 300, 100 steps, and Gauss 8x32 rays.

| issue | our tests | literature (IDORT 4.1.1/4.1.4, J21 4.2/4.3, JSD12 6.5, DSJ12 5.2) | consequence |
|---|---|---|---|
| **Injection** | No boundary beam. A hot emitter (ba*/xb20: a 5x5-cell rectangle, ρ = 128, κ = 1, τ_cell = 1, T held at 1) inside a collimator channel made of two absorbing walls (ρ = 1e4, τ_cell = 78, 3 cells thick), pointed at 0° or 20°. cyl/shd3b: emitting disks | Intensity on ONE discrete ray, set in the ghost cells of a boundary. One cell wide (crossing beams) or the whole face (shadow). The direction is a member of the angle set | The beam's angular spread and intensity are set by the channel geometry and by the under-resolved emitter limb, not by an input. The "beam" is a cone of about ±4.6° (5-cell channel, 0.3 long), so the beam width at x = 0.5 is physical spreading, not numerical diffusion |
| **Emitter energetics** | T is re-imposed after every stage, but the implicit E-solve eliminates the gas locally. **Inside each step the emitter cools.** Per step it emits c κ ρ a T^4 dt ≈ 300, about 1.6x its thermal energy ρT/(γ-1) = 192 | Beams are boundary intensities; the gas is fixed | The E-solve sees a cooler emitter than the formal solution and the "exact" reference do. Rerun with γ = 1.000001 (heat capacity x 1e6, everything else identical): beam E/J_exact goes from 1.08 / 0.49 / 0.75 (ba0 / ba20 / xb20, tJe; this reproduces the NOTE-thinfix numbers exactly) to **1.85 / 0.87 / 1.34**; cen goes from 7.97 / 1.52 / 1.79 to 12.7 / 2.30 / 2.60. So about 40 % of the reported tJe beam deficit was emitter cooling, not transport. The same mechanism makes `gas_feedback = false` NOT a static medium inside the solve (Sect. 4.4) |
| **"Exact" reference** | collref/cylref: opaque surfaces at the geometric boundary (τ → ∞), I = aT^4 there, continuous angles (NPHI = 4096 in-plane azimuths, 48 polar nodes), infinite in z. **Evaluated only for x ≥ 0.45** (XLO), beyond the channel mouth | Analytic straight rays from a boundary intensity that is part of the discretised problem | Not exact for the problem the code solves: (a) the emitter limb has τ_cell = 1, not a sharp surface; (b) the emitter cools inside the step (above); (c) continuous angles against 512 discrete rays. Ray effects ("spokes", DSJ12) are therefore counted as errors of the moment scheme. (d) The channel interior and the emitter are not checked at all (the blank panel x < 0.45) |
| **What is tested** | Moment E against exact J. The formal solution itself already misses: J_sc/J_ex = 0.65 / 0.15 / 0.21 (ba0 / ba20 / xb20). Moment E/J_sc = 1.45 / 3.36 / 3.83 | IDORT/J21 compare the transport solution with the analytic one. With beams on grid diagonals (45°, square cells) the formal solution is exact | **"E ≈ J_exact" in the old tests was two compensating errors:** the SC loses 35-85 % of a 20° beam, and the moment E overshoots J_sc by 1.5-3.8x (6.8x with a static emitter). The face flux cannot be judged from E/J_exact. In the literature geometry J_sc is exact (Sect. 4.1: J_sc L1 = 0 for the 4-ray crossing beams, 1.1 % for the shadow), so every remaining error is the moment solve |
| **Angles** | Gauss 8x32 about x1. No in-plane ray at 20° or 45°: the nearest in-plane directions are 11.4°, 26°, 40.3°, 53.7°. Beam directions are not members of the set | IDORT: in-plane rings (4 for crossing beams, 2x12 for the shadow); J21: 4 quadrant rays; the beam direction is always in the set | A beam not on a ray is smeared by the angular interpolation of the SC sweep. The strong ray patterns in the J_sc maps are the 32 azimuths of a compact source (DSJ12's "spokes"), not noise |
| **Resolution** | 128 cells across, 5-cell beams, 20° | 64x256 and 256x1024 (crossing beams), 512x256 (shadow) | Comparable; resolution is not the issue |
| **c dt/dx** | ≈ 300 | J21 crossing beams 256 (mine); IDORT/J21 shadow in steady mode (≫ 1); JSD12 shadow one step of ≈ 500 crossings | Fine; the literature tests are also in the implicit regime |
| **Diagnostic of the face flux?** | Weakly. Emitter cooling, channel collimation, SC smearing, ray effects and reference error all enter at the 30-80 % level | The crossing beams with 4 rays isolate the moment face flux completely; the shadow does too with ring angles | Use the literature geometries (Sect. 5) |
| **Other** | shd3b: 6 absorbing disks along the bottom plus an emitter and an absorber. The shadow edges are judged at the 10-90 % level from a reference with continuous angles | IDORT/J21: geometric umbra/penumbra of ±15° beams | The geometric shadow edges (exact field, mid-level crossings): x = 0.4: inner wedge \|y\| = 0.032, outer edge 0.179; x = 0.5: 0.060 and 0.205. With ring angles the code's J_sc and E get the plateaus to 0.3 % and the edge positions to within one cell |

**About the figures in /raven/ptmp/jinma/beamplots_1011:**
- The exact panel covers only x ≥ 0.45 by construction (collref XLO).
- The ray patterns in the J_sc maps are the discrete-ordinate ray effect of a small hot source seen with 32 azimuths. They are not a bug; DSJ12 Sect. 5.2 warns about exactly this.
- In shd3b, every arm and J_sc sit about 0.6x below "exact" in the lit region (cut at x = 0.65). That is the same SC/angle loss, so the 0.184 / 0.295 "exact" edge widths cannot be reached by any face flux.


---

## 3. What was added to AthenaK (branch littests-1011, local)

Commits, in order (all local on `littests-1011`; NOT pushed):

**6c5a7260 `implicit_beam_n`** (rad_m1.hpp, rad_m1_vet.cpp, rad_m1_implicit.cpp; read only when given, so the default is bitwise unchanged).
- Up to 2 beams through the x1min face: `implicit_beam<b>_{e, ang, y0, y1}` give the energy density, the in-plane angle (degrees) and the x2 patch.
- Each beam is put on ONE up-going SC ray; in 2-D also on its mu3 mirror, so F3 = 0.
  - The ray's entering intensity is frac E_b/(n w_r), where frac is the patch coverage of the cell. It is added in the VetShortChar boundary branch.
  - The moment solve takes the same rays' face flux, F_in = c Σ w μ1 I, as an imposed flux on those columns (Marshak elsewhere), both in the row and in the face-flux kernel.
- Requires vet_sc, implicit_bc_x1min = marshak, and a single MeshBlock.

**e64b93a7 `vet_ring = N`** (multiple of 4).
- The IDORT/J21 2-D ring: θ = (l + 1/2) 2π/N, equal weights, with ray r + nray/2 the x1 mirror of ray r.
- `vet_ring_groups = 2` adds the tilted ring sqrt(1/3)(cos, sin), ±sqrt(2/3), so f = 1/3 for an isotropic field.
- N = 4 gives exactly ±45° and ±135°; N = 12, 36, 60 contain ±15° and ±45°.
- Same commit, shadow pgen: `shadow_dsq` (the D^2 exponent of IDORT/J21) and `shadow_peq` (T = t0 ρ0/ρ, pressure equilibrium).

**e2622521 `gas_vx`** for the `thick_pulse` and `equil` pgens (dynamic diffusion, moving equilibrium).

**5586161e** The beam ray is the one nearest the beam direction in 3-D. The first version picked rays by in-plane angle and could choose a ray with mu3 = 0.99.

**Binaries** (in /raven/ptmp/jinma/littests_1011/bin, built with gcc/13 openmpi/5.0, OpenMP):
- athena_cpu_e64b93a7: crossing beams, ring shadow, figures.
- athena_cpu_e2622521: stage 2.
- athena_cpu_5586161e: Gauss-angle beam/shadow reruns.

The vet_sc / face-flux code is identical in all three to 42aa6bbc (vetacc-1011).

**Run-mode settings that were necessary** (not code):
- `hydro/gamma = 1.000001` for the static shadow and the beams. It makes the gas heat capacity effectively infinite (see 4.4).
- `time/evolution = kinematic` + `hydro/rsolver = advect` for static media that need T fixed and a large dt (non-LTE, static diffusion). Otherwise the hydro sound-speed CFL sets dt.


---

## 4. Reproduced results (Part 3)

Arms, with keys exactly as in xtf/arms/battery_arms.sh and NOTE-thinfix:
- **cen:** `closure=vet_sc vet_tensor=full implicit_flux=central`.
- **hr:** + `implicit_flux_faces=all implicit_flux_beam=halfrange implicit_flux=blend implicit_blend=idort_f tau0=1 alpha=1 flo=0.3 fhi=0.6 implicit_blend_xthin=30`.
- **tJe:** hr keys with `implicit_blend_xthin=0 implicit_thin_switch=fse implicit_thin_corr=hcs implicit_thin_recon=jse implicit_thin_bq=lin`.

All runs use backward Euler, CPU, one MeshBlock and OpenMP.

### 4.1 Crossing beams in vacuum (IDORT Fig 2 / J21 Fig 3). FAIL, all arms. Figure: fig/cb_maps.png

**Set-up** (runs/cb*):
- Rotated so that the beam axis is the SC sweep axis: paper y → x1 ∈ [-2, 2] (256 cells), paper x → x2 ∈ [-0.5, 0.5] (64 cells, periodic).
- Beams of E_b = 1 through one cell at x2 = ∓0.1, at ∓45°.
- ρ = 1e-8, c = 1e3, c dt/dx = 256 (J21's dt = 4e-3), 20 steps (20 box crossings).
- Exact = cell-averaged straight strips.
- Metric: beam energy per row (Σ_x E over exact; 1 means energy conserved along the beam), beam peak and FWHM at y = 0.65.

| arm / angles | row energy at y = -1.5 / 0 / 1.9 | peak (exact 0.5-1) | FWHM (cells) | background E | J_sc row energy, peak |
|---|---|---|---|---|---|
| cen, ring 4 | 0.74 / 0.74 / 0.74 | 0.074 | 3.0 | 0.019 | 1.000, 1.000 |
| hr, ring 4 | 0.000 / 0 / 0 | 0 | - | 3e-12 | 1.000, 1.000 |
| tJe, ring 4 | 0.017 / 0 / 0 | 0 | - | 8e-7 | 1.000, 1.000 |
| tJe, ring 12 / 36 | 0.049 / 0.101 at y = -1.5, 0 beyond | 0 | - | - | 1.000 |
| cen / hr / tJe, ring 4, **c dt/dx = 0.5** (1100 steps) | 0.61 / 0.000 / 0.017 at y = -1.5 | 0.072 / 0 / 0 | | | 1.000 |
| same at c dt/dx = 4 / 32 | cen 0.61 / 0.69 at y = -1.5 (0.72 at the top); hr 0; tJe 0.017 then 0 | as above | | | 1.000 |
| cen / hr / tJe, ring 4, 256x1024 | cen 0.72 (all of it in a uniform floor, E max 0.039); hr and tJe 0 | 0.04 / 0 / 0 | | 6e-3 / 0 / 0 | 1.000 |
| cen / hr / tJe, **Gauss 8x32** (production; the beam goes on the nearest ray, 40.2°) | 2.7-5.2 / 0.29-0.43 / 0.24 | 0.14 / 0.02 / 0.02 | - / 15 / 11 | 0.14 / 9e-3 / 2e-3 | 1.000, 0.084 |

**Reading:**
- The 4-ray formal solution is **exact**: the beams sit on cell diagonals and bilinear interpolation is exact there. So every error is the moment solve.
- hr and tJe **lose the beam within 5-10 cells**: row energy drops by e every ~5-7 rows (see the right-hand panel of the figure).
- The energy does not go anywhere physical. Near the boundary E piles up (row 0 holds 2.7x the injected), spreads sideways into cells where J_sc = 0, and leaves back through the x1min Marshak face.
  - Mechanism, read off the dumps: the half-range face coefficients come from the SC half-range sums.
  - In a cell with J_sc = 0 next to a 1-cell beam, the face has "no formal solution" (M1HrFace sets the face to central). The E that numerical diffusion puts there therefore has no upwind outflow, and in the transparent box the central coefficient c^2 dt/dx^2 is huge.
  - This happens at c dt/dx = 0.5 too, so **it is not an implicit-regime effect**.
- cen keeps the right row energy only in the average. It is the FLD-like behaviour IDORT warns about: the box fills to 0.019 E_b everywhere and the beams sit on top at 0.07.
- With the production Gauss angles, the beam rides a ray 4.8° off and J_sc spreads it to 8 % peak. cen then **creates** energy (row energy 2.7-5x): the E-only operator piles up at c dt/dx = 256 with no outflow.
- The literature (IDORT, J21 first order) shows beams that broaden slowly but keep their energy and cross cleanly; both are 1-cell beams with peak about 0.5-0.8.

### 4.2 Static shadow (IDORT Fig 5 / J21 Fig 4). PASS for hr and tJe. Figures: fig/shd_r12g2.png, fig/shd_g8x32.png

**Set-up** (runs/shd*):
- 512x256, ±15° beams over the whole x1min face, 200 each (400 total).
- σ = ρ^5.5, ρ from the IDORT D^2 tanh ellipse, P = 1.
- c = 1.9e5, c dt/dx = 1e4, 10 steps.
- **γ = 1.000001** so the gas is static inside the solve (Sect. 4.4).
- Exact = attenuated straight beams, integrated along the ray with 8 sub-steps per cell. It contains no thermal part (S = ρ^-4 ≤ 1 against 400).
- Levels at x = 0.4: exact lit 157.6, wing 78.8, inner wedge 157.6, umbra 0.

| arm / angles | L1(E - exact)/exact | lit | wing | wedge | umbra at x = 0.15 (exact 0) | outer edge 10-90 % (cells; exact 2) |
|---|---|---|---|---|---|---|
| J_sc (any arm), ring 12x2 | 0.011 | 158.1 | 79.3 | 157.1 | 0.5 | 16 |
| **cen**, ring 12x2 / 36x2 / 60x2 | 0.082 / 0.121 / 0.082 | 171.5 | 142.9 | 107.2 | 36.9 | 41 |
| **hr**, ring 12x2 (36, 60 the same) | 0.012 | 158.2 | 79.4 | 150.9 | 1.4 | 16 |
| **tJe**, ring 12x2 (36, 60 the same) | 0.011 | 158.2 | 79.3 | 165.3 | 0.22 | 16 |
| cen / hr / tJe, Gauss 8x32 (the beam goes on the 11.4° ray) | 0.099 / 0.029 / 0.030 | 186 / 160 / 160 | 150 / 80.5 / 80.4 | 79 / 121 / 108 | 23 / 1.0 / 0.2 | 40 / 15 / 14 |

Edge positions at mid-level, x = 0.4: exact 0.032 / 0.179; tJe 0.029 / 0.179; J_sc 0.032 / 0.182. At x = 0.5: exact 0.060 / 0.205; code 0.057 / 0.207. That is within one cell (dy = 0.0023).

**Reading:**
- With angle sets that contain the beam direction, hr and tJe reproduce the paper's shadow: umbra, penumbra and the half-intensity wings, with L1 = 1.1 %, the same as the formal solution itself.
- The remaining error is the SC bilinear smearing of a 15° beam: the edge is 16 cells wide instead of 2. That is the DSJ12 linear-interpolation diffusion; IDORT's DOM upwind is similar.
- tJe overshoots the inner wedge by 5 % (165 against 157.6); hr undershoots by 4 %.
- **cen fills the umbra** (37 against 0) and smears the edges over 41 cells. That is the J21 "FLD would fill the shadow" failure.
- The production Gauss set cannot carry a ±15° beam (the nearest ray is at 11.4°). Everything then shifts, by L1 3 % for hr and tJe.
- **Before the γ fix** (γ = 5/3 with `gas_feedback = false`), all arms gave E 2.5x J_sc in the lit region and a spike of 2700 in the umbra (L1 = 1-3.6). See 4.4.

### 4.3 Old beam tests rerun with a static emitter (runs/old_*; NOTE-thinfix inputs, nlim 100)

| test | J_sc/J_ex | tJe E/J_ex, γ = 5/3 → 1.000001 | cen E/J_ex, γ = 5/3 → 1.000001 | tJe E/J_sc (γ = 5/3 → 1.000001) |
|---|---|---|---|---|
| ba0 | 0.655 | 1.083 → **1.848** | 7.97 → 12.7 | 1.45 → 2.47 |
| ba20 | 0.151 | 0.485 → **0.874** | 1.52 → 2.30 | 3.36 → 6.06 |
| xb20 | 0.209 | 0.754 → **1.344** | 1.79 → 2.60 | 3.83 → 6.83 |

The γ = 5/3 tJe numbers reproduce NOTE-thinfix exactly. See Sect. 2.

### 4.4 Two mode problems found on the way (not face-flux issues)

1. **`gas_feedback = false` is not "static gas" inside the implicit solve.**
   - The row eliminates the gas locally, so the absorbed radiation heats the gas within the step; the heating is then discarded.
   - Slab test: a 15° beam, σ = 1 (runs/shd/*_slab*). With γ = 5/3, E/exact grows from 1.0 to 2.45 across τ = 1 (cen: 2.8 to 4.7). With γ = 1.000001 it is 0.999-1.008 for hr and tJe, and 1.09-1.31 for cen.
   - J_sc is right in both cases, because the SC uses the step-start T.
   - Any test with a held-temperature object (all `cyl` tests, the atmosphere at large dt) is affected.
   - Suggested fix (not done): with gas_feedback = false, freeze T in the row (de/dT → ∞).
2. **Dynamic diffusion with `gas_feedback = false` advects E at (4/3) v.**
   - Centroid at t = 4 / 16 is 5.33 / 21.33 (wrapped), against v t = 4 / 16; L1 ≈ 2.
   - The comoving work term v·∇P is dropped with the feedback. With `gas_feedback = true` the centroid is exact (4.0009 / -3.9999 wrapped).

### 4.5 One-dimensional and local tests

| test | arms | result | literature target | verdict |
|---|---|---|---|---|
| thermal eq. 1 (T, E) = (1, 100), κ = 100, c = 100, dt = t_th/20 | all identical | E_final 96.795, T 3.13663 (exact 96.795, 3.13663); max \|E/E_ex - 1\| over time 1.1e-5, \|T/T_ex - 1\| 2.3e-4. At J21's dt = 1e-3 (one step): 96.799 / 3.134 | equilibrium (mine) 96.8 / 3.137 | **PASS** |
| thermal eq. 2 (100, 1), κ = 1 | all identical (hr shows a spurious v_x = 0.021) | 145.788 / 3.4747 (exact 145.788 / 3.4748) | 145.8 / 3.475 | **PASS** (hr: spurious velocity in a uniform medium, not understood) |
| thermal eq. 3, v = 3, c = 10 | all (and implicit_vimp) | E 1.0325, T 0.976, v 3.0011, F 4.14 (≈ 4/3 v E) | E 1.13, v 2.956 (J21) | **FAIL.** The radiation momentum F/c^2 = 0.041 is created but not taken from the gas, so 0.13 of kinetic energy is never converted |
| non-LTE atmosphere, ε = 1e-1 ... 1e-8, 1280x4, Gauss 8x32, c dt/dx = 1e6, 60 steps (200 gives identical results), kinematic gas | hr, tJe | L1 1.1 % (ε = 0.1) down to 0.01 % (ε = 1e-8). Surface E/E_an = 1.000-1.005. Max relative deviation 10 % at τ ≈ 0.3 for every ε and with nmu = 2. This is the Eddington analytic solution against exact transfer | DOM matches the analytic at the J21/IDORT plot level | **PASS** |
| | cen | surface E 6-15 % high, max 23 % | | marginal |
| static diffusion, 256x4, c = 10, ρκ_s = 4e4, √(Dt) = 0.13, 0.18, 0.22, kinematic gas | all identical | L1 0.10 / 0.11 / 0.32 %, peak 1.001 | analytic Gaussian | **PASS** |
| dynamic diffusion, 1280x4, v = 1, c = 1e3, κ_s = 4e4, hydro CFL 0.4 (c dt/dx ≈ 170), gas_feedback = true | all identical | position exact; peak 0.871-0.875 of analytic at t = 4-16; L1 0.115 | analytic (J21 needs α = 20 in 2-D) | **PASS** for position, peak 13 % low (an extra diffusion set in by t = 4, constant afterwards) |
| linear waves, c = 10, 1-D on 2-D, N = 16-256, exact eigenmode, one period, nt = 2N (c dt/dx ≈ 3.9), Gauss 4x8 | cen, tJe | L1(ρ) order 1.9-2.1 for every (P, τ). The only exception is P = 10, τ = 1, with 1.77 then 1.30 at 256 (floor 3.8e-4). ω matches the quasi-static VET dispersion to < 1e-3 | IDORT: N^-2 down to 1e-6 | **PASS** |
| | hr | At P = 10 it does not converge. τ = 1: L1 = 3.2e-3 (cen) against **4.97e-3 / 4.92e-3 / 8.76e-3** for N = 64 / 128 / 256 (order -0.83). τ = 0.01: 2.24e-3 → 2.31e-3 (order -0.04). τ = 100: order 0.68 at 256. At P = 0.01 it is 2nd order | | **FAIL** at P = 10 (xthin = 30 is dt-dependent: X = c dt/dx/(1 + c dt χ) changes with N) |
| Marshak (T6 gate, S_N reference, 128 cells, κ = 300) | cen / hr / tJe at c dt/dx = 1 and 10 | L1(E) = 0.0135 / 0.0138 / 0.0139 (c dt/dx = 1); 0.0136 / 0.0140 / 0.0142 (10) | ≤ 0.02 | **PASS** |

Linear-wave table (ana/rw_table.txt; runs/rw/*/results.json). L1 is the error in ρ after one period, normalised by the amplitude, against the eigenmode of the code's continuum model (quasi-static VET, `vetqsc_d`). ω is shown at N = 256, or the largest N run. P = aT^4/P_gas (IDORT's ℙ):

```
P        tau    arm  | L1(N=16) L1(N=32) L1(N=64) L1(N=128) L1(N=256) | orders | w_re(256) w_im(256) | ref vetqsc_d | transp
0.01     0.01   cen  | 9.51e-02 3.80e-02 1.06e-02 2.79e-03 7.02e-04 | 1.32 1.84 1.93 1.99 | 8.11212 -0.00059 | 8.11156 -0.00055 | 8.11156 -0.00054
0.01     0.01   hr   | 9.51e-02 3.80e-02 1.06e-02 2.80e-03 7.04e-04 | 1.32 1.84 1.93 1.99 | 8.11213 -0.00059 | 8.11156 -0.00055 | 8.11156 -0.00054
0.01     0.01   tJe  | 9.51e-02 3.80e-02 1.06e-02 2.79e-03 7.02e-04 | 1.32 1.84 1.93 1.99 | 8.11212 -0.00059 | 8.11156 -0.00055 | 8.11156 -0.00054
0.01     1      cen  | 9.18e-02 3.63e-02 1.00e-02 2.49e-03 5.87e-04 | 1.34 1.86 2.01 2.09 | 8.11602 -0.04180 | 8.11559 -0.04173 | 8.11228 -0.04143
0.01     1      hr   | 9.18e-02 3.63e-02 1.00e-02 2.54e-03 6.74e-04 | 1.34 1.86 1.98 1.91 | 8.11630 -0.04167 | 8.11559 -0.04173 | 8.11228 -0.04143
0.01     1      tJe  | 9.18e-02 3.63e-02 1.00e-02 2.49e-03 5.87e-04 | 1.34 1.86 2.01 2.09 | 8.11602 -0.04180 | 8.11559 -0.04173 | 8.11228 -0.04143
0.01     100    cen  | 9.38e-02 3.67e-02 9.99e-03 2.46e-03 5.67e-04 | 1.35 1.88 2.02 2.12 | 8.11247 -0.00676 | 8.11190 -0.00673 | 8.11134 -0.00663
0.01     100    hr   | 9.38e-02 3.67e-02 9.99e-03 2.47e-03 5.73e-04 | 1.35 1.88 2.02 2.11 | 8.11248 -0.00677 | 8.11190 -0.00673 | 8.11134 -0.00663
0.01     100    tJe  | 9.38e-02 3.67e-02 9.99e-03 2.46e-03 5.67e-04 | 1.35 1.88 2.02 2.12 | 8.11247 -0.00676 | 8.11190 -0.00673 | 8.11134 -0.00663
10       0.01   cen  | 6.31e-02 2.16e-02 5.53e-03 1.31e-03 2.92e-04 | 1.55 1.97 2.08 2.16 | 7.98863 -0.53319 | 7.98801 -0.53314 | 7.99079 -0.52043
10       0.01   hr   | 6.32e-02 2.17e-02 5.93e-03 2.24e-03 2.31e-03 | 1.54 1.87 1.40 -0.04 | 7.99605 -0.53536 | 7.98801 -0.53314 | 7.99079 -0.52043
10       0.01   tJe  | 6.31e-02 2.16e-02 5.53e-03 1.31e-03 2.92e-04 | 1.55 1.97 2.08 2.16 | 7.98863 -0.53319 | 7.98801 -0.53314 | 7.99079 -0.52043
10       1      cen  | 3.90e-02 1.22e-02 3.20e-03 9.35e-04 3.80e-04 | 1.67 1.94 1.77 1.30 | 6.23523 -0.87618 | 6.23399 -0.87560 | 6.28768 -0.83695
10       1      hr   | 3.94e-02 1.28e-02 4.97e-03 4.92e-03 8.76e-03 | 1.62 1.37 0.01 -0.83 | 6.25279 -0.90489 | 6.23399 -0.87560 | 6.28768 -0.83695
10       1      tJe  | 3.90e-02 1.22e-02 3.20e-03 9.35e-04 3.80e-04 | 1.67 1.94 1.77 1.30 | 6.23523 -0.87618 | 6.23399 -0.87560 | 6.28768 -0.83695
10       100    cen  | 6.05e-02 1.51e-02 3.70e-03 9.14e-04 2.26e-04 | 2.00 2.03 2.02 2.01 | 15.31402 -0.43669 | 15.31514 -0.43666 | 14.38498 -0.44140
10       100    hr   | 6.04e-02 1.49e-02 3.38e-03 4.82e-04 3.01e-04 | 2.02 2.14 2.81 0.68 | 15.31655 -0.43696 | 15.31514 -0.43666 | 14.38498 -0.44140
10       100    tJe  | 6.05e-02 1.51e-02 3.70e-03 9.14e-04 2.26e-04 | 2.00 2.03 2.02 2.01 | 15.31402 -0.43669 | 15.31514 -0.43666 | 14.38498 -0.44140
```

Note: the exact time-dependent transfer root (`transp`) differs from the quasi-static VET root by up to 6 % in Re ω at P = 10, τ = 100. That is the model difference between VET and DOM, not a discretisation error. I did not check whether `transp` picked the right root there.

### 4.6 Irradiated cloud (IDORT Fig 14 / J21 Fig 5). Qualitative. Figure: fig/cloud_rho.png

**Set-up:**
- The shadow with hydro and gas_feedback = true, a = 1e-7, beams 2e-5 each.
- c = 1.9e5, dt = 3e-4 (IDORT's dt_max), 533 steps, 36-ray ring.
- Wall time: cen 476 s, hr 110 s, tJe 342 s.

**Results:**
- hr and tJe agree (centroid of ρ > 3 at x = 0.009 / 0.016 / 0.027 at t = 0.06 / 0.10 / 0.16).
- The lit face is flattened and pushed downstream; the shadowed side stays cold (min T 0.1).
- cen pushes the cloud 3x less (0.003 / 0.005 / 0.009) and heats it from all sides, because its umbra is filled.

There is no published number to compare; the morphology matches the papers' qualitatively. The background also heats (T up to 1.5), because τ ≈ 1 absorbs about half the beam.


---

## 5. Recommendation

1. **Replace ba0 / ba20 / xb20 / cyl / shd3b by the literature geometries.** All are cheap (seconds to a minute on 72 cores):
   - (a) **Crossing beams**, IDORT 4.1.1 / J21 4.2: 64x256, `vet_ring = 4`, `implicit_beam_n = 2`, c dt/dx = 256 (and 0.5).
     - J_sc is exact, so it measures the moment face flux alone.
     - Pass criterion: beam energy per row ≥ 0.9 of injected at the top, peak ≥ 0.3 E_b at y = 0.65, background ≤ 1e-3 E_b.
     - Every current arm fails. **This is the test the new face flux must pass.**
   - (b) **Static shadow**, IDORT 4.1.4 / J21 4.3: 512x256, `vet_ring = 12, vet_ring_groups = 2`, `hydro/gamma = 1.000001`.
     - Pass: L1 ≤ 2 % against the attenuated-beam exact field, umbra ≤ 1 % of lit, edge positions within 1 cell.
     - hr and tJe pass, cen fails.
   - (c) **The same two with the production Gauss 8x32 set**, as a robustness check. Here the beam direction is not representable, so compare with J_sc, not with exact.
2. **Keep the 1-D/local gates:** non-LTE atmosphere (kinematic gas), static and dynamic diffusion (with gas_feedback = true), thermal equilibrium 1-2, linear waves (cen/tJe second order), Marshak.
3. **Drop or fix the emitter-in-a-channel tests.** If an internal emitter is wanted, make the gas static inside the solve (γ → 1, or a code fix for gas_feedback = false). Compare with the code's own J_sc and with an exact reference computed on the code's ray set, not with continuous angles.
4. **Bugs and defects to follow up** (outside this task):
   - (i) the moment face flux in cells with J_sc = 0 next to a beam (4.1);
   - (ii) gas_feedback = false still heats the gas inside the solve (4.4.1);
   - (iii) the work term is lost with gas_feedback = false (4.4.2);
   - (iv) no gas deceleration in moving thermal equilibrium (4.5);
   - (v) hr's spurious velocity in a uniform medium (eq2) and its loss of convergence on linear waves at P = 10.


---

## 6. Files

- Inputs:
  - inp/cb.athinput, cbg, cbh (crossing beams); shd, shdg (shadow); cloud;
  - eq1-3, eq3v; nltek1-8 (nlte/nlteg: earlier, invalid modes); sdifk (sdif/sdifc: invalid, hydro-CFL bound); ddif, ddiff; marsh;
  - old_{ba0,ba20,xb20,shd3b}.
- Run scripts: run1.sh <test> <arm> <tag> [keys], arms.sh, pool.sbatch / poolg.sbatch (xargs pool), rw.sbatch; job lists jobs_s*.txt.
- Runs: runs/<test>/<arm>_<tag>/ (in.athinput, cmd.txt, run.log, time.txt, rc.txt, bin/ or tab/); runs/rw/* (linear waves).
- Analysis (ana/):
  - cb.py, cbsum.py (crossing beams);
  - shd.py, shdsum.py (shadow; exact cached in shd_exact_512.npz);
  - s2.py (eq / nlte / sdif / ddif);
  - rw/run_radwave.py, rw/radwave_disp.py (copied from tests_m1/runs_5f_h2fast, with an `arm=` key), rwtab.py;
  - fig_cb.py, fig_shd.py, fig_1d.py, fig_cloud.py;
  - mkin.py (athinput editor).
- Figures (fig/, white background, < 120 KB each): cb_maps.png, shd_r12g2.png, shd_g8x32.png, nlte.png, diffusion.png, thermal_eq.png, cloud_rho.png.
- Literature sources (arXiv LaTeX): /tmp/claude-28895/.../scratchpad/lit_* (session scratch; not kept).

