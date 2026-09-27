# m1-wedge (09-26): implicit M1 + vet_col on a gravity-bearing spherical-polar rad-hydro wedge

Task: `docs/handover/TASK-2026-09-26-m1-real-wedge.md` (rt-integration 206318c9). Branch `m1-wedge`
from `c1dc8226` (the base of every bitwise gate below). Run tree, scripts, inputs and logs:
`/viper/ptmp2/jinma/sprhd_0926` (RESUME.md there).

## 1. Audit: what the sp M1 path supports when the gas moves (file:line at c1dc8226)

Verified claims of the TASK appendix, with the lines that carry them.

| item | status | where |
| --- | --- | --- |
| radiation force on the gas, per x1 face, half to each cell | yes | `rad_m1_implicit.cpp:8463-8513` (dm1/dm2/dm3), write-back `:8589-8596` |
| energy exchange from the assembled row (e_gas + E conserved) | yes | `rad_m1_implicit.cpp:8448-8458` (`:8455`) |
| gas work W = vbar.dm, removed from E | yes | `rad_m1_implicit.cpp:8514-8531` |
| F in physical (r, theta, phi) components = hydro momentum basis | yes | `rad_m1_sph.cpp:26-31` |
| comoving cell flux F = F0 + (v + v.D) E with the stage velocity (vimp) | yes | `rad_m1_implicit.cpp:8538-8573`; sp rows `ImplicitVimpBuild :5756` |
| enthalpy flux (implicit_enthalpy = plm) with A_f/V_i on sp | yes; not exact on a stretched r grid | `rad_m1_curvilinear_design.md` sect. 8 |
| etotgrav (rho Phi inside u(IEN)) | yes | `rad_m1_opacity.cpp:89`, `rad_m1_implicit.cpp:8446` |
| force_reference = wb_arad (residual momentum, full work) | yes | `rad_m1_coupling.cpp:38-58`, `rad_m1_implicit.cpp:8482` (dmref) |
| a generic point-mass source | **no** | `srcterms.cpp:40-44` has const_accel only; pgens supply gravity |
| sp_cart_polar_momentum / sp_cart_all_momentum with M1 | **no** (never read by rad_m1) | `grep sp_cart src/rad_m1` empty |
| a user callback for the implicit x1 face BC | **no** | `SetImplicitX1BC` (`rad_m1_implicit.cpp:5501`) sets a type + flux only |
| the momentum deposit on sp | plain 1/2-1/2 face split, no area weight: O(dr^2) vs the cell force | `rad_m1_implicit.cpp:8463-8480` |

So the coupling needed nothing; the gap was gravity, an IC, walls and the table hand-over, which
are problem-generator work.

## 2. Route: (b), as a new file

`src/pgen/tests/rad_m1_wedge.cpp`, `<problem>/m1_test = sph_wedge` (pgen `rad_m1_beam`, the
PROBLEM-less binary), dispatched from `RadiationM1Tests2`, declared in `pgen.hpp`. Reasons:
red_giant.cpp (4900 lines) is built around two-stream/conduction/MLT with active defaults (sponge on,
required opac_table/teff/ptop, its own IC march), so an M1 mode would have to neutralise dozens of
switches; the M1 test pgen already owns every M1 hook (force reference, face-flux IC, M1 ghosts, T-S4),
and what it lacked is small and copied from proven code (RedGiantGravity's WB branch,
box_convection's table reader).

What the new setup does (all keys new, `wg_*`; nothing else changes):
- point mass `wg_gm`; gravity source = RedGiantGravity's form (plain `-rho g` + work, or with
  `wellbalance_dynamic + wb_x1` the background's area-weighted pressure drop);
- `wg_ic = grey`: grey (opacity const) atmosphere in hydrostatic + radiative equilibrium, RK4 from
  the top face inward, Eddington closure, `E_top = F_top/(c q)`, ideal gas (test A);
  `wg_ic = file`: `r rho eint [E]` column (any gas-only EOS);
- `wg_phi_eff` (default = force_reference == wb_arad): the x1 WB pair gets
  `Phi_eff = Phi - int a_ref dr`, `a_ref = kappa_t F/c` of the IC, and the M1 coupling applies only
  the residual (`SetForceReference`, face form) -- the spherical analogue of box_convection's
  `wb_phi_eff + force_reference = wb_arad`;
- `<rad_m1>/opacity = table`: `wg_opac_table` + `wg_planck_table` handed over with
  `SetOpacityTables` and the units cross-check (box_convection's code);
- closed x1 walls (`ix1_bc = ox1_bc = user`): ghosts = IC profile scaled by the adjacent active
  cell's ratio (rho, eint), v_r mirrored; M1 ghosts copy (inner) / dark (outer); the implicit
  solve uses its own face BCs (flux in, Marshak out);
- `wg_seed` (deterministic eint seed on global indices), `wg_spot_*` (isochoric hot spot);
- user history: comoving face luminosities (first face, middle, tau = wg_tau_int, Marshak face),
  L_in, E_rad, e_gas and interior (tau >= wg_tau_int) mass, KE, radial KE and momentum, sum p dV.

CUDA-safe: no lambdas inside kernels, host-built columns deep-copied to the device, no host reads of
device Views, no class members inside kernels (all captured by value first).

Results: sections 3-7.

## 3. Gate (A): static grey point-mass atmosphere (hydrostatic + radiative equilibrium, free gas)

Input `tests_m1/runs_6a_m1wedge/inp/gateA.athinput`: GM 1, r 1..1.25 (96 cells), wedge
theta pi/2 +- 0.1, phi 0..0.2 (16x16), ideal gas, kappa_t 1e4 (half absorption), c 300, a 3e5,
F_in 0.009 -> Gamma = kappa F/(c g) = 0.3 everywhere, tau(r_in) 4.4e3, P_rad/P_gas 0.43 deep,
tau = 1 at r = 1.203; vet_col, hesdirk2, mg_gc, WB polytropic + phi_eff + wb_arad.
Interior = tau >= 1 (r <= 1.203). CPU, 4 ranks, t = 6 (~3.5 sound crossings, ~8 diffusion times of
the deep layer), `cpu/A1`:

| arm | Mach_rms interior t=1 / 3 / 6 | L_bot/L_in | L_top/L_in (range) | NON-CONVERGED |
| --- | --- | --- | --- | --- |
| r3 (cfl 0.3) | 1.33e-4 / 1.22e-4 / 4.8e-5 | 1.00039 | 0.9991-1.0002 | 0 |
| h9 (cfl 0.9) | 1.33e-4 / 1.25e-4 / 4.9e-5 | 1.00039 | 0.9991-1.0002 | 0 |
| v4 (vet_col_every 4, lag) | 1.33e-4 / 1.22e-4 / 4.8e-5 | 1.00039 | 0.9991-1.0002 | 0 |
| h9v4 | 1.33e-4 / 1.26e-4 / 4.9e-5 | 1.00039 | 0.9991-1.0002 | 0 |
| nowb (plain -rho g, full force) | 8.0e-5 / 6.2e-5 / 3.7e-5 | 1.00032 | 0.9997-1.0001 | 0 |

- The flow is a decaying radial ringing of the initial truncation mismatch, not a growing mode;
  interior max Mach <= 5e-3 (r 1.19-1.20, tau 1-3), top max 1.3e-2 (`binv.py`).
- Resolution (`cpu/Ares_d2`, 4x4 lateral, t = 3): nx1 48 / 96 / 192 gives Mach_rms 4.0e-4 / 1.2e-4 /
  5.7e-5 and L_bot/L_in - 1 = 1.3e-3 / 3.9e-4 / 1.05e-4: the imbalance is truncation error, ~2nd order.
- The WB pair with Phi_eff is not needed for this smooth profile (nowb is as quiet); it is kept on
  for the He wedge, as in the box.
- dt is set at r_in (radial cell / c_s), never at the top.

## 4. The He wedge (the He-box FeCZ star on a sphere)

`tests_m1/runs_6a_m1wedge/ic/build_ic_sph.py`: the V3edd box star (g0 3.98107e5, F 2.475202e15) with
M = 3.15 Msun -> R0 = sqrt(GM/g0) = 3.2405e10 cm; r = R0 + z over the box range z = -1.8149e8 ..
4.129e7 (84 cells as the box), re-solved as a spherical grey RE + HSE column (Eddington closure, the
box's python EOS port and Rosseland table): tau(r_in) 99.3, Gamma max 0.751 at z = -9.23e7 (box
0.7507), tau = 1 at r = 3.23971e10. The FeCZ is 0.7 % of R here, so the wedge is the box width
(2 deg x 2 deg, 32x32 or 64x64); a 10-20 deg wedge of this star would have 60:1 cells.
Input `inp/he.athinput`: general table EOS gas only, lhllc, etotgrav, WB polytropic + phi_eff,
force_reference wb_arad, opacity table, vet_col, hesdirk2, mg_gc (levels 1), gas Newton + EOS cache.

### 4.1 CPU gates (`cpu/gate_d2`, 4 ranks)
- restart bitwise: He 20 + 20 vs 40 cycles, and grey (A) likewise: last rst payload identical, every
  restarted hst row verbatim in the straight run.
- 2 vs 4 ranks (He, 40 cycles): user hst max relative difference 2.2e-9 (the tiny KE columns),
  <= 3e-12 elsewhere (round-off of the global reductions).
- 300 cycles (t = 48 s): L_bot/L_mid/L_int/L_top all 1.00000-1.00003 of L_in; Mach_rms (tau >= 1)
  2.6e-4 -> 6.9e-4; NON-CONVERGED 0; no NaN/FATAL.

### 4.2 GPU gates (MI300A, apudev, HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1; `gpu/g1_d2`, job 11984593)
- restart bitwise on 2 GPUs (200 + 200 vs 400 cycles): TRUE (33/33 hst rows verbatim).
- 1 vs 2 GPUs, 400 cycles: user hst max rel. difference <= 4e-11, 4.6e-9 in KE, 2.7e-6 in the
  (near-zero) net radial momentum.
- relaxation, 2 GPUs, 84x32x32, grid-scale eint seed 1e-3, 34 500 cycles (t = 5 550 s ~ 1.5 box
  turnovers, 7.5 min): L_top/L_in = 1.000000 +- 2e-6 (800 s means), L_bot 1.000005; dt 0.161 s
  set at r_in (never at the top). The flow is a RADIAL ringing, Mach_rms (tau >= 1) 9.1e-4 ->
  1.4e-3 (800-s means, slowly growing, max 2.1e-3); the lateral Mach decays 9e-6 -> 3.7e-6: the
  grid-scale seed is damped by radiative diffusion and no convection starts in this window
  (bin max |v| 1.4e4 cm/s interior, 3.9e4 at the top).
- `gpu/g4_n2` (job 11984805, apu 1 h, 84x64x64, box-like seed 1e-2 in one lateral mode k = 4 over
  tau 5-50): SUBMITTED, not analysed (see RESUME.md for the analysis commands).

### 4.3 Lever timings (`gpu/g2_d2`, job 11984595; 84x64x64, 2 GPUs, fresh start, cycles 20-120,
two interleaved repeats; ms per cycle, rep1/rep2)

| arm | ms/cycle | wall per sim-second |
| --- | --- | --- |
| ref (mg_gc, team_red, predictor step, one_pass 8, EW 1e-2, chunk auto, every 1, cfl 0.3) | 18.15 / 18.01 | 0.112 |
| implicit_precond = rbgs_fwd | 20.09 / 19.78 (+10 %) | 0.124 |
| vet_col_chunk = 4 | 18.24 / 18.02 (=) | 0.113 |
| implicit_op_team_red = false | 18.02 / 18.14 (=) | 0.112 |
| predictor none + one_pass 0 | 17.20 / 16.83 (-6 %) | 0.106 |
| implicit_one_pass = 0 | 17.93 / 17.87 (-1 %) | 0.111 |
| implicit_lin_ew_max = 0 | 17.99 / 17.98 (=) | 0.112 |
| vet_col_every 4 + time2_vet_col lag | 16.67 / 16.61 (-8 %) | 0.103 |
| cfl 0.9 | 18.52 / 18.33 (+2 %) | 0.038 (2.9x cheaper per sim-s) |

Caveat: measured from the initial state (a near-static column, 1.1 inner BiCGStab iterations per solve, max 5),
where the predictor has nothing to predict; the ranking must be re-measured from an evolved
(convecting) state, e.g. the g4 restart. NON-CONVERGED 0 in all 18 runs.

## 5. Test (C): a 5 % hot spot in the gate-(A) atmosphere (moving gas; `gpu/g3_d2`, job 11984594)

96x32x32, spot at r = 1.1 (tau ~ 1e3), w = 0.03, isochoric T -> T(1 + 0.05 G); t = 0.6 (~2 spot
diffusion times, ~5 sound crossings of the spot). Reference: cfl 0.1. L1 over the domain of the
temperature excess (vs the spot-free twin) and of the momentum field, relative to the reference:

| arm | dT, t = 0.2 / 0.4 / 0.6 | rho v, t = 0.2 / 0.4 / 0.6 |
| --- | --- | --- |
| cfl 0.3 | 0.85 % / 1.6 % / 0.46 % | 1.5 % / 1.3 % / 0.06 % |
| cfl 0.9 | 2.2 % / 3.1 % / 2.7 % | 1.1 % / 1.0 % / 0.32 % |
| cfl 0.3 + every 4 lag | 0.76 % / 1.6 % / 0.24 % | 1.5 % / 1.3 % / 0.06 % |
| cfl 0.9 + every 4 lag | 1.5 % / 1.9 % / 2.0 % | 1.1 % / 1.0 % / 0.32 % |
| cfl 0.3, 1 GPU | identical to cfl 0.3 (2 GPUs) at the printed precision | same |

Spot excess energy (gas + radiation + KE) agrees within 1 % across arms at every time. Verdict on
moving gas: cfl 0.9 costs ~2-3 % in the temperature field against ~1 % at 0.3 (2.9x cheaper per
simulated second); vet_col_every = 4 with lag adds no measurable error (slightly smaller here).
No analytic (Mihalas) damping-rate comparison was made.

## 6. Existing problems stay bitwise (vs the base c1dc8226)
- GPU (`gpu/gb`, jobs 11984680/2/3; hipgate_0926 inputs, 2 GPUs, fresh 12 cycles + restart 6 + 6):
  T-S4 sp wedge (none), He box (box_convection), dhj WASP-121b 1x (deep_hot_jupiter_rt): fresh
  BITWISE (5/5, 15/15, 9/9 files), every restart BITWISE, for base and new (b899f03c).
- CPU: T-S4 (vetevery_0926 wedge, 32x32, 12 cycles) base vs d2 identical (all bin/rst).

## 7. Evolved-state results (09-26 evening; binary n3 = 317d4051 = f8f7232f + open top/sponge)

### 7.1 Convection onset (`gpu/g4_n2`, job 11984805; 84x64x64, k = 4 seed 1e-2 over tau 5-50, t = 27 590 s)
`scripts/conv.py` (horizontal means per dump; v_MLT = 1.86e4 cm/s, the box value at z = -1.02e8):
- velocities grow to v_r' (radial minus horizontal mean) = 1.6-2.1 v_MLT and v_lat = 0.9-2.3 v_MLT
  at t = 24-27.6 ks (0.01 v_MLT at 4 ks), largest in the upper half (z = -1.3e8 .. 0);
- but the convective flux is ZERO: F_conv/F (gas enthalpy + 4/3 E_r fluctuation fluxes) = 1e-9 .. 1e-7,
  sign-alternating in time and depth, and corr(s', v_r') flips sign between dumps. These are waves
  (the growing radial/non-radial acoustic motion of 7.3), not convection. For scale, rho v_MLT^3 is
  only 2e-9 F here: MLT convection in this Fe bump is very inefficient, so an F_conv/F test cannot
  detect it; the sign-coherent s'-v_r' correlation is the discriminant, and it is absent.
- In every 84x32x32 arm (g7/g9) the lateral Mach grows too (e-folding 4.8-8.1 ks from 3e-6 after
  ~15 ks), with the sponge as well; whether this becomes convection needs a longer sponge run.

### 7.2 Lever timings on MOVING gas (`gpu/g5_n3` job 11989483, `gpu/g5b_n3` job 11990268)
Restart of g4 hw.00002.rst (t = 27 590 s, Mach_rms 1.1e-2), 2 GPUs apudev, 180 cycles/arm, cycles
20-180 timed, arms interleaved, 2 (g5) / 3 (g5b) repeats, same binary. `scripts/lev.py`.

| arm | ms/cycle | wall per sim-s | vs ref | Picard/solve | BiCGStab inner/solve |
| --- | --- | --- | --- | --- | --- |
| ref (mg_gc, predictor step, order 2, one_pass 8, EW 1e-2, chunk auto, team_red, every 1, cfl 0.3) | 24.1 | 0.150 | 0 | 2.00 | 7.28 |
| implicit_predictor = none | 30.9-31.2 | 0.193 | +28-29 % | 3.00 | 7.40 |
| implicit_predictor_order = 1 | 24.8 | 0.155 | +3.0 % | 2.01 | 7.81 |
| implicit_one_pass = 0 | 24.1 | 0.150 | -0.1 % | 2.10 | 6.67 |
| implicit_precond = rbgs_fwd | 22.0 | 0.137 | -8.5 / -8.7 % | 2.00 | 8.30 |
| vet_col_chunk = 4 | 24.3 | 0.152 | +1.0 % | 2.00 | 7.28 |
| implicit_op_team_red = false | 24.2 | 0.151 | +0.7 % | 2.00 | 7.28 |
| implicit_lin_ew_max = 0 | 25.4 | 0.159 | +5.6 % | 2.00 | 8.42 |
| cfl 0.9 | 34.9 | 0.073 | -52 % per sim-s (+45 % per cycle) | 2.83 | 10.5 |
| vet_col_every 4 + lag | 26.4 | 0.165 | +9.7 % | 2.00 | 9.99 |

NON-CONVERGED 0 in all 32 runs. The near-static "predictor off -6 %" is REFUTED: on moving gas the
step predictor saves one Picard iteration per solve (3.00 -> 2.00, +28 % without it). From the static
IC the lagged state is already the answer (1.1 inner it/solve), so the predictor's extra work bought
nothing; once the gas moves, the extrapolated increment is a far better first guess than the lagged
state. Also reversed: vet_col_every 4 (+10 %: the lagged tensor costs 2.7 inner it/solve more than
rebuilding it) and mg_gc vs rbgs_fwd (rbgs_fwd -8.7 % here: more iterations, 8.3 vs 7.3, but cheaper
ones; on the T-S4 radiation wedge mg_gc was -18 %, so the sp default is not changed; the He input
may name rbgs_fwd). No default changed.

### 7.3 The radial ringing (`gpu/g7_n3` jobs 11989485/6, `gpu/g7b` 11990284, `gpu/g9_n3/sp` 11990710; 84x32x32, grid seed 1e-3)
- Period 59.2 s in every closed-wall arm (FFT of <v_r>, hst every 2 s). Acoustic estimate from the IC
  (Gamma_1 of gas + radiation): 2 int dr/c_s = 65 s = the closed-closed fundamental of the DOMAIN.
  The star's own fundamental is ~ 2 pi sqrt(R^3/GM) = 1 800 s and cannot live in a 0.7 %-of-R
  wedge. So this is a box mode of the closed walls, not the stellar kappa-mechanism pulsation.
- Growth, amplitude e-folding over t = 0-7.6 ks: cfl 0.3 (w3) 9.9e3 s (0.6 % per period);
  cfl 0.6 (w6) 2.6e3 s (2.3 % per period), i.e. 3.9x faster for 2x dt. A kappa-driven mode grows at a
  dt-independent rate; this growth is numerical anti-damping that scales ~ dt^2, the same sign of dt
  dependence as the He-box KE growth (kedt_0926; no operator hunt here). w3 saturates at
  |<v_r>|/c_s = 1.8e-3 after ~15 ks; w6 at 2.8e-3 had stalled to ~0.7 s/cycle after 20 000 cycles
  (1-GPU side-by-side launch, see below).
- open top (wg_bc_top = open, cfl 0.3, 36 ks): no 59-s peak, |<v_r>|/c_s decays 2.6e-4 -> 1.0e-4;
  but the column settles: <v_r>/c_s = -1e-4 steady, interior mass +2.8 %, E_rad +4.3 %. Not usable.
- top sponge (wg_sponge_rate = 0.02 /s above r(tau = 1), w = ((r - r0)/(r_top - r0))^2; cfl 0.3,
  25.7 ks): the mode decays, |<v_r>|/c_s 1.65e-4 -> 3.7e-5 and stays; Mach_r (tau >= 1) flat 9.1e-4
  (wall: 6.4e-3); L_top/L_in 1.000000, M_int drift 9e-4, dt 0.161 s unchanged (still set at r_in);
  cost +0.7-1.8 % per cycle (g8 job 11990669, 400-cycle A/B on 2 GPUs and on each GPU).
  Now the DEFAULT (7.6).
- The g7b sponge arm first ran at 0.57 s/cycle, and w6 slowed likewise after 20 000 cycles; both were
  the 1-rank GPU-1 half of a side-by-side launch on an apu node. g8 measured the same arms at
  12.4 ms/cycle on either GPU with identical iteration counts, so it is the launch, not the code:
  use 2-rank runs.

### 7.4 Convergence order (coordinator task 5; `gpu/g6_n3`, job 11989484; `scripts/order.py`)
Gate-(A) grey atmosphere + a smooth hot spot (amplitude 1e-2, w = 0.04 at r = 1.1), production sp
settings (vet_col, hesdirk2, time2_vet_col predict, vimp default), t = 0.2 (~0.5 spot diffusion times).
SPACE, cfl 0.1, levels 48x16x16 / 96x32x32 / 192x64x64, each vs the next finer level restricted
(2x2x2 volume average), relative L1 of the spot perturbation (spot minus spot-free twin):

| var | e(48-96) | e(96-192) | order |
| --- | --- | --- | --- |
| T | 8.2e-3 | 2.4e-3 | 1.80 |
| rho | 1.9e-2 | 5.9e-3 | 1.69 |
| E | 1.1e-2 | 2.7e-3 | 2.07 |
| v | 1.8e-2 | 4.3e-3 | 2.11 |

(Raw fields: T 1.94, rho 1.99, E 1.97; the raw v of the spot-free twin is noise, so its order (1.09) is
not meaningful.)
TIME, 96x32x32, vs cfl 0.025, relative L1 of the perturbation, local order p between neighbours:

| var | every 1: cfl 0.8 / 0.4 / 0.2 / 0.1 | p (0.4-0.2-0.1) | every 4 + lag: 0.8 / 0.4 / 0.2 / 0.1 | p |
| --- | --- | --- | --- | --- |
| T | 7.0e-3 / 8.1e-3 / 3.4e-3 / 1.4e-3 | 1.26, 1.30 | 7.8e-3 / 7.2e-3 / 2.8e-3 / 1.1e-3 | 1.37, 1.34 |
| rho | 4.1e-3 / 1.2e-3 / 4.0e-4 / 1.4e-4 | 1.59, 1.52 | 4.2e-3 / 1.3e-3 / 4.5e-4 / 1.7e-4 | 1.54, 1.44 |
| E | 1.9e-3 / 1.5e-3 / 4.8e-4 / 1.6e-4 | 1.67, 1.62 | 2.0e-3 / 1.5e-3 / 4.9e-4 / 1.6e-4 | 1.65, 1.58 |
| v | 8.4e-3 / 3.6e-3 / 1.6e-3 / 6.5e-4 | 1.22, 1.27 | 1.6e-2 / 7.2e-3 / 3.3e-3 / 1.5e-3 | 1.13, 1.12 |

Stiffness: c rho kappa_P dt at the spot = 2.0e3 (cfl 0.1) .. 1.6e4 (cfl 0.8) (domain 0.06 .. 1.1e5);
diffusion number over the spot width c dt/(3 rho kappa_t w^2) = 0.006 .. 0.05. Verdict: SPACE second
order (1.7-2.1). TIME between first and second order, 1.2-1.6, in the strongly stiff regime (c rho kappa
dt >> 1), consistent with the H-ESDIRK2 stage-order-1 reduction of the radiation-only sp tests
(1.4 at c rho kappa dt = 16); cfl 0.8 is not yet in the asymptotic range (T error non-monotone).
vet_col_every 4 + lag: same T/rho/E errors, v error 2.0-2.3x larger, order ~1.1 (first order in the
lagged tensor, as expected).

### 7.5 CUDA read-through of rad_m1_wedge.cpp and restarts
- Read against the TASK rules: every kernel captures only locals (Views copied from namespace globals
  or members: crho = wg_rho_, uh = ph->u0, gm = wg_gm_, the new sponge/open-top scalars); the two
  device helpers are KOKKOS_INLINE_FUNCTION free functions; no generic or nested lambdas; host reads
  go through create_mirror_view + deep_copy (tables, column, T/kappa); mb_bcs is read through its
  d_view (the mesh keeps it synced); the namespace-scope Views are released in pgen_final_func. No
  finding. Not compiled with nvcc here (no CUDA on viper).
- Time2RstSet: the wedge (hesdirk2 + predictor) writes and reads the slope/predictor channels. The
  DeepCopyAcross fix came with the merge of f8f7232f (m1-port 48321155); GPU restart of the EVOLVED
  g4 state (g5 rs_str 250 cycles vs rs_rr 100 + 150): 3/3 rst payloads identical, hst rows verbatim.

### 7.6 Follow-ups (09-26 night)
**rbgs_fwd as the sp default: NOT flipped.** A confirmation from the evolved SPONGE state (g9 hw.00004.rst,
84x32x32, t = 25.7 ks, 400 cycles, 3 interleaved repeats on 2 GPUs + a 1-GPU pair; `gpu/g10_n3`, job
11991033) gives ms/cycle 9.25 (mg_gc) vs 8.92 (rbgs_fwd), -4.0 % per sim-s on 2 GPUs (rep spread ~1 %);
1 GPU 8.96 vs 8.58, -4.3 %. Inner iterations per solve 1.16 (mg_gc) vs 1.20, Picard 1.00. The gain is
under the 5 % bar, so the flip (commit badf65fd, gates not run) is parked on local branch
m1-wedge-sp-rbgs. The -8.7 % of 7.2 was measured on the faster-moving closed-wall state (Mach 1e-2).

**Top sponge ON by default in m1_test = sph_wedge** (736c49bf): wg_sponge_rate = 0.02 above
wg_sponge_r0 = the IC's tau = wg_tau_int radius; wg_sponge_rate = 0 restores the old scheme.
Gates (`gpu/g13` job 11991209, `gpu/g15` job 11991353, `gpu/gb2` jobs 11991211-3/11991355, vs base
b5 = 5d6a68df): new default == old + explicit 0.02, and new + wg_sponge_rate 0 (+ wg_wall_zero_flux false)
== old default: all bin/rst files and the old hst columns identical; restart of the new default bitwise;
He box, dhj and T-S4 (sph_atm) fresh + restart bitwise (the T-S4 hydro.hst differs only by rows appended
when b5 was run twice).

**Zero-flux walls** (0d4a27f8, wg_wall_zero_flux, default true). The closed x1 walls were LEAKY: the
ghosts carry the scaled initial profile, not the mirror image of the edge cell, so the Riemann flux through
the wall face carried mass (L x 1: +1.0e-4 of the mass in 1.4 ks; the super-Eddington test below: -17 %
in 2.5 ks). The earlier open-top +2.8 % in 36 ks was INFLOW through the top (the copied edge velocity
pointed inward). Fix: after the update, the wall faces' mass, transverse-momentum and energy fluxes are
taken back out of the edge cells (the pressure flux stays); wg_bc_top = open now mirrors an inward edge
velocity and keeps only outward mass flux. Result (84x16x16, `gpu/g15`): wall dM/M0 = -3e-13 in 1.4 ks,
open top -1.1e-5 once (initial settling) then constant; no inflow. hst gains M_tot, Mdot_top, Mdot_bot
(Riemann face fluxes of the last stage, before the correction).

**Open top for a super-Eddington envelope** (`gpu/g16_n6`, job 11991354): He wedge 84x16x16, F_in x 1.5
(Gamma_max of the IC 0.75 -> 1.13 in the Fe bump), 2.8 min per arm, 2 GPUs:

| arm | t reached | dM/M0 (0.65 / 1.3 / 2.0 / 2.6 ks) | dt min [s] | Mach (tau >= 1) end | L_top/L_in end | ms/cycle |
| --- | --- | --- | --- | --- | --- | --- |
| open, no sponge | 2 618 s | -0.38 / -0.41 / -0.43 / -0.44 | 0.1135 | 4.6e-3 | 1.0000 | 10.3 |
| open + sponge | 2 719 s | -0.38 / -0.41 / -0.42 / -0.43 | 0.1170 | 4.6e-3 | 1.0000 | 9.9 |
| wall + sponge | 2 502 s | 2e-13 (closed) | 0.1006 | 3.5e-2 | 0.9975 | 10.0 |
| wall | 2 351 s | 1e-13 (closed) | 0.0995 | 4.9e-2 | 1.0031 | 10.6 |

NON-CONVERGED 0 in all arms. With the open top the super-Eddington layer is expelled cleanly: 38 % of
the mass in the first 650 s, then a declining wind (Mdot M0/Mdot ~ 2e4 -> 3e4 s in the earlier leaky run;
-3 % of M0 over the last 2 ks here), the integrated top face flux equals the mass lost, no reflected
pulse (interior Mach falls to 4.6e-3), dt is not limited at the top (min 0.114 s, set at r_in, vs 0.100 s
closed), and the Marshak M1 face holds L_top/L_in = 1.0000 with outflowing gas. Closed or sponged, the
expelled gas piles up under the lid (top density x 24 in the leaky run) and the interior stays 8-10x more
agitated. Cost: same per cycle (10.3 vs 10.0-10.6 ms), cheaper per sim-s (larger dt). Not yet a steady
wind in 2.6 ks; the top-density floor and a longer run decide wind vs inflated envelope.

### 7.7 Super-Eddington wedge to steady state (`gpu/g17_n7`, jobs 11991958 + 11991959, 2 x 4 h apu)
He wedge 84x16x16, F_in x 1.5 (IC Gamma > 1 at z = -1.27e8 .. -7.3e7 cm, tau 18-55, below the IC
photosphere at z = -8e6), open (outflow-only) top, no sponge, 1 GPU per arm on one node (per-step
--gres=gpu:1 --mem=90G; an srun step without them waited for the GPUs of the first step).
New option wg_bc_bot = reservoir (2c1825e1, default wall; bitwise when off: `gpu/g18b` job 11991989,
12/12 files + hst vs b7 = 041fac8f, wall and open top): the inner ghosts hold the IC rho, eint with the
edge cell's radial velocity (floating). `ana/se.py <run> <window>`: windowed Mdot, mass, L_top/L_in,
horizontal-mean profiles and their drift.
Time scales of the super-Eddington layer (IC): thickness 5.4e7 cm, sound crossing 9 s, thermal
(energy content / F) ~1 s; flow time = domain transit 2.2e8 cm / v_r: W2 ~150 s (v_r 0.6-3e6 cm/s),
W1 ~4e4 s (v_r ~5e3 cm/s).

| | W1 closed bottom (t = 545 ks) | W2 reservoir bottom (t = 353 ks) |
| --- | --- | --- |
| Mdot_top (wedge, 1.16e-3 sr) | 2.6e13 -> 3.5e12 g/s over 50-500 ks, ~ t^-0.9, still -11 % per 50 ks | 2.760e17 +- 0.1 % (50-ks windows, 0-350 ks) |
| mass | -57 % by 50 ks, -67 % at 545 ks, still -0.4 % per 50 ks | +12.2 % (fills, then flat to 1e-4); Mdot_bot = Mdot_top to 5e-5 |
| L_top/L_in | 0.99999 | 0.99539 (0.46 % of L lifts the gas) |
| profile drift per 50 ks | rho 1.2 %, T 0.5 %, v_r 1 %, Gamma 3 % | rho 2e-4, T 7e-5, v_r 1e-3, Gamma 2e-4 |
| photosphere | sinks z = -8e6 -> -7.1e7 (deflation) | z = +3.5e7 (inflated, near the top) |
| Gamma below the photosphere | < 0.87 everywhere | 0.57-0.90 everywhere; no Gamma > 1 layer left |
| max Mach_r / dt min / NON-CONVERGED | 0.03 / 0.1135 s / 0 | 0.62 / 0.1133 s / 0 |

Verdicts. W2 reaches a steady outflow within ~10 ks (criteria met: Mdot constant to 0.1 % over > 2000
flow times, Mdot_bot = Mdot_top, L_top steady, profiles frozen), but it is NOT a radiatively driven wind:
after the adjustment Gamma < 1 everywhere, the flow is subsonic throughout (Mach 0.62 at the top,
v = 3.2e6 << v_esc = 1.6e8 cm/s), and it is pushed by the reservoir pressure through the open lid.
Its rate, 2.39e20 g/s/sr (4.8e-5 Msun/yr over 4 pi), is 0.80 of the photon-tiring limit L/(GM/R)
(3.0e20 g/s/sr) and 90x a continuum-driven estimate Mdot ~ (Gamma - 1) L/(c c_s) = 2.6e18 g/s/sr, i.e.
set by the two boundaries, not by the physics. The flow starts at the base (tau 85), below the
photosphere. W1 has no steady state in finite time: the super-Eddington layer is expelled in the
first ~2 ks (the launch is below the photosphere, tau 18-55), then the deflated envelope (Gamma < 0.9)
keeps evaporating through the open top with Mdot ~ t^-0.9. The Gamma > 1 flagged in the top 2-4 cells
(tau ~ 0) is the free-streaming E gradient, not a force balance. Nothing numerical limits either arm:
no floor reached (rho_min 1.2e-10 / 7e-8 vs dfloor 1.6e-12), dt set at r_in (0.113-0.18 s), no
NON-CONVERGED. A physical continuum-driven wind needs the domain through the sonic point (and to
several R for escape), with a base that supplies mass at the IC pressure: that is the next test.

## 8. What is next for a gravity-bearing He wedge
1. Convection onset with the top sponge on (the g4 velocities were the closed-wall mode's waves):
   a long 84x64x64 k = 4 run with wg_sponge_rate = 0.02, judged by the sign of corr(s', v_r').
2. The dt-dependent anti-damping of the acoustic mode (7.3) goes to the kedt_0926 operator hunt.
4. The extended 4 Msun presupernova star (FeCZ 0.64-0.97 R, where a 10-20 deg wedge is natural)
   needs an MLT-based column: its radiative column is super-Eddington in the Fe bump.
5. Caltech H200 build of rad_m1_wedge.cpp (written to the CUDA-safe rules; not compiled with nvcc).
