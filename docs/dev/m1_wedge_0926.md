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

## 7. What is next for a gravity-bearing He wedge
1. Convection onset: g4 (box-like k = 4 seed) -> KE, F_conv/F, L_top vs L_in against the box.
2. Re-time the levers from the evolved g4 state (the IC timings above are not representative).
3. The radial ringing grows slowly (9e-4 -> 1.4e-3 in 6 000 s): decide whether it is the physical
   Fe-bump kappa-mechanism or the closed walls (an open/damped outer wall, as box bc_mode 3 / red
   giant open BCs, is the test).
4. The extended 4 Msun presupernova star (FeCZ 0.64-0.97 R, where a 10-20 deg wedge is natural)
   needs an MLT-based column: its radiative column is super-Eddington in the Fe bump.
5. Caltech H200 build of rad_m1_wedge.cpp (written to the CUDA-safe rules; not compiled with nvcc).
