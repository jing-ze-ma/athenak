# runs_5q_sporder2b: second order in space AND time for implicit M1 + VET on the sp wedge

- **Date.** 2026-09-25, viper. **Branch** `m1-sp-order2b` from rt-integration fe12a518
  (rt-integration 0c44f3cb, m1-coarse2, merged in afterwards), worktree
  `/viper/ptmp2/jinma/wt_sporder2b`. **Work dir** `/viper/ptmp2/jinma/sporder2b_0925`
  (`bin/`, `st/` study runs, `gate/`, `gpu*/`, `gates/`, logs).
- **Binaries.** `ref` = git archive fe12a518; `new` = snapshots of this branch
  (`scripts/snap.sh`, `scripts/build.sh`; gcc 14 + openmpi 5, Release, MPI; GPU rocm 6.3,
  gfx942, `HSA_XNACK=1`, `HSA_NO_SCRATCH_RECLAIM=1`).
- **Driver.** `st.py` (cases, arms, parallel runs, per-variable self-convergence),
  `study_final.txt` / `study_gas.txt` (the (case, arm) lists), `scripts/study.sh`,
  `scripts/evalall.sh`. Results: `RESULTS_radial.txt`, `RESULTS_gas.txt`,
  `RESULTS_lateral.txt`.

## 1. New sp defaults (keys; a restart whose file lacks a key keeps the old value)

All are GetOrAdd on the spherical-polar wedge with the default `!restart_run` (the
vet_col_surface_q / implicit_vimp convention; the resolved value is echoed and written
to the next restart). Cartesian meshes read them only when named: unchanged.

| key | new sp default | old | what |
| --- | --- | --- | --- |
| `implicit_marshak_face` | `linear` with eddington / vet_sc / tau / vet_col | `cell` | runs_5o. `cell` stays with the lagged m1 closures: the runs_5h `pp_np` atmosphere (closure m1, closure_lag step, implicit_cfl 1e6) went from 30 round-off NON-CONVERGED solves to 54 with Picard blow-ups (resid 5.9e28) |
| `vet_col_order2` | `true` | `false` | runs_5o |
| `vet_col_fk_min` | `0` | `1/3` | the f_K clamp (sect. 3) |
| `vet_col_reflect_top` | `true` | `false` | mirrored top under a reflecting outer x1 (sect. 4) |
| `vet_col_surface_face` | `true` | `false` | Marshak q = H(face)/J_face with the linear face E (sect. 3) |
| `time2_vet_col` | `predict` | `lag` | the vet_col tensor of the hesdirk2 stages (sect. 2) |
| `time2_vstage` | `true` | `false` | implicit_vimp: stage velocity in the lagged v-terms, gas work in the E row (sect. 5) |
| `implicit_precond` | `mg_gc` (+ `implicit_mg_levels = 1`) where `rbgs_fwd` was the default | `rbgs_fwd` | runs_5p_coarse2 |

## 2. Time: the vet_col tensor at t^{n+1}

Both hesdirk2 implicit stages sit at c = 1 (the explicit part is Heun), so D(U^n) is an
O(dt) lag: the step is first order in time with vet_col (pulse, time-only refinement at
nx1 = 128, `pulse_u_T`: lag 1.09 1.06 1.03). Arms (`time2_vet_col`):

- `predict`: ONE build per step at P = U^n + dt K1 (radiation), Heun predictor +
  dt K1 (gas; the hydro u0 at that point is already the old vector, so g dt K1 is added),
  inner face flux F^n + dt K1_F. K1 is restart state: restarts stay bitwise.
- `rebuild`: predict for stage 1 plus a second build at the stage-1 solution Y1 (E, T
  saved at the end of the stage-1 solve).
- `extrap`: D^n + (dt/dt_prev)(D^n - D^{n-1}) (DIAGNOSTIC; history not in the restart).

| pulse_u_T (time only) | L1 (levels 1..16) | orders |
| --- | --- | --- |
| lag | 5.56e-2 2.61e-2 1.25e-2 6.12e-3 | 1.09 1.06 1.03 |
| predict | 2.08e-2 4.98e-3 1.20e-3 2.93e-4 | 2.06 2.05 2.03 |
| rebuild | 1.83e-2 4.37e-3 1.06e-3 2.62e-4 | 2.07 2.04 2.02 |
| extrap | 2.64e-2 6.57e-3 1.61e-3 3.95e-4 | 2.01 2.03 2.02 |

Rebuild is 12 % more accurate but costs a second build (+18 % per cycle on the He wedge
GPU, sect. 7); extrap is not restart-safe. **predict is the default.**

## 3. The f_K clamp at 1/3, and the surface q

- **Why it existed.** Copied from vet_sc, whose projection takes chi = the largest
  eigenvalue of D (>= 1/3 by construction) or n.D.n about the flux. vet_col takes
  K_rr/J about r_hat, which is physically in [0, 1]: it drops below 1/3 wherever the
  intensity is limb-brightened (behind an outgoing radiating front, inside an emitting
  shell); D = diag(f, (1-f)/2, (1-f)/2) stays realizable for 0 <= f <= 1. The clamp
  put a kink into f_K where K/J crosses 1/3.
- **Fix.** `vet_col_fk_min` (0 on sp; the flux-axis option keeps 1/3). NON-CONVERGED 0
  everywhere.
- **Surface q.** With the linear Marshak face E, q = H(face)/J(top cell) mixes a face
  and a cell quantity. `vet_col_surface_face`: J extrapolated to the face the same way
  (r^2 J linear from the top two shells, limited), so F = c H(face) (E/J)_face.
  pulse_u (space-time) at 256 -> 512: cell q 1.90 / 1.84 (E / F), face q 2.08 / 2.07.
- **The test itself.** The 5o pulse used e_out = 1e-2 with a Marshak top: F = 0 at the
  top is incompatible with c q E, and that start-up transient runs inward at c sqrt(f)
  as a grid-scale front (Eddington too: 2.00 1.99 1.41 at e_out 1e-2, 2.01 1.99 2.00 at
  1e-8). All Marshak-top cases here use e_out = 1e-8.
- **Rays.** vet_col_ncore = 64 in the study (base_rad). ncore 8 caps pulse_u at 512
  (1.90 1.97 1.14 vs 1.93 2.06 2.08).

## 4. Reflecting outer x1 under vet_col

Specular reflection at the top sphere keeps the impact parameter p, so per ray
I_in(top) = b + a I_in: b the outgoing intensity at the top face of a sweep from a
vacuum top, a the round-trip transmission (product of the segment factors; 0 for core
rays that start from the diffusion intensity). `vet_col_reflect_top`: a first sweep gives
b and a, a second sweep runs from I_in = b/(1 - min(a, vet_col_reflect_amax = 0.99)).
Exact for the steady transfer problem, both kernels (team and column), 2x the build cost
when (and only when) the outer x1 is reflecting. Limitation: in a nearly transparent
closed cavity (a -> 1) the intensity has memory the instantaneous formal solution cannot
have; the cap bounds I_in there.

| case (space-time) | vacuum top (old) | mirrored top |
| --- | --- | --- |
| lat_th_mid (lateral theta mode, rho kappa 1, reflecting r), E | 0.62 0.71 | 2.01 2.00 |
| same, F1 | 0.61 0.66 | 2.00 1.88 |
| rtop_u (radial pulse, reflecting r), E | 1.95 2.06 2.08 | 1.96 2.06 2.08 |

## 5. implicit_vimp: stage-time consistency (user addition)

Every convergence case runs with the production defaults (hesdirk2 + vimp where valid:
gas_feedback and dbg_gas_force on) and a vimp-off control. Two O(dt) lags were found
and fixed (`time2_vstage`, sp default):

1. **The derived cell flux** F = F0 + a E took a = v + v.D of the STAGE-START velocity;
   the stage's radiative kick is missing (O(dt)). Now a at the velocity the write-back
   gives the gas (and E0 - E at v_old + dv^k inside the Picard loop). rsw_u_T F order
   1.00 -> 1.99.
2. **The gas work** v.dm of the radiative kick was subtracted from E AFTER the solve, a
   splitting error in every stage. Now it is in the E row (its own kernel,
   `ImplicitWorkRow`, lagged at the pass's dv^k; the write-back removes only
   work - w^k). mvrs_u_T (passive E in a moving scattering gas) E 1.25 1.14 1.08 ->
   2.00 2.00 2.01.

The tensor inside the vimp rows (a = v + v.D) is tau_ten: time2_vet_col = predict puts
it at t^{n+1} as well. Without vimp the radiation-dominated moving slab is first order
in every variable (rsw_u_T novimp 1.12 1.06 1.03).

## 6. Combined space-time convergence (dx and dt refined together)

L1_V(q_n - R q_2n)/L1_V(|q_n - <q_n>|); radial n = 32..512 (4 x 4 angular), lateral
s = 1..8 (n1 = n_lat = 16 s); `_T` = time only (nx1 = 128, dt halved, 5 levels).
`old` = the new binary with every new key at its old value (bitwise the fe12a518
behaviour, sect. 8). Last three orders:

| test | var | old | new |
| --- | --- | --- | --- |
| vet_col radiating pulse, rho kappa 1, Marshak top, uniform | E / F | 1.21 0.81 0.87 / 1.28 0.93 0.86 | 1.93 2.06 2.08 / 1.87 2.05 2.07 |
| same, stretched r | E / F | 1.51 1.23 1.00 / 1.52 1.25 0.98 | 1.61 1.98 2.07 / 1.53 1.97 2.08 |
| free-streaming shell (rho kappa 1e-3), reflecting, uniform | E / F | 0.87 0.89 0.87 / 0.94 0.91 0.88 | 1.95 2.09 2.08 / 2.00 2.07 2.07 |
| same, stretched | E / F | 1.08 0.96 0.93 / 1.09 1.03 0.95 | 1.62 2.03 2.08 / 1.88 2.04 2.09 |
| free-streaming shell leaving through a Marshak top | E / F | 1.01 0.96 0.98 / -0.12 0.38 0.66 | 1.91 2.07 2.08 / 2.08 2.05 2.08 |
| pulse between reflecting walls | E / F | 1.19 0.86 0.75 / 1.25 0.89 0.74 | 1.96 2.06 2.08 / 1.91 2.04 2.07 |
| vet_col atmosphere transient (T-S4, c = 100), uniform | E / F | 1.91 1.80 1.79 / 1.69 1.49 1.32 | 1.99 1.99 1.99 / 1.97 1.98 1.99 |
| same, stretched | E / F | 1.96 1.81 1.87 / 1.76 1.59 1.36 | 2.00 1.98 2.00 / 1.98 1.97 2.01 |
| lateral theta mode on the atmosphere (s = 1..8) | E / F1 / F2 | 1.92 1.91 / 1.76 1.67 / 1.58 1.58 | 1.98 1.99 / 1.96 1.96 / 1.83 1.96 |
| lateral phi mode on the atmosphere | E / F1 / F3 | 1.92 1.91 / 1.86 1.69 / 1.72 1.64 | 1.98 1.99 / 1.99 1.97 / 1.71 1.94 |
| lateral theta mode, reflecting shell | E / F1 | -0.14 0.75 / -0.30 0.62 | 2.01 2.00 / 2.00 1.88 |
| stiff coupled, rho kappa_P = 1 (vet_col) | E / T / Eg | 1.81 1.65 1.37 / 1.89 1.95 1.91 / 1.85 2.12 2.22 | 2.04 2.01 1.97 / 2.02 1.89 2.12 / 1.89 1.93 2.20 |
| stiff coupled, rho kappa_P = 100 .. 1e6 | E, F, T, gas | 1.5-2.5 (pre-asymptotic, same as new) | same (identical to 3 digits) |
| moving gas v_r, gas carries E (vimp on) | E / F / T / M | 2.08 1.96 1.32 / 2.36 2.24 1.95 / 2.0 / 2.2 | 1.98 1.99 2.05 / 2.47 2.32 1.88 / 1.99 1.98 2.01 / 2.40 2.33 2.16 |
| radiation-dominated moving slab (vimp on) | E / F / T / M / Eg | 1.96 2.02 2.02 / 1.90 1.91 1.85 / 1.96 2.02 2.02 / 1.93 1.98 1.99 / 1.77 1.88 1.96 | 1.96 2.02 2.01 / 1.94 1.98 1.99 / 1.96 2.01 2.01 / 1.93 1.98 1.98 / 1.77 1.88 1.95 |
| same, vimp OFF (control) | E / F / T | | 2.36 2.28 1.54 / 2.22 2.15 1.52 / 2.33 2.30 1.55 |
| lateral flow v_theta + theta mode (vimp on / off) | E / F2 / T / M2 | 1.73 1.92 / 2.26 2.19 / 1.73 1.90 / 2.38 2.28 | 1.72 1.91 / 2.25 2.22 / 1.73 1.90 / 2.38 2.28 (off: identical) |

Time only (`_T`, nx1 = 128):

| test | var | old | new |
| --- | --- | --- | --- |
| vet_col pulse | E / F | 1.05 0.90 1.12 / 1.12 1.09 1.17 | 2.06 2.05 2.03 / 2.04 2.04 2.03 |
| free-streaming shell, Marshak top | E / F | (lag) 1.07 1.04 1.03 | 2.00 2.03 2.02 / 2.02 2.03 2.02 |
| atmosphere transient | E / F | 1.09 1.08 1.08 | 1.41 1.48 1.55 (1.57 1.74 at levels 32, 64; errors 1e-7 .. 5e-10, the Picard floor) |
| stiff, rho kappa_P = 1 | E / T / M | 1.78 1.57 1.49 / 1.72 1.57 1.40 / 1.65 1.38 1.22 | 1.96 1.99 2.00 / 1.99 2.00 2.01 / 1.91 1.96 1.98 |
| stiff, rho kappa_P = 100 | E / F / T | 2.10 2.16 2.17 / 1.92 1.86 1.70 / 2.05 2.09 2.09 | 2.09 2.13 2.15 / 2.03 2.08 2.11 / 2.04 2.08 2.10 |
| stiff, rho kappa_P = 1e4 and 1e6 | all | 2.00-2.01 | 2.00-2.01 (be: 0.9-1.1) |
| radiation-dominated moving slab | E / F / T | 1.96 1.90 1.80 / 1.01 1.00 1.00 / 1.92 1.89 1.83 | 2.01 2.00 1.98 / 2.01 2.00 1.99 / 1.95 1.97 1.97 |
| same, vimp off | all | | 1.03-1.14 |
| moving gas carrying E | E / F | 0.81 1.00 0.94 / 2.55 0.95 2.15 | 1.89 1.98 1.98 / 1.94 1.99 1.98 |
| passive E in moving scattering gas | E / F | (vstage off) 1.25 1.14 1.08 / 1.31 1.20 0.90 | 2.00 2.00 2.01 / 2.01 2.01 2.02 |

(The stiff kappa_P = 1 time-only row first showed 1.40 / 1.01: the predict state added
dt K1 to a gas state that already carried (1-g) dt K1, fixed.)

## 7. Stiff order reduction

Stiff cases: grey absorption rho kappa_P = K, gas at radiative equilibrium with the
pulse (gas_teq), arad = 1e8, rho = 270 (heat capacities comparable), c dt rho kappa =
cfl dx K from 0.02 (K = 1) to 1.6e4 (K = 1e6, coarsest level). **No order reduction**
in E, F, T or the gas in any case: 2.00-2.15 in time for K = 100 .. 1e6, 1.96-2.01 for
K = 1. H-ESDIRK2 has stage order 1 (runs_5k_h2acc: Prothero-Robinson reduction), but
our stiff mode (E - a T^4 relaxation) has a slow manifold that moves with the diffusion
and the error it leaves is O(dt^2 / (K c dt)): the stiffer, the smaller. Where it would
appear: a FORCED stiff mode whose equilibrium moves on the step time scale (a boundary
or source that drives E - a T^4 at rate ~1/dt) - none seen in our regimes. A fix would be
a stage-order-2 ESDIRK (e.g. ESDIRK3 with an explicit first stage, 3 implicit solves);
not implemented.

## 8. Gates (`scripts/gate_cpu.sh`, `gate_final.log`; `scripts/gate_gpu.sh`)

- **Cartesian** (A: be named 10 cases, B: key absent 7 cases incl. Milne vet_col): CPU
  17/17 BITWISE vs fe12a518; GPU slab / box3d_nd be / box3d_nd hesdirk2+vimp 3/3
  BITWISE (a first build with the work term inside the assembly kernel changed the GPU
  code generation of the Cartesian kernels, hst tot-E 1e-14 at cycle 1; the term moved
  to its own kernel, bitwise again).
- **Old sp behaviour named** (CO: be named 14 cases, FO: key absent 7 cases, every new
  key at its old value, implicit_precond = rbgs_fwd where mg_gc would resolve, found by
  a 1-cycle probe, `scripts/pcfix.py`): 21/21 BITWISE (`gate_coFo.log`).
- **Old restart** (O: sph_sym_gas_rst restarts written by fe12a518 - vet_col, vet_col +
  reflecting top, Eddington - continued by ref and new; O2: sp_sph_atm_vc, the
  implicit_precond default): 4/4 BITWISE; the keys come back at their old values.
- **Restart with the new defaults** (R: the three cases, R2: mg_gc): straight vs restart
  at the midpoint, rst payloads 4/4 BITWISE.
- **NON-CONVERGED**: 0 in every new-default run; the pre-existing round-off-floor ones
  (resid ~1e-12 at tol 1e-12: A_pp_np 30, rw 327, sp_sph_fs 60/119, sym_gas 26) are
  identical in ref and new (C_symg with the new Marshak face: 27 vs 26, all at resid
  1.8e-12).
- C / F with the keys absent differ from ref exactly where a Marshak face, vet_col or
  mg_gc is involved (sd, sds, sym, symg, mk, vatm, vsym, vsymg, vstr); the rest
  BITWISE. O2/R2 rst headers: the old restart continues with implicit_precond = rbgs_fwd,
  vet_col_order2 = 0, marshak cell, time2_vet_col lag; a new run writes mg_gc / 1 /
  linear / predict.
- `tests_m1/gates/gates.py`: GATES: PASS. `tst/test_suite/rad_m1`: 3 passed (unchanged
  expectations: those tests are Cartesian).
- cpplint (repo filters) clean on the changed lines; flake8 (90 columns) clean.

## 9. GPU (apudev job 11978848, 2 x gfx942, `gpu/log.out.11978848`)

- Bitwise ref vs new: slab (1 rank, be), box3d_nd be (2 ranks), box3d_nd key absent
  (hesdirk2 + vimp), He wedge with the old keys named: 4/4 BITWISE.
- Correctness: T-S4 vet_col wedge n = 64, 100 steps, new defaults, 1 GPU vs 1 CPU rank:
  max dE/E 7.5e-10, dF 4.6e-10 (vacuum top); reflecting top 3.3e-10 / 5.4e-11 (stage
  linear tolerance 1e-9); identical Picard statistics, NON-CONVERGED 0.
- Cost, He wedge 96 x 128 x 128 (hesdirk2, vet_col, 16 blocks, 2 GPUs, 40 steps, same
  binary, interleaved, 2 repeats):

  | arm | wall (s) | ms/cycle | vet_col build | inner its/solve |
  | --- | --- | --- | --- | --- |
  | old keys (fe12a518 behaviour, rbgs_fwd) | 1.243 / 1.246 | 31.1 | 3.77 ms | 10.60 |
  | new defaults but implicit_precond = rbgs_fwd | 1.306 / 1.317 | 32.8 (+5.4 %) | 4.78 ms | 11.19 |
  | **new defaults (mg_gc, levels 1)** | 1.035 / 1.045 | **26.0 (-16 %)** | 4.80 ms | 1.97 |
  | new defaults + time2_vet_col = rebuild | 1.219 / 1.209 | 30.4 | 4.78 ms (x2 per step) | 1.74 |

  NON-CONVERGED 0 in all eight. The order-2 work costs +5 % (vet_col_order2's build,
  +0.6 inner its/solve); predict costs nothing over lag; mg_gc more than pays for it.

## Files

`st.py`, `solib.py` (from runs_5o), `study_final.txt`, `study_gas.txt`, `RESULTS_*.txt`,
`scripts/`: `snap.sh`, `build.sh`, `study.sh`, `evalall.sh`, `gate_cpu.sh`,
`gate_gpu.sh` (+ `gate_gpu2.sh`, `gate_gpu3.sh`: the tight-tolerance check and the
bisect of the GPU code-generation difference), `gcmp.py`, `rstcmp.py`, `addkeys.py`,
`pcfix.py`.
