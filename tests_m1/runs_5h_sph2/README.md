# runs_5h_sph2: hesdirk2 and implicit_vimp on the spherical-polar wedge and with vet_col; new defaults

- **Branch.** `m1-sph2` from rt-integration d9583121; worktree `/viper/ptmp2/jinma/wt_sph2`.
- **Work dir.** `/viper/ptmp2/jinma/sph2_0924`: `bin/`, `cpu/` (every run below), `gpu/`
  (job 11965525), `inp/`, `scripts/` (copied to `scripts/` here), the `*.log` result files.
- **Binaries (md5).**
  - ref = the m1def2 binaries of fe147c9f (d9583121 differs only in a wrapped comment),
    symlinked as `bin/athena_ref_{box_cpu,none_cpu,box_gpu}`.
  - new = snapshot of this branch: none cpu dc2c2feb, box cpu e9f58130, none gpu ff351811,
    box gpu e91b26d3. The committed source differs from the snapshot in comments only
    (one wrapped line, the Time2VetStart header).
  - `bin/athena_v2_none_cpu` (e6810580): the snapshot plus an experimental vet_col tensor
    extrapolation that is NOT committed (section 3). It was used, with extrapolation off,
    for `scripts/ss2.sh`, `tv3.sh`, `tvx.sh`.
  - gcc 14 + openmpi 5, Release, MPI; GPU rocm 6.3, gfx942, apudev, `HSA_XNACK=1`,
    `HSA_NO_SCRATCH_RECLAIM=1`.

## 1. Changes

| file | change |
| --- | --- |
| `rad_m1_sph.cpp` | `SphericalS1Check` no longer refuses `time_scheme != be` or `implicit_vimp` on the wedge. |
| `rad_m1_implicit.cpp`, `ImplicitVimpBuild` | sp overwrite blocks (`if (sph)`, Cartesian text untouched): the Jacobian rows of the face fluxes use `dxface` (x1, x2, x3) and, for a chi(f) closure, the S2 integrating factor `M1SphDrr` on x1 faces; the enthalpy-flux coefficients and their lagged right-hand side are recomputed with `dt A_f / V_i` per face (the Cartesian face loop runs, its `jr` part is put back to `jr0`, and an sp loop replaces it). |
| `rad_m1_implicit.cpp`, sp x1 row (`m1_impl_asm`) | The sp overwrite now adds the vimp terms (`J1M`, `JD`, `J1P`, `JRHS`), which it had dropped. It also uses `M1EnthCorrT` in a hesdirk2 stage solve under vimp, as the Cartesian row does. |
| `rad_m1_implicit.cpp`, `ImplicitSolve` | vet_col: the stage-1 solve builds the tensor and the surface q once per step at U^n (`Time2VetStart`); the stage-2 solve keeps them. A backward-Euler step builds them as before. |
| `rad_m1_time2.cpp` | `Time2Init` accepts vet_col (it sets `tau_closure`; the plain tau closure stays refused) and allocates `vet_now`/`vet_opac` for it. `Time2VetStart` calls `VetColBuild` for vet_col: E^n, T^n and the opacities of U^n, like vet_sc. No extrapolation: `Time2VetExtrapolate` runs for vet_sc only. |
| `rad_m1_implicit.cpp`, `ImplicitInit` | The default rule (below). |

Nothing in the hesdirk2 stage machinery itself carries geometry. The old vector
(`start + t2inc`) enters the face and row equations through `f0x*n` and `EN`, which the S1/S2
sp rows already use. The slopes are differences of physical-component states. The gas
momentum slope is the face-mean deposit, which is the same on sp.

**Defaults** (`ImplicitInit`, extending d9583121):
- `time_scheme` defaults to `hesdirk2` wherever:
  - transport = implicit;
  - the mesh is Cartesian **or the spherical-polar wedge** (no cubed sphere, no polar
    boundary);
  - the closure is not tau (**vet_col included**, on sp and on Cartesian);
  - `<time>/integrator = rk2`;
  - the run is not a restart.
- Elsewhere the default stays `be`. A restart whose file lacks the key keeps `be`.
  Explicit input overrides.
- `implicit_vimp` follows the unchanged d9583121 rule. It is on where the resolved scheme is
  hesdirk2 and its preconditions hold: multi-D, bicgstab, nghost >= 2, and hydro with
  coupling, gas_feedback and dbg_gas_force. That now includes the wedge.
- Nothing that ran before becomes fatal. The wedge requires rk2 and implicit transport
  already, and every closure accepted on sp is accepted by hesdirk2.

## 2. Gates

### Bitwise (ref vs new)

| gate | cases | result | log |
| --- | --- | --- | --- |
| A: Cartesian, `time_scheme = be` named, CPU | box 1 rank; box 2 ranks overlap+faces; slab 1; slab vet_sc full 2; slab2d_nd 2; box3d_nd 2; cart_sym, marshak_cart, rw_cart, pp_np (sp input) | 10/10 BITWISE | `gate_AB.log` |
| B: Cartesian, key absent (the current hesdirk2 + vimp defaults), CPU | slab2d_nd 2; box3d_nd 2; box3d_nd vet_sc 2; cart_sym; marshak_cart; rw_cart | 6/6 BITWISE | `gate_AB.log` |
| C: sp S1/S2/S5 + Cartesian vet_col (Milne) with `be` named | the 15 runs_5e/5g gate_sp cases | 15/15 BITWISE | `gate_CDE.log` |
| D: an old sp restart (ref, no key, i.e. be) continued by ref and by new | sph_sym_gas_rst, 4 ranks | BITWISE (26 files); new = ref straight run (rst) | `gate_CDE.log` |
| GPU, job 11965525, 2 GPUs | slab 1 rank be; box3d_nd 2 ranks be; + vet_sc; box3d_nd key absent (hesdirk2 + vimp) | 4/4 BITWISE (10 files each) | `gpu/log.out.11965525` |
| `tests_m1/gates/gates.py` (new box binary) | 14 pairs | GATES: PASS | `gates_eval.txt` |
| `tst/test_suite/rad_m1` (ATHENAK_M1_DATA = faces_0924/m1data) | 3 | 3 passed (1512 s) | `pytest.log` |

### hesdirk2 on the wedge (new binary)

**T-S1 transient** (`val_T1.log`). Setup:
- sph_sym: stretched r, tau 2000, Gaussian shell, 4 ranks, t = 100;
- dt from implicit_cfl 100/50/25/12.5; the reference is hesdirk2 at 3.125;
- tolerances 1e-12 / 1e-14.

Rel. L1 of E:

| closure | be | hesdirk2 | hesdirk2 orders |
| --- | --- | --- | --- |
| eddington | 3.0e-3 .. 3.8e-4 (1.00) | 8.3e-5 .. 1.3e-6 | 1.99 2.00 2.03 |
| m1 | 3.0e-3 .. 3.8e-4 (1.00) | 8.3e-5 .. 1.3e-6 | 1.99 2.00 2.03 |
| vet_col | 3.0e-3 .. 3.8e-4 (1.00) | 8.4e-5 .. 1.4e-6 | 1.98 1.98 1.99 |

**T-S5 radwave along theta at R = 100** (`val_T5.log`). Setup:
- 4x64x4, 2 periods, dt from hydro cfl 0.4/0.2/0.1/0.05;
- the reference is hesdirk2 at 0.0125;
- rel. L1 of rho - 1 over the amplitude.

| arm | errors | orders |
| --- | --- | --- |
| eddington, be | 5.6e-2 .. 1.1e-2 | 0.94 0.85 0.60 |
| eddington, hesdirk2 + vimp (the default) | 1.75e-3 .. 2.5e-5 | **2.03 2.03 2.08** |
| eddington, hesdirk2 without vimp | 2.8e-2 .. 8.7e-3 | 0.86 0.55 0.29 |
| m1 (offdiag lagged), be | 8.2e-2 .. 1.5e-2 | 0.8 |
| m1, hesdirk2 + vimp | 1.5e-2 .. 2.0e-3 | 1.06 0.55 1.24 |
| m1, the same wave on the CARTESIAN slab (`val_T5C.log`) | be 8.2e-2 .. 1.5e-2; hesdirk2 1.5e-2 .. 2.2e-3 | hesdirk2 1.04 0.58 1.15 |

- The M1-closure order is the same on the Cartesian slab to two digits, so it is a property
  of the chi(f) closure (lagged per pass), not of the sp geometry.
- Without vimp, hesdirk2 is not second order, as runs_3x found on Cartesian meshes. So this
  test also validates the sp vimp rows.

**Steady states unchanged** (`val_SS.log`, `val_SS2.log`, `val_SS3.log`). The metric is
max |h2 - be| / max |be| of E and F1.
- T-S1 steady diffusion: 4.6e-13 uniform, 8.3e-11 stretched.
- T-S2 free streaming at tol 1e-10: 1.4e-11. Both runs had 0 NC and 0 fallbacks.
  - At the input's tol 1e-12 every be solve is NON-CONVERGED: it stalls at 2e-11, as
    runs_5b noted.
  - Under hesdirk2 every stage is then refused and redone with backward Euler. The result
    is identical to be.
- T-S4 M1 (thermal, gas held): E 4.9e-7, F1 7.9e-7.
- T-S4 **vet_col pure scattering: 4.9e-10**.
- T-S4 vet_col thermal with the gas held does not agree at tolerance: 1.9e-4 at cfl 0.3 and
  5.3e-5 at cfl 0.075, in E, largest in the thin top.
  - Both states are steady to 1e-9.
  - The difference is O(dt), and it disappears without the thermal source. It is the timing
    of the source temperature: under hesdirk2 the tensor sees T^n; under be it sees the T
    after the hydro step.
  - The held gas is not a hydro steady state: the hold keeps e_int while the hydro step
    heats it. So T^n and the post-hydro T differ by one step of that heating.
- The same pattern holds for the stretched-r vet_col atmosphere (5.9e-5 at cfl 0.3) and the
  Cartesian Milne slab (3.5e-5).

**T-sym** (`val_SYM.log`, the default hesdirk2, 4 ranks, 500 steps, E changes by 0.80):
- spread(E)/E is 4.8e-13 Eddington, 4.6e-13 vet_col and 4.0e-14 m1 (tol 1e-12, lin 1e-14,
  line preconditioner);
- sym_gas at c = 100, vimp on by default: E 2.1e-14 / 2.4e-14, rho 4.3e-15 / 5.0e-15,
  v_r 3e-15.
- Each of the three cold-gas sym runs has one stage fallback at cycle 1:
  `min(E,T) = -1`, the gas energy check with T_gas = 1e-16. The step is redone with
  backward Euler.

**Restart with hesdirk2 on sp** (E, `gate_CDE.log`). sph_sym_gas_rst (Eddington and vet_col),
4 ranks, straight vs restart at 250: the rst payloads are BITWISE (2/2 each).

### hesdirk2 + vet_col

- **T-S4 steady state**: see above. It matches at tolerance without the thermal source.
  With the held-gas thermal source it differs by O(dt).
- **Second order in time**: T-S1 with vet_col, orders 1.98-1.99.
- **Tensor lag** (`val_TV2.log`, `val_TVX.log`, `val_TV3.log`). Setup: the T-S4 atmosphere,
  pure scattering, relaxing from the Eddington start to t = 6.4.
  - hesdirk2 errors are 2.4e-8 .. 1.8e-9, 20x below be (4.8e-7 .. 2.2e-8). Both orders are
    1.2-1.3.
  - An extrapolated chi and surface q (`D* = D^n + (dt/dt_prev)(D^n - D^{n-1})`,
    experimental binary v2) give 2.2e-8 .. 1.7e-9, the same orders. At the finest dt they
    differ from D^n by 3e-11.
  - The m1 closure gives the same orders.
  - So the tensor lag is not what sets the order here. The reduced order is common to the
    chi(f) and VET closures in this thin-top transient (probably order reduction in the
    stiff initial layer). The extrapolation was NOT committed; both stages use D^n, as
    vet_sc does by default.
- **Thin-cell seed test** (`val_THIN.log`; the runs_5b Finding-1 reproducer, kappa_s 25,
  seed 1e-6, 400 steps): the spread of E is 0 (below float32) at cycles 200 and 400.
  This holds on the wedge (radial and flux axis) and on the Cartesian twin. There were no
  fallbacks.

### vimp on the wedge

- **Risk 6** (`val_R6.log`). Setup: static pure-scattering M1 atmosphere,
  force_reference = wb_arad, gas free for 100 steps after 1000 held steps.

| arm | rms abs(v) | max abs(v) | max abs(T/T0 - 1) |
| --- | --- | --- | --- |
| be (reproduces runs_5b: 2.46e-5) | 2.46e-5 | 5.8e-5 | 1.7e-4 |
| be + vimp | 2.46e-5 | 5.8e-5 | 1.7e-4 |
| hesdirk2 | 7.5e-6 | 1.8e-5 | 4.9e-5 |
| hesdirk2 + vimp | 7.5e-6 | 1.8e-5 | 4.9e-5 |

  - No arm gets worse; hesdirk2 is 3.3x quieter.
  - There were 0 NC and 0 vimp positivity fallbacks.
  - vimp changes nothing visible here: the gas is static, and m1 uses the upwind enthalpy.
- **He wedge grid** (`val_HE.log`). Setup: hewedge at 96x32x32, stretched r, 4 ranks, gas
  free, wb_arad, 60 steps. Arms: be; hesdirk2 + vimp with vet_col; and the same with
  Eddington.
  - NON-CONVERGED = 0 in every arm.
  - 0 stage fallbacks, 0 vimp positivity fallbacks.
- T-S5 (above) is where vimp matters: it is needed for the second order.

## 3. GPU cost (job 11965525, 2 GPUs, apudev, same binary, interleaved, 2 repeats)

The He wedge grid is 96x128x128 with stretched r: 16 blocks of 96x32x32, gas held,
vet_col, 40 steps at the same dt.

| arm | time for 40 steps (s) | inner iterations per solve | solves |
| --- | --- | --- | --- |
| be | 2.014 / 1.904 | 27.2 | 40 (51 passes) |
| hesdirk2 | 1.921 / 1.909 | 10.7 | 79 (101 passes) |

- Per step hesdirk2 costs the same as be: **1.00x**. It takes two stage solves, but at
  g dt each with time2_lin_tol_fac 10, so the Krylov work is lower.
- The tensor is built once per step in both: 3.46 ms per build, 40 builds.

## 4. Still refused on sp

- the poles (polar boundary, a theta range reaching 0 or pi);
- the cubed sphere and the theta stretch;
- closure = vet_sc and tau;
- implicit_offdiag = operator;
- transport = explicit;
- more than one MeshBlock along x1;
- SMR/AMR;
- implicit_halo_mpi = true;
- implicit_flux != central;
- implicit_recon != dc;
- implicit_trans_limit != none;
- dbg_tensor.

hesdirk2 still refuses closure = tau (Cartesian) and integrators other than rk2.

## HANDOVER

State: m1-sph2 = rt-integration d9583121 + one commit (code, this README, `scripts/`), not
merged. Build trees are deleted; the binaries are kept in `/viper/ptmp2/jinma/sph2_0924/bin`.

Open points:
1. With a thermal source and the held gas, vet_col's source temperature differs between the
   schemes (T^n vs the post-hydro T), which moves the held-gas steady state by O(dt).
   - If this matters for production, the stage-1 build could take the stage-start T
     instead of T^n.
   - The runs of section 2 give the numbers to compare against.
2. The M1-closure radwave and the thin-top relaxation are about first order under
   hesdirk2, on Cartesian meshes as on the wedge. The errors are still 5-20x below be.
3. The one cold-gas stage fallback at cycle 1 (T_gas = 1e-16) is harmless and not
   sp-specific.
