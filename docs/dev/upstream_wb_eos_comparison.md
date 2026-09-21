# Upstream vs. ours: well-balanced schemes and general EOS

Read-only comparison, 2026-09-21. Our branch: `rt-integration` @ `8d09252b`.
Upstream: `origin` = IAS-Astrophysics/athenak; `origin/main` @ `c5a0d7f9`;
`origin/wb-scheme` @ `ee085659` (2026-08-26, Jacob Fields, 44 commits ahead of
`origin/main`, merge-base `432e06b8`). Every claim below is quoted with a
file:line on the named ref; anything not read directly is marked *(inferred)*.

## Part 1 — upstream well-balanced scheme (`origin/wb-scheme`)

### Scope

Branch diff (`git diff --stat 432e06b8 origin/wb-scheme`, 22 files, +2574/-43) is
mostly *not* well-balancing: it also carries an orthonormal-tetrad unit test, a
transforming HLLE and a new HLLD solver for DynGRMHD. The well-balanced part is
confined to:

- `src/dyn_grmhd/rsolvers/wb_equilibrium.hpp` (new, 530 lines) — the scheme;
- `src/dyn_grmhd/dyn_grmhd_fluxes.cpp` — call sites;
- `src/dyn_grmhd/dyn_grmhd.cpp` / `.hpp` — parameters, source-term order;
- `src/dyn_grmhd/dyn_grmhd_util.hpp` — the interface FD/interpolation operators.

It applies **only to `src/dyn_grmhd` — GRMHD on a dynamical (Z4c/ADM) spacetime**.
The lapse `adm.alpha` and the spatial metric `adm.g_dd` are the only "gravity" it
knows; there is no Newtonian-gravity or fixed-metric (`cartesian_ks`) path, and
nothing in `src/hydro`, `src/mhd` or `src/srcterms` is touched. It is equally
valid for a static or an evolving metric, because the equilibrium is rebuilt from
the current ADM variables every stage *(inferred from the call site being inside
`CalcFluxes`)*.

### Equilibrium model

The assumed equilibrium is **static GR hydrostatic balance**, i.e. the
momentum-constraint balance

  ∂_i (α √γ P) = − √γ (e) ∂_i α + (metric-derivative term),

discretised cell-by-cell. There is **no column integration, no reference state,
no cache**: for each cell *i* the equilibrium pressure at its two interfaces is
obtained *algebraically* from the cell's own `(n, P, T, Y)` and the local metric,
`wb_equilibrium.hpp:42-86` (`PressureEquilibrium`):

```
e0 = eos.GetEnergy(n0, T0, Y);                              // line 47
q0 = g^{ab}_0 (g_ab,±1 − g_ab,0);                           // lines 59-61
pe = (α0 √γ0 P0 (1 + q0/4) − 0.5 √γ0 e0 (α_{±1} − α0)) / (α_{1/2} √γ_{1/2});
```

so the balanced interface pressures are `pr_i = pe(−)`, `pl_{i+1} = pe(+)`. This
is a *local, sub-cell* hydrostatic reconstruction in the Käppeli–Mishra spirit,
but the sub-cell equilibrium is neither isentropic nor isothermal: it is
"constant `(n, T, Y)` of the cell centre", with the total energy density `e`
supplying the source weight. If `pe ≤ 0` the scheme **falls back to donor cell**
(`pr_i = pl_{i+1} = P0`, return `false`, lines 63-66, 79-82) and the caller leaves
the ordinary reconstruction untouched (`if (balanced)` guards, lines 215-218 etc.).

Only **pressure** is well-balanced. Density, velocity, B and the scalars go
through the ordinary reconstruction.

### Role of the EOS / "temp saved during C2P"

`e0 = eos.ps.GetEOS().GetEnergy(n0, T0, Y)` needs a temperature. Commit
`52bc5fe7` ("Use temp saved during C2P inside equilibrium") replaced an in-kernel
`GetTemperatureFromP` call with a read of the `DvceArray5D temperature` array
that `DynGRMHD` already fills during con2prim (`dyn_grmhd.cpp:205`,
`Kokkos::realloc(temperature, nmb,1,ncells3,ncells2,ncells1)`); the WB kernels
take it as `const DvceArray5D<Real> &temp` and read `temp(m,0,k,j,i)`
(`wb_equilibrium.hpp:188`, `197`, `335-337`). This is a pure cost/robustness
optimisation — one fewer EOS inversion per cell per direction. It also means the
WB scheme inherits whatever EOS policy is compiled in (ideal, piecewise
polytrope, CompOSE, hybrid — the instantiation list at
`dyn_grmhd_fluxes.cpp` tail).

### Reconstruction coupling

- **DC**: `BalancePressureDCX1/2/3` (`wb_equilibrium.hpp:184`, `228`, `272`)
  simply *overwrite* `ql(IPR)`, `qr(IPR)` with the equilibrium values.
- **PLM and everything above DC** (`plm`, `ppm4/ppmx`, `wenoz`): after the normal
  reconstruction has run, `BalancePressureX1/2/3` (`wb_equilibrium.hpp:316`,
  `388`, `461`) *replace* the pressure component by an **equilibrium-preserving
  reconstruction of the deviation**:

```
PLM(q(i-1) - p_im1, 0.0, q(i+1) - p_ip1, ql(IPR,i+1), qr(IPR,i));  // line 375
ql(IPR,i+1) += pl_ip1;  qr(IPR,i) += pr_i;
```

  i.e. the equilibrium is extrapolated from cell *i* to its neighbours
  (`PressureEquilibriumExtrap`, lines 94-176, one fixed-point step per neighbour,
  lines 155 and 168), the deviation δP = P − P_eq is limited with the standard
  PLM limiter *with δP_i ≡ 0 by construction*, and P_eq is added back at the
  interfaces. Note the consequence: **whatever high-order method was selected,
  pressure ends up piecewise-linear** in the WB path. The `if(balanced)` guard
  means a cell that fails the positivity test keeps its original high-order
  pressure states.

### Source-term discretisation

The GR source term is unchanged in form — a cell-centred, non-conservative
`T^{μν} ∂_a (α, β^i, γ_ij)` contraction in `AddCoordTermsEOS`
(`dyn_grmhd.cpp:568-...`), using centred FD `Dx<NGHOST>` of order 2/4/6. What the
branch adds is that the order is now an *input parameter* rather than tied to
`nghost`: `dyn_grmhd.cpp:405-415` maps `source_order = 2/4/6 → Dx<2>/Dx<3>/Dx<4>`
(commits `1730918f`, `f2a73eba`). The matching *interface* half-difference
operators `HalfDifferenceInterface<width,dir,sign>` live in
`dyn_grmhd_util.hpp:66-87` (width 1/2/3 ⇒ 2nd/4th/6th order), and the WB call
sites pick `width = source_order/2` (`dyn_grmhd_fluxes.cpp`, the `switch
(source_order_)` blocks). So the equilibrium condition is discretised to the
*same order as the source term* — that is the actual well-balancing requirement.

### The "hybrid method" and the "new formulation"

- **Hybrid** (`c2e9bfd2`, "Updated extrapolation, enabled hybrid method"): the
  earlier attempt had a fully bespoke `PiecewiseLinearWBX1/2/3` living in
  `src/reconstruct/plm.hpp` and `wenoz.hpp`. `c2e9bfd2` commented those out and
  switched to the *hybrid* arrangement described above — run the standard
  reconstruction first, then correct only `IPR` by the equilibrium + deviation
  PLM. `57a59951` then deleted the dead `wenoz.hpp` WB code (−58 lines) and the
  dead branches in `dyn_grmhd_fluxes.cpp` (−107 lines), and also changed the
  WB-DC path (commit message: *"WB-DC changed; needs testing"*).
- **New formulation** (`ee085659`, branch tip, 2026-08-26, *"Started working on
  new formulation"*): adds, but does not yet use, **interface interpolation**
  operators `InterpToInterface<width,dir,sign>` (`dyn_grmhd_util.hpp:89-113` and
  `148-...`, 2nd/4th/6th order) and volume-averaged derivative operators
  `Dxvol<NGHOST>` in `src/utils/finite_diff.hpp:67-...`. Nothing calls them yet —
  the work is unfinished. *(inferred from grep: no call sites on the branch.)*

### Input parameters

All in `<mhd>` (`dyn_grmhd.cpp:155-192`):

| parameter | default | meaning |
| --- | --- | --- |
| `well_balanced` | `false` | enable the scheme |
| `source_order` | `2*(nghost-1)`, minus 2 if `fofc=true` | FD order of source + equilibrium (must be even, ≥2, ≤ max) |
| `rsolver = hlld` | — | new HLLD; WB warns if not used |
| `monitor_failures` | `false` | count HLLD failures (HLLD only) |

### Stated limitations

- `dyn_grmhd.cpp:156-160`: *"Well balancing may not work without HLLD"* — a
  warning, not an error, if `rsolver != hlld`.
- `dyn_grmhd.cpp:163-175`: `well_balanced` + `fofc` together cost 2 ghost zones,
  so `source_order` is capped; too few ghosts is a fatal error.
- Commit `aff7ba96`: *"Works in 1D for DC/PLM, 2D/3D and WENOZ have problems"*;
  `3b79952f`: *"Bug fixes. Scheme works in 3D now."*
- `57a59951`: WB-DC *"needs testing"*.
- FOFC cannot use the new solvers: it falls back to `hlle_dyngr`
  (`dyn_grmhd_fluxes.cpp`, the `FOFC<...>` dispatch).
- Metric index orders are permuted in the x2/x3 kernels "to reduce floating-point
  symmetry errors" (`wb_equilibrium.hpp:417`, `490`).

### Tests

**None.** `tst/` is untouched by the branch; no `<mhd>/well_balanced` appears in
any `inputs/` file (`inputs/dyngr/whisky_tov.athinput`, `mag_tov.athinput` are the
intended targets but leave it at the default `false`). The only new input file is
`inputs/unit_tests/tetrad_test.athinput`, for the unrelated tetrad test. The
`dyngr_tov` pgen's ghost-zone potential loop was commented out in `3b79952f`
(`src/pgen/dyngr_tov.cpp:215`, `:357`) — i.e. the branch is a work in progress.

## Part 2 — upstream general EOS (`origin/main` @ `c5a0d7f9`)

### `src/eos/primitive-solver/` — an embedded copy of Fields' PrimitiveSolver library

19 files, `namespace Primitive`. `EOS<EOSPolicy, ErrorPolicy>` (`eos.hpp:66`)
multiply inherits an EOS policy and an error policy; the wrappers do
code-units<->eos-units conversion around each policy call (`eos.hpp:157-400`).
The required policy API is listed at `eos.hpp:13-26`
(`TemperatureFromE/FromP`, `Energy`, `Pressure`, `Entropy`, `Enthalpy`,
`MinimumEnthalpy`, `SoundSpeed`, `SpecificInternalEnergy`,
`Minimum/MaximumPressure`, `Minimum/MaximumEnergy`). Shared state
(`eos_policy_interface.hpp:19-44`): `n_species, mb, min_n/max_n, min_T/max_T,
min_Y[]/max_Y[]`.

**EOS policies (4).**
- `IdealGas` (`idealgas.hpp:21`) — Γ-law; `SetGamma` clamps to (1.00001, 2.00001)
  (`:118-121`), which is itself buggy at Γ=2 (see `origin/hotfix/dyngr-gamma`).
- `PiecewisePolytrope` (`piecewise_polytrope.hpp:34`) — up to `MAX_PIECES` cold
  pieces + an ideal thermal part `gamma_thermal`; init from the input block or a
  file (`:154`, `:170`).
- `EOSCompOSE<LogPolicy>` (`eos_compose.hpp:32`) — 3-D tabulated; `LogPolicy` is
  `NormalLogs` or `NQTLogs` ("not-quite-transcendental" bit-twiddling log2/exp2,
  `logs.hpp:54,68,198-211`), selected by `<mhd>/use_NQT`.
- `EOSHybrid<LogPolicy>` (`eos_hybrid.hpp:32`) — 1-D cold table (uniform in
  log n_b) + Γ_th ideal thermal piece.

**Error policies: only `ResetFloor`** (`reset_floor.hpp:25`). There is no
`DoNothing` policy on main; `<mhd>/dyn_error` must be `reset_floor`
(`dyn_grmhd.cpp:79-105` fatals otherwise).

**Units**: `unit_system.hpp:19` (`c,G,kb,Msun,MeV` + length/time/density/mass/
energy/pressure/temperature/chemical-potential factors), factories `MakeCGS,
MakeGeometricKilometer, MakeGeometricSolar, MakeNuclear, MakeMKS` (`:101-174`).

### Where it is usable — dynamical GR only

`git grep -l primitive_solver_hyd.hpp origin/main -- src` returns only
`src/dyn_grmhd/*`. Despite its name, `PrimitiveSolverHydro`
(`src/eos/primitive_solver_hyd.hpp:44`) reaches straight into
`pmy_pack->pmhd->{nmhd,nscalars,bcc0,bfc}` (`:228`, `:313`, `:404`) — it is
DynGRMHD-only. Newtonian and SR/static-GR hydro/MHD keep the classic
`EquationOfState` hierarchy: `<hydro>/eos` ∈ {`ideal`, `isothermal`}
(`hydro.cpp:49-66`) and `<mhd>/eos` ∈ {`ideal`, `isothermal`}
(`mhd.cpp:74-106`); for dynamical GR the `<mhd>/eos=ideal` path installs the
placeholder `NoOpDynGRMHD` (`mhd.cpp:79-81`). The real selector is
`<mhd>/dyn_eos` ∈ {`ideal`, `piecewise_poly`, `compose`, `hybrid`}
(`dyn_grmhd.cpp:79-105`); 6 EOS × `ResetFloor` are explicitly instantiated
(`dyn_grmhd.cpp:689-722`).

**So: upstream has no general/tabulated EOS for Newtonian hydro or MHD at all.**

### Table format and interpolation

Not HDF5. Tables go through `TableReader::Table::ReadTable`
(`src/utils/tr_table.cpp:55-145`): an ASCII header (`<metadata>`, `<scalars>`,
`<points>`, `<fields>`) followed by raw binary `double` arrays. A CompOSE table
must be converted offline *(inferred)*.

Axes (`eos_compose.cpp:46-99`): 3-D, `nb` × `yq` × `t`, with `nb` and `t` stored
as log2 and `yq` linear; **the table must be uniformly spaced** in log n_b, log T
and y_q (`eos_compose.hpp:15-16`) — the inverse spacings are taken from the first
two points (`:74,84,96`). Tabulated fields (`eos_compose.hpp:39-47`): log2 P,
entropy/baryon, μ_b, μ_q, μ_l, log2 e, c_s, built from CompOSE Q1..Q7 and cs2
(`eos_compose.cpp:101-185`). Storage is `m_table(iv, in, iy, it)`.

Interpolation is **trilinear** (`eos_compose.hpp:337-400`). Inverse lookup
`temperature_from_var` (`:403-487`) is bilinear in (n, y_q) at fixed T index, then
a **bisection over the T index** plus a linear interpolation in log T (`:485`);
out of range it clamps to `min_T`/`max_T` (`:445-449`) and otherwise asserts
(`:455-467`). `m_min_h` is found by a full table scan at T_min with an NQT-aware
safety factor (`eos_compose.cpp:195-251`, PR #736).

### c2p algorithm

`primitive_solver.hpp:13-14` states it follows **Kastaun, Kalinani & Ciolfi, PRD
103, 023018 (2021)**: a 1-D root solve in **μ = 1/(hW)**.
- Bracket: `mul = 0`, `muh = 1/h_min`; if `r² > h_min²` the upper bound is
  tightened by a safe Newton solve on `UpperRoot` (tol 1e-10) and nudged by
  `(1+1e-10)` (`:456-480`). Failure → `Error::BRACKETING_FAILED`.
  `CheckDensityValid` (`:289-328`) uses `W_max = sqrt(1 + r²/h_min²)` and the EOS
  density bounds.
- Solve: **Anderson–Björck modified false position**,
  `root.FalsePosition(RootFunction, mul, muh, mu, tol, ...)` (`:495-496`) — not
  Brent/TOMS748 (`numtools_root.hpp` also has `Chandrupatla` and `NewtonSafe`,
  unused for μ).
- Tolerance/iterations: `<mhd>/c2p_tol` (default 1e-15), `<mhd>/c2p_iter`
  (default 50) (`primitive_solver_hyd.hpp:139-140`); relative convergence
  (`numtools_root.hpp:72,103`). `SolverResult::iterations` reports the *cap*, not
  the count (`primitive_solver.hpp:497-499`) — a live wart.
- Failure: `HandleFailure` (`:279-285`) → `ResetFloor::FailureResponse` sets
  atmosphere and re-densitizes via `PrimToCon`. Errors enumerated in
  `ps_error.hpp:12-23`; printed up to `<mhd>/c2perrs` (default 1000) times; in the
  floors-only (FOFC) pass a failure instead raises the FOFC flag
  (`primitive_solver_hyd.hpp:453-454`).

### Temperature and composition

`MAX_SPECIES = 3` (`ps_types.hpp:17`). **Y_e is an ordinary passive scalar**:
`Y[s] = cons[CYD+s]/cons[CDN]` (`primitive_solver.hpp:347-351`), stored at index
`nhyd + s` in the generic MHD scalar slots — there is no dedicated `IYE`.
`EOSCompOSE::SetNSpecies` accepts only 1 (`eos_compose.hpp:313-322`).

Temperature is **not** an evolved primitive. It is cached per cell each ConsToPrim
in `DynGRMHD::temperature` (`dyn_grmhd.hpp:104`, filled at
`primitive_solver_hyd.hpp:522`) and used by the coord-source and FOFC kernels —
and, on `origin/wb-scheme`, by the WB equilibrium. It is *not* used to warm-start
the next c2p; PrimToCon recomputes T from P (`:190`, `:270`) *(inferred)*.

### Floors

From `<mhd>` in `primitive_solver_hyd.hpp:131-157`: `dfloor` (→ `n_atm`),
`tfloor` (→ `T_atm`), `dthreshold` (→ `n_threshold`, default 1.0), `gamma_max`
(default 50 → `v_max`), `max_bsq` (default 1e6), `s<N>_atmosphere`.
`ResetFloor::PrimitiveFloor` (`reset_floor.hpp:35-52`) resets to atmosphere **and
zeroes the velocity** below `n_atm*n_threshold`; `ConservedFloor` (`:55-73`) uses
`tau_floor = e(D/mb, T_atm, Y) - D + 0.5 B²` (`eos.hpp:472-473`);
`MagnetizationResponse` (`:76-88`) rescales b if `b² > max_bsq` and returns
`CONS_ADJUSTED` (sign fixed in `df04c52a`, PR #780). Defaults
`fail_conserved_floor = fail_primitive_floor = false`, `adjust_conserved = true`.

### Recent EOS activity upstream

`df04c52a` (#780) MagnetizationResponse sign fix; `06672426` (#774) PWP zero
ε-offset; `55e7ec5e` (#749) modified energy limiting; `1b1e4e39` (#736)
conservative `h_min` for tabulated EOS. Branches:
`origin/projects/eosevolve` (1 commit `94e3001a`, +22 lines in
`src/pgen/dyngr_tov.cpp`: a user source term that relaxes the ideal-gas Γ from
`eos_soft_gamma_0` to `eos_soft_gamma_1` over `problem/eos_soft_tau` — a
time-dependent EOS demo, no framework change);
`origin/tov-refactor-tabeos` (7 commits, `+865/-76`, the `src/utils/tov/*` TOV
solvers and the 4-D table layout — already on main via a squashed PR, so the
branch is stale, *inferred*); `origin/hotfix/dyngr-gamma` (2 commits, one-line
`idealgas.hpp:119` clamp fix for Γ=2, not merged).

## Part 3 — ours (`rt-integration` @ `8d09252b`)

### Well-balanced scheme

**Scope**: Newtonian **hydro and MHD**, mirrored implementations
`src/hydro/hydro_wellbalance.cpp` and `src/mhd/mhd_wellbalance.cpp`
(`SetWbBackgroundPressure`, `RemoveWbFlux`, `AddWbVar`/`RemoveWbVar`,
`BuildWBCache`, plus `EnableWBEffectivePotential`/`SetWBEffectivePotential` at
`hydro_wellbalance.cpp:356`, `:384`). No relativistic path. Two schemes share the
`wb_` prefix: **static WB** (`wellbalance_static`, a frozen pgen-supplied
background `u0wb/w0wb/w0facewb`, `hydro.hpp:233-235`) and **dynamic WB**
(`wellbalance_dynamic`), the Käppeli–Mishra deviation scheme, which is the one to
compare with upstream.

**Equilibrium**: no stored column either — the background is rebuilt *locally per
cell* by walking the potential, `src/utils/wb_background.hpp` (541 lines, design
note at `:8-61`). The closure is selectable, `WBOption` (`src/athena.hpp:84`):
`isodensity | isothermal | isentropic | adaptive | polytropic | isentropic_dt |
adaptive_fast`; `adaptive*` compares `|dlnT|` with `ds/c_v` across the stencil
(`wb_background.hpp:172-182`). `WBAdvance()` (`:274-415`) integrates one segment in
Φ: isodensity closed form (`:352`), isothermal frozen-coefficient exponential
(`:371`), polytropic midpoint with T linear in Φ (`:312-351`), `isentropic_dt`
adiabat in (d,T) via Γ₁ (`:287-311`), `isentropic` exact enthalpy projection
h+Φ=const by Newton (`WBEnergyFromEnthalpy`, `:210-241`).
`WBBackgroundStencil()` (`:424-539`) returns `im1, imh, i, iph, ip1` by walking
**two segments per side** and does exactly **one temperature inversion per
stencil** (`:447`), seeded from the c2p cache. Ideal-gas closed-form fast path in
`Hydro::getWBerho` (`hydro.hpp:1147-1300`). `WBGuard()` (`:88-109`) flattens the
stencil to the anchor if any walked state is non-finite, non-positive, or outside
1e±6 of the anchor.

**Effective potential**: `phicc_wb` / `phi_wb_x1f` (`hydro.hpp:186-187`) are by
default *shallow copies* of the true `phicc0`/`phi0.x1f` (`hydro.cpp:214-217`), so
the default is bitwise the true potential; `EnableWBEffectivePotential()` allocates
separate Views so Φ_eff = Φ − ∫a_rad dz can be used
(`<problem>/wb_phi_eff`, `wb_arad_force`, `wb_arad_file`,
`src/pgen/box_convection.cpp:1767-1777`; filled at `:2942-2948`). Centrifugal
potential is not WB-specific: `deep_hot_jupiter_rt.cpp:233,722,2433` adds it to
`phicc0`, so WB sees it only via the default aliasing *(inferred)*.

**Cache**: `wbq0` is `(nmb, 15, k, j, i)` = 5 stencil points × (d, e, p)
(`hydro.cpp:209`, `mhd.cpp:240`), built by `BuildWBCache`
(`hydro_wellbalance.cpp:286-332`) over `i ∈ [is-1, ie+1]`. Refresh policy
(`hydro_fluxes.cpp:170-180`, `mhd_fluxes.cpp:176-179`): rebuild if
`wb_cache_every <= 0`, or `!wb_cache_built`, or `stage == 1 && ncycle %
wb_cache_every == 0`. So the reference state is **time-dependent and lagged by up
to `wb_cache_every` cycles** — the known f-mode-growth trap (memory:
`fmode-wb-cache-culprit`).

**Reconstruction**: deviations `q1 = q − q0` for `IEN` (always) and `IDN` (if
`wb_rho`), PLM-limited, with `q0_iph/q0_imh` added back
(`WbLocalPiecewiseLinearX1`, `hydro.hpp:390-412`); fallback to plain PLM of the
full state at a local extremum or a negative interface value (`:409-411`). Under a
general EOS the **pressure channel** is deviation-reconstructed too
(`WbPiecewiseLinearDerX1`, `hydro.hpp:658`) because the solvers consume
`wder(IDPR)` (`hydro_fluxes.cpp:239-243`). **Only PLM**: with WB on, the
`ReconstructionMethod` switch is bypassed (`hydro_fluxes.cpp:204-233`) — `dc`,
`ppm4`, `ppmx`, `wenoz` are silently ignored, with no warning. Non-uniform-grid
`PLM_nonuniform` (`hydro.hpp:788-806`) on spherical polar and cubed sphere.

**Source term**: the gravity source *is the background's own pressure difference
across the cell*, from the same cache (`wb_background.hpp:46-53`). Cartesian:
`bdt*(pr − pl)/dzc` (`box_convection.cpp:3312-3327`); curvilinear: area-weighted
face form `bdt*(area1(i+1)*(pr−p) + area1(i)*(p−pl))/volume`
(`red_giant.cpp:3405`). `<problem>/wb_grav_source = background | plain`
(`red_giant.cpp:1566-1576`). Static WB uses `−g*(ρ − ρ_bg)` (`:3378-3380`). The
curvilinear *geometric* source has its own WB form `SrcTermsCurvilinearWB`
(`coordinates.cpp:1350`, geometry cached by `BuildWBGeometry()`, `:1221`).

**Geometries**: Cartesian (`WbLocalPiecewiseLinear*`), spherical polar and cubed
sphere (position-aware `GridPiecewiseLinear*`, honouring `wb_rmin/wb_rmax`). On a
cubed sphere `RemoveWbFlux` multiplies the angular face background pressures by
that face's `sin_xi`/`sin_eta` (`hydro_wellbalance.cpp:108-119`). Guards: missing
or unknown `wb_option` aborts (`hydro.cpp:220,236-241`), no `wb_direction` aborts
(`:243-249`), static+dynamic together aborts (`red_giant.cpp:2807`). No guard for
WB + non-PLM reconstruction, WB + AMR, or WB + cubed sphere.

**Parameters** (`<hydro>`, same names in `<mhd>`, `hydro.cpp:175-252`):
`wellbalance_dynamic` (false), `wb_x1/x2/x3` (false), `wb_rho` (false),
`wb_cache_every` (int, 0 = every stage), `wb_rmin`/`wb_rmax` (0 = off),
`wb_option` (**required**, no default), `wellbalance_static`,
`wellbalance_static_reconst`, `etotgrav`.

**Restart**: the cache is *not* restart state; `ncycle` is, which caused the
silent-unbalanced-restart bug (`docs/handover/NOTE-2026-09-09-wb-restart-bug.md`,
KE_r ×614 on red_giant). Fixed as `5c0b98e4` with `wb_cache_built`
(`hydro.hpp:199`, `mhd.hpp:257`), plus invalidation when a restart file is written
(`src/outputs/restart.cpp:930-936`). The potential itself is not restart state —
the pgen must refill it.

**Tests**: `tst/test_suite/hydro/test_hydro_fofc_wb_cpu.py` +
`tst/inputs/fofc_wb.athinput` (pgen `wb_atm`), and
`inputs/tests/wb_column.athinput` / `wb_column_radcond.athinput` (pgen
`src/pgen/wb_column.cpp`). Production: `red_giant_*`, `he4_presn_*`,
`deep_hot_jupiter*`.

### General / tabulated EOS

`<hydro>/eos = general` (`hydro.cpp:66`), `<mhd>/eos = general` (`mhd.cpp:85`),
refused under SR/GR; classes `GeneralHydro`/`GeneralMHD`
(`src/eos/general_hyd.cpp`, `general_mhd.cpp`). A `<units>` block is mandatory
(`general_hyd.cpp:41-48`). Second switch `general_eos = gamma | table`
(default `gamma`, `src/eos/eos.cpp:200`) — `gamma` runs the whole general path
against an exactly known answer and is the bitwise regression basis.

**Table** (`src/eos/eos_table.{hpp,cpp}`): axes uniform in `x = log10 ρ[cgs]`,
`y = log10 T[K]` (`eos_table.hpp:9-11`); four surfaces log10(e_gas/ρ),
log10(p_gas/ρ), μ, log10(n_e/n_tot), each with value + ∂x + ∂y + ∂x∂y
(`ITE=0, ITP=4, ITMU=8, ITXE=12`, `:64`), stored `tbl(ny, nx, 16)` with the
variable index fastest (`:121-130`). Defaults `eos_logd_min=-14`,
`eos_logd_max=2`, `eos_logt_min=1.5`, `eos_logt_max=8`, `eos_dlog=0.05`,
`eos_dlogt = 0.2*eos_dlog` (`eos_table.cpp:175-188`).
Interpolation is **bicubic Hermite on tabulated analytic derivatives**
(`HermitePatch`, `:166`; `Interpolate`, `:224`), so `chi_rho`, `chi_T`, `c_v` and
`Γ₁ = chi_rho + p chi_T²/(ρ T c_v)` are exact *on the interpolant* (`:26-35`) —
exactly what the WB walk needs. Off-table: clamped patch coordinates plus an
explicit linear term (C1 extrapolation, `:37-44`).

**Composition**: analytic **Saha built at run time** on the host, once per node
(`src/eos/eos_composition.hpp`, 543 lines): H2/H/H+/He/He+/He++/e− plus a metal
donor, H2 dissociation equilibrium, closed by charge neutrality as a 1-D root find
on n_e (`:19-24`); energy = translational + H2 rot/vib + chemical (`:26-31`). No
Coulomb, degeneracy or pressure ionization (`:33-36`). Parameters `eos_xh`
(0.7381), `eos_yhe` (0.2485), `eos_a_metal` (16), `eos_h2`, `eos_ionization`,
`eos_metal_*` (`eos_table.cpp:85-103`).

**Radiation pressure**: `<block>/eos_radiation` (default false) adds `aT⁴` and
`aT⁴/3` analytically *outside* the table (`eos_table.hpp:291-298`), with an
optional density/temperature-gated taper `eos_rad_rho_hi/lo`, `eos_rad_t_hi/lo`
(`eos_table.cpp:108-109,138-139`, `rad_taper::WeightGated`, `eos_table.hpp:99-113`)
so the aT⁴ terms can hand over to a transport module; bitwise identical when off.

**T inversion and caches**: `SolveLog<MODE>` (`eos_table.hpp:372-404`) —
safeguarded Newton + bisection on the logs, 80 iterations, tol 1e-13 (double);
mode 0 = e at fixed ρ, 1 = p at fixed ρ, 2 = p at fixed T solving for log ρ.
`EvalEOnly` (`:422`) reads only the energy patch during the solve. `ClampLogT`
(`:592-604`) holds the answer inside the table and reports the clamp. Caches:
per-cell `wtemp` (`hydro.hpp:116`, `hydro.cpp:138`, `mhd.cpp:166`) written by
ConsToPrim and read back next stage as the **warm start**
(`general_hyd.cpp:123-124`, `:160-162`, 2-3 Newton steps), plus derived
`wder(IDPR) = p` and `wder(IDG1) = Γ₁` (`general_hyd.cpp:157-159`) so the Riemann
solvers never call the EOS. Fused accessor
`TemperaturePressureGamma1(d,e,tguess,t,p,g1,tclamped)` (`src/eos/eos.hpp:474`);
`ThermoAt`/`Enthalpy`/`ChiT`/`Cv`/`Gamma1` all take an optional known T so the WB
walk does not re-solve (`eos.hpp:315-520`).

**Floors** (`eos.cpp:25-140`): `dfloor`, `dfloor_keep_velocity`,
`dfloor_keep_temperature`, `vceil`, `vceil_thermalise`, `pfloor`, `tfloor`,
`eos_floor_consistent`, `efloor_from_ekin`, `efloor_as_tfloor`, `sfloor`
(refused with `general_eos = table`, `eos.cpp:206-215`), `defer_cons_floors`;
`tfloor_kelvin` (`:136-178`) with a `tfloor_set` sentinel so a restart's dumped
parameter list is not mistaken for a user value. `efloor_as_tfloor` raises
`tfloor` to `TableTempMin()` (`:220-241`).

**Restart**: the table is *not* in the file — rebuilt at setup from the input
(`eos_table.cpp:9-13`). The **caches are**: `wtemp`, `wder(IDPR)`, `wder(IDG1)`
(`src/outputs/restart.cpp:887-921`, sizes `:427-430`), and
`Hydro::c2p_freeze_derived` (`hydro.hpp:200-205`) makes the first post-restart
ConToPrim use `ConsToPrimFrozen` (`src/eos/general_hyd_frozen.cpp:52`) because the
table inversion is not idempotent (~5e-14 drift in 24% of cells, `:9-30`).

## Part 4 — comparison

| | upstream `wb-scheme` | ours (`rt-integration`) |
| --- | --- | --- |
| Target equations | DynGRMHD only (Z4c/ADM metric); no Newtonian path | Newtonian hydro **and** MHD; no relativistic path |
| Equilibrium model | static GR hydrostatic balance from the cell's own `(n,T,Y)`; no closure choice | choice of 7 closures (`isodensity … adaptive_fast`); plus a separate frozen-background "static WB" |
| How computed | algebraic, per cell, per stage, from `α`, `γ_ij` | walked in Φ, two segments/side, 5-point stencil, cached in `wbq0` |
| Reconstruction variables | pressure only (δP, with δP_i ≡ 0) | `IEN` always, `IDN` if `wb_rho`, `IDPR` under general EOS |
| Source discretisation | unchanged `T^{μν}∂_a g` FD, order 2/4/6 via `source_order`; WB stencil matched to it | source **is** the cached background Δp across the cell (Cartesian and area-weighted curvilinear forms) |
| Multi-D | x1/x2/x3 all supported; "works in 3D" since `3b79952f` | x1/x2/x3 selectable; only x1 gets cache/Φ_eff/rmin-rmax |
| EOS coupling | any PrimitiveSolver policy; uses the c2p `temperature` array | general tabulated EOS is a first-class partner: one inversion per stencil seeded by `wtemp`; ideal-gas closed forms as a fast path |
| Reference in time | recomputed every stage (no lag) | lagged by `wb_cache_every` cycles (0 = every stage) |
| Geometry | Cartesian only (the metric carries the curvature) | Cartesian, spherical polar, cubed sphere (incl. `sin_xi/sin_eta` face weights) |
| Reconstruction orders | DC path + a PLM correction applied on top of any method | PLM only; the method switch is bypassed |
| Tests | none | `test_hydro_fofc_wb_cpu.py`, `wb_column*.athinput` |
| Status | WIP; "needs testing"; new formulation started, unused | in production (red giant, He4 presupernova, hot Jupiter) |

| | upstream general EOS | ours |
| --- | --- | --- |
| Where usable | DynGRMHD only | Newtonian hydro/MHD only (refused under SR/GR) |
| Regime | nuclear (CompOSE, piecewise polytrope, hybrid), Y_e | stellar/planetary gas (H/He/metal Saha), no Y_e |
| Table axes | log n_b × y_q × log T, uniform, custom binary | log ρ × log T, uniform, built at run time |
| Interpolation | trilinear + T-index bisection | bicubic Hermite with tabulated analytic derivatives |
| Inversion | Kastaun μ = 1/(hW), Anderson–Björck false position, tol 1e-15, 50 it | safeguarded Newton+bisection on log T, tol 1e-13, 80 it, warm-started from `wtemp` |
| Radiation pressure | inside the table (CompOSE) | analytic `aT⁴` outside the table, optionally tapered |
| Caches | `temperature` (diagnostic; not a warm start) | `wtemp` + `wder(IDPR, IDG1)`, restart-persisted, frozen first post-restart c2p |

### What upstream does that we could adopt

1. **Matching the equilibrium stencil to the source-term FD order**
   (`source_order` ↔ `HalfDifferenceInterface<width>`): our WB is exact by
   construction for the cached background, but our reconstruction is hard-wired to
   2nd-order PLM. Upstream's `width = source_order/2` machinery is the right shape
   if we ever want 4th/6th-order WB.
2. **The "apply WB as a correction on top of whatever reconstruction ran"**
   pattern: it would let us keep `ppmx`/`wenoz` for the non-WB channels instead of
   silently disabling them (our `hydro_fluxes.cpp:204-233` bypass).
3. **Positivity fallback per cell returning a bool** (their `balanced` flag) is
   cleaner than our flatten-the-stencil `WBGuard`, and composes with a fallback to
   the *original* high-order states rather than to the anchor.
4. The PrimitiveSolver **unit-system class** and the `NQTLogs` fast log/exp are
   reusable independently of GR.

### What we do that upstream lacks

- WB for Newtonian gravity at all; WB on curvilinear and cubed-sphere grids; WB in
  MHD; a selectable equilibrium closure; an effective potential including radiative
  acceleration; a cached background with an explicit refresh cadence; a
  restart-correct cache; WB regression tests.
- A general EOS usable outside GR: Saha composition built at run time, bicubic
  Hermite with exact thermodynamic derivatives, analytic radiation pressure with a
  taper, warm-started log-T inversion with restart-persisted caches and the frozen
  first post-restart c2p.

### Incompatibility / merge surface

Files touched by **both** since `merge-base(HEAD, origin/wb-scheme) =
a1b2afab`, with (added/deleted) line counts ours | theirs:

| file | ours | theirs |
| --- | --- | --- |
| `src/driver/driver.cpp` | 349/39 | 10/0 |
| `src/dyn_grmhd/dyn_grmhd.cpp` | 40/40 | 83/11 |
| `src/mhd/mhd.cpp` | 286/6 | 7/0 |
| `src/mhd/mhd_tasks.cpp` | 228/20 | 1/1 |
| `src/parameter_input.cpp` | 51/9 | 1/1 |
| `src/pgen/tests/gr_bondi.cpp` | 0/6 | 14/6 |

The two well-balanced schemes do **not** collide in code: theirs lives entirely in
`src/dyn_grmhd/`, ours in `src/hydro`, `src/mhd`, `src/utils/wb_background.hpp`,
`src/coordinates`. The real conflicts are elsewhere (our `driver.cpp` and
`mhd*.cpp` changes for the RT/M1 work). The one genuine **namespace** collision is
the parameter name: upstream uses `<mhd>/well_balanced` (bool) for the GR scheme
while we use `<mhd>/wellbalance_dynamic` + `wb_option` for the Newtonian one — they
would coexist but confusingly. `<mhd>/source_order` is free on our side.
