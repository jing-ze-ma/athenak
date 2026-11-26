# Workflow of `<rad_vet>`: what runs, in what order, and where to read it

Open **`src/rad_vet/rad_vet_tasks.cpp`** first: the task lists and the implicit stage solve
(`StageSolve`) as a driver of named steps. The equations are in [SCHEME.md](SCHEME.md), the keys in
[KEYS.md](KEYS.md). Paths below are relative to `src/rad_vet/`.

## 1. One time step (`Driver::Execute`, `driver/driver.cpp`; `rad_vet_imex.cpp`)

```mermaid
flowchart TD
  A[dt_retry checkpoint] --> B[rotation step, auto-BE, slope correction]
  B --> C{"Time2Active? (hesdirk2, valid K1)"}
  C -- yes --> D[Time2Stash] --> E[hydro stage s = 1, 2]
  E --> F["Time2Stage: Time2GasTime, Time2FormStage, implicit stage (dt_sub = g dt)"]
  F --> G{failed?}
  G -- "no, s = 1" --> E
  G -- "no, s = 2" --> K[Time2Accept]
  G -- yes --> R[Time2Restore: back to U^n] --> H
  C -- no --> H[hydro stages 1, 2] --> I["BackwardEuler: implicit stage (dt_sub = dt), stores K1"]
  I --> L{dtr_fail?}
  K --> L
  L -- yes --> M[checkpoint back, smaller dt] --> B
  L -- no --> N[ledgers, outputs, restart]
```

## 2. One implicit stage (the `vet_*` task lists, `rad_vet_tasks.cpp`)

```mermaid
flowchart TD
  T0[InitRecv] --> T2["Opacity (rad_vet_opacity)"] --> S0
  subgraph S["StageSolve"]
    S0[SolveBegin] --> S1["StageClosure (rad_vet_closure)"] --> S2["OldVector, FaceCoefficients"]
    S2 --> S3["PassSetup, FirstPass (rad_vet_newton)"] --> P1
    P1["PassClosure: closure, opacity, donors of the iterate"] --> P2["PassRows: right-hand side, operator (rad_vet_rhs, _solver)"]
    P2 --> P3["PassSolve: E (rad_vet_solver)"] --> P4["PassUpdate: gas T, face F, momentum, drag (rad_vet_newton, _flux)"]
    P4 --> P5{PassSettled?}
    P5 -- no --> P1
    P5 -- yes --> V1["SolveEnd, VerdictSolve"] --> W["UpdateState = write-back (rad_vet_update)"]
    W --> V2[VerdictState: floors, gas eint] --> C2["ClosureAfterSolve: ray moments of the new state"]
  end
  C2 --> T3["GasUpdate: floor caps, gas ghosts, primitives (rad_vet_gas)"] --> T4[ghosts of E, F] --> T5[ClearSend, ClearRecv]
```

The face flux solve (`FluxSolve`, `rad_vet_flux.cpp`): F^H and the dissipation -> line ends across
MeshBlocks (`rad_vet_flux_schur.cpp`) -> x1-line solve -> shared faces and panel seams (`rad_vet_flux_seam.cpp`)
-> bound |F| <= c E and the steady balance (`rad_vet_flux_bounds.cpp`).

## 3. The closure (`rad_vet_closure*.cpp`, `rays/`)

```mermaid
flowchart TD
  A["BuildClosure | ClosureExtrap | ClosureIterate"] --> B["PredictState: E, T, opacity at t + dt"]
  B --> C[TraceRays]
  C -- Cartesian --> D["ScCartBuild (rays/sc_cart_sweep)"] --> E["ClosureFromRays"]
  C -- "spherical polar, cubed sphere" --> F["ScSource (rays/sc_build, sc_source)"] --> G{rotating ring?}
  G -- no --> H["ScSweep (rays/sc_sph_sweep) | ScCsSweep (rays/sc_cs_sweep)"] --> I["ClosureFromSc (closure_rays), face sums (closure_faces)"]
  G -- yes --> J["RingMoments (closure_ring): sweep two sets, shift, mix"]
  E --> K[ClosureFill] --> L[ghosts]
  I --> K
  J --> K
  A -. "extrapolated: no sweep" .-> K
```

## 4. The linear solve for E (`rad_vet_solver.cpp`)

```mermaid
flowchart TD
  S["PreconditionerSetup: face weights, rows, lines, coarse levels"] --> R[r = B - A X]
  R --> P[p] --> M1["y = M^-1 p: V-cycle of x1-line sweeps"] --> A1["v = A y: ghosts, rows, cross terms"] --> B1[alpha, s]
  B1 --> M2["z = M^-1 s"] --> A2["q = A z"] --> C[omega, X, r] --> T{"converged (true residual)?"}
  T -- no --> P
  T -- breakdown --> R
  T -- yes --> E[iteration count]
```

## 5. Files

| file | what |
|---|---|
| `rad_vet_tasks.cpp` | task lists, `StageSolve` |
| `rad_vet_imex.cpp`, `rad_vet_fct.cpp` | IMEX stages; limited face energy fluxes |
| `rad_vet_opacity.cpp`, `rad_vet_gas.cpp`, `rad_vet_newdt.cpp` | opacity; gas after the solve; time step |
| `rad_vet_closure.cpp`, `_closure_rays`, `_closure_faces`, `_closure_ring` | closure: state and choice; cells; faces; rotating sets |
| `rad_vet_oldvector.cpp`, `rad_vet_rhs.cpp`, `rad_vet_newton.cpp` | stage start; right-hand side; passes |
| `rad_vet_solver.cpp` | operator, preconditioner, BiCGStab |
| `rad_vet_flux.cpp`, `_flux_schur`, `_flux_bounds`, `_flux_seam` | face fluxes |
| `rad_vet_update.cpp`, `rad_vet_ledger.cpp` | write-back; L_z and energy ledgers, reports |
| `rad_vet_bvals.cpp`, `rad_vet_reduce.hpp` | ghost exchange, global sums |
| `rad_vet.cpp`, `rad_vet_params*.cpp`, `rad_vet_init.cpp`, `rad_vet_geometry.cpp` | constructor; keys; arrays; geometry |
| `rad_vet_steady.cpp`, `rad_vet_restart.cpp` | steady start; restart and dt_retry |
| `rad_vet.hpp`, `rad_vet_defs.hpp` | the class; slots and small structs |
| `rad_vet_args.hpp`, `_faces`, `_rows`, `_momentum`, `_flux`, `_moments` `.hpp` | device functions |
| `rays/sc_cart_*` | short characteristics, Cartesian |
| `rays/sc_sph_*` | short characteristics, spherical polar |
| `rays/sc_cs_*` | short characteristics, cubed sphere |
| `rays/sc_*` (others) | shared by both curved meshes: directions, tables, source, build, set-up, band ring |

Every file starts with what it does, from which files it is called and into which it calls.
