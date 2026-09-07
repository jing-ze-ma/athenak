# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Session handover (read this first on a new machine)

Work moves between machines (viper, orion) that cannot see each other's run directories or
home directories. The state of the project at the last handover, what is running where, and what
to do first live in **`docs/handover/HANDOVER-<date>.md`** (newest date wins). The assistant's
memory notes as of that handover are copied under `docs/handover/claude-memory-<date>/`; on a
machine whose own memory directory is empty or older, read `MEMORY.md` there first and copy the
directory into the local memory directory so recall works. Analysis scripts that the notes refer
to are in `docs/handover/scripts/`. Never write into `run/` (165 GB of output).

## What this is

AthenaK: block-based AMR astrophysical fluid dynamics + numerical relativity code, written in C++17 on top of
[Kokkos](https://kokkos.org/) for performance portability (CPU / NVIDIA / AMD / Intel GPUs). It is a rewrite of
Athena++, not a port; only a subset of Athena++ features exist. Documentation lives on the
[wiki](https://github.com/IAS-Astrophysics/athenak/wiki).

`kokkos/` is a git submodule (`git submodule update --init --recursive` if missing). In-source builds are
rejected by CMake.

## Build

```bash
cmake -B build                      # configure (Release by default)
cd build && make -j                 # -> build/src/athena
```

CMake options (all off/default unless set):

| Option | Meaning |
| --- | --- |
| `-D Athena_ENABLE_MPI=ON` | MPI parallelism (also sets `Kokkos_ENABLE_MPI`) |
| `-D Athena_ENABLE_OPENMP=ON` | OpenMP |
| `-D Athena_SINGLE_PRECISION=ON` | `Real` = float instead of double |
| `-D PROBLEM=<name>` | compile `src/pgen/<name>.cpp` as the user problem generator (see *Problem generators*) |
| `-D Athena_ENABLE_RNS=ON` | link the RNS rotating-neutron-star library |
| `-D CMAKE_BUILD_TYPE=Debug` | also turns on Kokkos debug + View bounds checking |

GPU builds go through Kokkos flags, e.g.
`cmake -B build -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE80=On -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper`.
Some problem generators pull in external libraries (`z4c_two_puncture` → TwoPunctures+GSL, `z4c_spectre_bbh` →
SpECTRE `BundledExporter` via `-D SPECTRE_ROOT=`, `lorene_bns`, `elliptica_bns`, `sgrid_bns`); see the top-level
`CMakeLists.txt`.

**New source files must be added explicitly to `src/CMakeLists.txt`** — there is no globbing.

## Run

```bash
./build/src/athena -i inputs/tests/linear_wave_hydro.athinput      # run
./build/src/athena -i <file> mesh/nx1=128 time/tlim=0.5            # override any <block>/param
./build/src/athena -r <rstfile>            # restart (input embedded in restart file)
./build/src/athena -n -i <file>            # parse input and dump parameters, then quit
./build/src/athena -m -i <file>            # write mesh structure and quit
./build/src/athena -c                      # show build configuration
./build/src/athena -d <dir> -t hh:mm:ss     # run directory, wall-clock limit
```

Input files are the Athena++ `<block>`/`param = value` format; see `inputs/` (organized by physics). Python
readers/plotters for outputs are in `vis/python` (`athena_read.py`, `bin_convert.py`, `plot_slice.py`, …).

## Tests and lint

`tst/run_test_suite.py` is the current driver (used by `.github/workflows/main.yml`); it rebuilds into
`tst/build`, runs pytest over `tst/test_suite/`, then cleans up. Requires `pytest`, `numpy`, `flake8`.
Always invoke from `tst/`:

```bash
cd tst
python run_test_suite.py --style                     # cpplint + custom C++ checks + flake8
python run_test_suite.py --cpu                       # all *_cpu tests
python run_test_suite.py --mpicpu                    # all *_mpicpu tests (adds Athena_ENABLE_MPI=ON)
python run_test_suite.py --gpu "-DKokkos_ARCH_AMPERE80=On -DCMAKE_CXX_COMPILER=$PWD/../kokkos/bin/nvcc_wrapper"
python run_test_suite.py --test test_suite/nr/test_nr_lwave1d_cpu.py    # single test file
```

Device selection is **encoded in the filename**: a new test must contain `_cpu`, `_mpicpu`, or `_gpu`, or it
will never run. Convergence tests hard-code per-configuration L1 error and error-ratio thresholds (see e.g.
`tst/test_suite/nr/test_nr_lwave1d_cpu.py`); changing an algorithm usually means updating those tables.
`tst/run_tests.py` is the older harness still referenced by `.gitlab-ci.yml`.

Style is enforced, and CI fails on any violation (`tst/test_suite/style/check_athena_cpp_style.sh`, config in
`CPPLINT.cfg`): Google C++ style with **90-column** limit, no tab characters, no `}}` on a single line, no
trailing whitespace, `#pragma` left-justified, and `src/` files must be mode 644 (non-executable). Python is
flake8 with the same 90-column limit (`setup.cfg`).

## Architecture

### Data model: everything is a Kokkos View, batched over MeshBlocks

`src/athena.hpp` is the foundation and worth reading first. It defines `Real`, the `DvceArray1D..6D` /
`HostArray*` / `DualArray*` / `ScrArray*` View aliases (all `LayoutRight`), the face- and edge-centered field
structs (`DvceFaceFld4D/5D`, `DvceEdgeFld4D`), the `par_for` / `par_for_outer` / `par_for_inner` loop wrappers,
and the multi-value reducer used by history outputs.

`par_for` flattens N-D index ranges onto a 1D `RangePolicy` (MDRange is deliberately not used).
`par_for_outer` + `par_for_inner` gives thread teams with scratch memory for the inner vector loop.

Fluid state is stored as **5D arrays indexed `(m, n, k, j, i)`** where `m` is the MeshBlock index *within a
pack* and `n` the variable. This is the central design decision: one kernel launch covers all MeshBlocks on a
rank, which is what makes GPUs efficient.

### Mesh hierarchy

- `Mesh` (`src/mesh/mesh.cpp`) — global grid: `mesh_size`, `mesh_indcs`, `mb_indcs`, boundary flags, dt, and
  the flags every kernel branches on (`one_d/two_d/three_d`, `multi_d`, `multilevel`, `adaptive`).
- `MeshBlockTree` (`meshblock_tree.cpp`, `build_tree.cpp`) — octree of `LogicalLocation`s; built either from
  scratch or from a restart file, then blocks are distributed with `load_balance.cpp`.
- `MeshBlockPack` (`meshblock_pack.cpp`) — the per-rank container that owns `MeshBlock`s, `Coordinates`, every
  physics module pointer, and the map of task lists. Physics code almost always works through `pmy_pack`.
- `MeshBlock` — per-block metadata (sizes, `NeighborBlock` connectivity as `DualArray`s), not per-block data
  arrays.
- `MeshRefinement` (`mesh_refinement.cpp`, `refinement_criteria.cpp`) — SMR/AMR: refinement flags,
  prolongation/restriction interpolation weights, MeshBlock creation/deletion, and load rebalancing.

Startup order in `src/main.cpp` matters and is documented there in numbered steps: `Mesh` ctor →
`BuildTreeFromScratch/Restart` → `AddCoordinatesAndPhysics` → `ProblemGenerator` → `Driver` + `Outputs` →
`Driver::Initialize/Execute/Finalize`. Anything holding a Kokkos View must be destroyed before
`Kokkos::finalize()`, which is why main deletes objects explicitly at the end.

### Physics modules are selected by input blocks

`MeshBlockPack::AddPhysics()` constructs modules purely from *which blocks exist in the input file*
(`pin->DoesBlockExist("hydro")`, `"mhd"`, `"radiation"`, `"z4c"`, `"adm"`, `"ion-neutral"`, `"particles"`,
`"turb_driving"`, `"units"`, …), and enforces the legal combinations (e.g. `<hydro>` + `<mhd>` together
requires `<ion-neutral>`). Modules that can be coupled (ion-neutral, radiation+fluid, Z4c+DynGRMHD) do *not*
assemble their own task lists; the coupling module does, via
`src/tasklist/numerical_relativity.cpp` for the NR side.

Each physics directory follows the same layout: `<mod>.cpp` (ctor, parameter reading, allocation),
`<mod>_tasks.cpp` (task assembly + thin task wrappers), `<mod>_fluxes.cpp`, `<mod>_update.cpp`,
`<mod>_newdt.cpp`, plus `rsolvers/` for Riemann solvers. Reconstruction is header-only in `src/reconstruct/`
(`dc/plm/ppm/wenoz`, selected by the `ReconstructionMethod` enum).

### Execution: Driver + TaskList

`Driver` (`src/driver/driver.cpp`) sets up the SSP-RK / ImEx integrator coefficients (`gam0/gam1/beta/delta`,
`a_twid`) and drives the main loop. Per cycle it calls `ExecuteTaskList` on task lists keyed by name in
`MeshBlockPack::tl_map`:

```
"before_timeintegrator"
  for stage in 1..nexp_stages:  "before_stagen" -> "stagen" -> "after_stagen"
  (optional RKG super-time-stepping sub-loop for resistivity: "*_rkg_*")
"after_timeintegrator"
```

The semantics of the names (documented in `src/hydro/hydro_tasks.cpp`): `before_stagen` = work that must
complete over *all* MeshBlocks before the stage (post MPI receives, reset comm flags); `stagen` = the stage
itself; `after_stagen` = work only valid once every block finished (clear MPI requests);
`before/after_timeintegrator` = operator-split physics.

`src/tasklist/task_list.hpp` implements the whole task machinery in one header. A `Task` is a
`std::function<TaskStatus(Driver*, int stage)>` plus a `TaskID` (a `std::bitset<64>`) and a dependency mask;
`DoAvailable()` runs whatever is unblocked, and tasks may return `TaskStatus::incomplete` to be retried (this
is how nonblocking MPI is overlapped with compute). Max 64 tasks per list.

A typical `stagen` chain (hydro): `CopyCons → Fluxes → SendFlux → RecvFlux → RKUpdate → SrcTerms →
[orbital advection] → RestrictU → SendU → RecvU → [shearing box] → ApplyPhysicalBCs → Prolongate →
ConToPrim → NewTimeStep`.

### Boundary values and communication

`src/bvals/` is the most intricate part of the code. `MeshBoundaryValues` is an abstract base with
cell-centered (`bvals_cc.cpp`) and face-centered (`bvals_fc.cpp`) derived classes, plus a separate path for
particles (`bvals_part.cpp`). Buffer index ranges are precomputed per neighbor for each of the same-level,
coarser, finer, prolongation, and flux-correction cases (`MeshBoundaryBuffer` in `bvals.hpp`;
`buffs_cc.cpp` / `buffs_fc.cpp` compute them). Packing/unpacking happens in device kernels into 2D Views
`(nmb, ndata)`; MPI sends the same device buffers. MPI tags are built by `CreateBvals_MPI_Tag(lid, bufid)`
using **receiver-side** lid/bufid, with `NUM_BITS_LID = 14` bits for the local ID — this caps MeshBlocks per
rank at 2^14 and exists because of MPI tag-space limits. `prolongation.cpp` / `prolong_prims.cpp` and
`flux_correct_{cc,fc}.cpp` handle level boundaries. Physical BCs per module live in `bvals/physics/`.

### EOS and conserved↔primitive inversion

`src/eos/` has one file per (system, EOS) pair: `ideal_hyd`, `isothermal_mhd`, `ideal_srhyd`, `ideal_grmhd`, …
The inversion kernels shared between them are headers (`ideal_c2p_hyd.hpp`, `ideal_c2p_mhd.hpp`).
`eos/primitive-solver/` is the more general (tabulated / piecewise-polytrope) solver used by the dynamical-GR
MHD module. C2P failures and floor applications are counted in `EventCounters` (see `mesh.hpp`) and reported
via the event-log output.

### Relativity

- `src/coordinates/` — `Coordinates` provides cell geometry and, for GR, the stationary metric
  (`cartesian_ks.hpp`); `adm.cpp` holds ADM variables, `excision.cpp` handles excised regions.
- `src/z4c/` — Z4c spacetime evolution (`z4c_calcrhs.cpp`, `z4c_gauge.cpp`, `z4c_update.cpp`), Weyl scalars and
  wave extraction, apparent-horizon dumps, compact-object tracking, CCE, and its own AMR criteria
  (`z4c_amr.cpp`). `z4c_macros.hpp` and `src/athena_tensor.hpp` provide the tensor abstractions.
- `src/dyn_grmhd/` — GRMHD on the dynamical Z4c spacetime; couples through
  `src/tasklist/numerical_relativity.cpp`, which merges Z4c and fluid tasks into one dependency graph.

### Other subsystems

- `src/particles/` — Lagrangian tracers and charged test particles; pushers + tasks + their own boundary path.
- `src/radiation/` — relativistic radiation transport on an angular geodesic grid (`src/geodesic-grid/`).
- `src/diffusion/` — viscosity, conduction, resistivity (the last with optional RKG super-time-stepping and its
  own task lists / CT update).
- `src/srcterms/` — generic source terms, ISM cooling, and the O-U turbulence driver.
- `src/shearing_box/` — shearing-periodic BCs and orbital advection for both cell- and face-centered fields.
- `src/outputs/` — `BaseTypeOutput` derived classes chosen by `file_type` in each `<output*>` block: formatted
  table, history, VTK (mesh and particles), binary, coarsened binary, PDF, Cartesian-grid and spherical-surface
  interpolated output, particle tracking, event log, and restart. `io_wrapper.cpp` abstracts POSIX vs MPI-IO;
  restarts support one-file-per-rank layouts (`rank_%08d/` subdirectories, detected from the path in `main.cpp`).
- `src/parameter_input.cpp` — the `<block>`/param parser; `GetOrAdd*` records defaults so the effective input
  can be written into restarts. Command-line `block/par=value` overrides everything.
- `src/utils/` — interpolation (`lagrange_interpolator`, `chebyshev`, `legendre_roots`), spherical harmonics,
  finite differences, TOV solver, tabulated-EOS readers, `show_config.cpp`.

### Problem generators

Two mutually exclusive mechanisms:

1. **Built in** — `problem/pgen_name` in the input file dispatches through the `if/else` chain in
   `ProblemGenerator::CallProblemGenerator` (`src/pgen/pgen.cpp`) to a method declared in `pgen.hpp` and
   implemented in `src/pgen/tests/`. Adding one means touching all three places.
2. **User file** — `cmake -D PROBLEM=<name>` compiles `src/pgen/<name>.cpp`, defines `USER_PROBLEM_ENABLED`,
   and `CallProblemGenerator` then calls **only** `ProblemGenerator::UserProblem` (the `pgen_name` dispatch is
   compiled out). This is how everything in `src/pgen/*.cpp` (tori, BNS, turbulence, Z4c punctures, …) is used.

A user problem enrolls optional hooks by assigning function pointers on the `ProblemGenerator` object:
`user_bcs_func`, `user_srcs_func`, `user_ref_func`, `user_hist_func`, `pgen_final_func` (called from
`Driver::Finalize()`, e.g. for convergence-error output via `OutputErrors`). The corresponding
`user_bcs`/`user_srcs`/`user_hist` flags are set from the input file, and the code fatals at startup if a flag
is set but the pointer was not enrolled.

## Conventions that matter when editing kernels

- Before a `par_for`, copy everything the lambda needs into local `auto`/`auto &` variables
  (`auto u0_ = u0; auto &mbsize = pmy_pack->pmb->mb_size;`). Never dereference `this` or a host-side pointer
  inside a device lambda — capturing a class member implicitly captures `this` and fails on GPUs.
- Loop bounds come from `pmy_pack->pmesh->mb_indcs` (`is/ie/js/je/ks/ke` for active cells, `cis/cie/…` for the
  coarse buffers used by AMR); dimensionality is handled by branching on `multi_d`/`three_d`, not by changing
  loop nests.
- Variable-index constants (`IDN`, `IM1`, `IEN`, `IBX`, `I00`…`I33`, `IPX`…`IPVZ`) and the algorithm enums
  (`ReconstructionMethod`, `TimeEvolution`, `BoundaryFlag`, `Hydro_RSolver`, `MHD_RSolver`) live in
  `src/athena.hpp` and the module headers — use them instead of literals.
- Reductions (history, dt) use `array_sum::GlobalSum` with `NREDUCTION_VARIABLES` slots; adding history
  variables may require raising that constant.

## Contributing

Feature branches + PRs to `main` (protected, ≥1 review required); direct pushes are rejected. New
functionality is expected to come with a regression test and wiki documentation. See `CONTRIBUTING.md`.
