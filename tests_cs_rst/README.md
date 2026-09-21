# tests_cs_rst — the general-EOS **CUBED-SPHERE** restart gate

Third and last instalment of the tabulated-EOS restart work: `6ae9ccbe` (hydro,
Cartesian), `0eb7e2c4` (MHD, Cartesian, gate in `tests_mhd_rst/`), and this directory,
which closes the hole those two left open and named:

> *Still not bitwise: the cubed sphere (`GnomonicEquiangleRaiseVel{,MHD}` re-solves the
> caches after c2p, hydro too)* — `0eb7e2c4`

## The hole

On every other grid the conserved-to-primitive conversion is `EquationOfState::ConsToPrim`
and nothing else, so freezing the thermodynamic cache inside it (`ConsToPrimFrozen`, the
two commits above) is enough for a restart to be a bitwise continuation.

On the **cubed sphere** the inversion is split in two, and the half that owns the general-
EOS cache is the second one:

| file:line | what it does |
| --- | --- |
| `src/hydro/hydro_tasks.cpp:966` (now `:977`) | calls `Coordinates::GnomonicEquiangleRaiseVel` right after `ConsToPrim` |
| `src/mhd/mhd_tasks.cpp:735` (now `:746`) | calls `Coordinates::GnomonicEquiangleRaiseVelMHD`, likewise |
| `src/coordinates/coordinates.cpp:803` / `:1012` | those two routines |
| `src/coordinates/gnomonic_raisevel.hpp:137` (`if (gen_)`) | **re-solves** `wtemp` and `wder` from the metric-corrected internal energy, warm started on the cached temperature |
| `src/coordinates/gnomonic_raisevel_mhd.hpp:94` (`if (gen_)`) | the same for MHD |

So on that grid the frozen `ConsToPrim` bought nothing: the cache the restart file had just
restored was overwritten a moment later by a fresh table inversion — started from a
*different* warm start than the straight run used (in the straight run the guess is the
temperature the ordinary `ConsToPrim` solved microseconds earlier from the **orthonormal**
internal energy; after a restart it is the converged value in the file). The root find
stops on a step size, so its answer depends on the guess it started from, and the two runs
enter the first post-restart cycle with different `p` and `Gamma_1`.

Two smaller defects of the same call sit next to it and are closed by the same change:

* those routines also **re-apply the deferred floors and the velocity ceiling**
  (`EOS_Data::defer_cons_floors`) to a conserved state that the file already holds in
  floored form;
* the `dfloor_keep_velocity` / `dfloor_keep_temperature` halves read `Hydro::dfl_fv`,
  which `GeneralHydro::ConsToPrimFrozen` does **not** write — so after a restart they ran
  against the *previous* cycle's density-floor record.

## What changed

| file | what |
| --- | --- |
| `src/coordinates/gnomonic_raisevel_frozen.cpp` (new) | `Coordinates::GnomonicEquiangleRaiseVel{,MHD}Frozen()`: the **algebraic half only** — metric velocity raise, metric kinetic energy, (MHD) cell-centred field in the orthonormal frame and its magnetic energy, `w0(IEN) = etot - ekin [- emag]` — in the same expressions and the same order as `GnomonicRaiseVelFloors` / `GnomonicRaiseVelMHDFloors`. No floor, no ceiling, no inversion; `u0`, `wtemp` and `wder` are not written. Its own translation unit for the reason given in `general_hyd_floors.cpp` (a second kernel in the same TU moves the default kernel's arithmetic by an ULP, and `coordinates.cpp` holds the default cubed-sphere path). |
| `src/coordinates/coordinates.hpp` | the two declarations |
| `src/CMakeLists.txt` | the new file |
| `src/hydro/hydro_tasks.cpp`, `src/mhd/mhd_tasks.cpp` | `const bool frozen = c2p_freeze_derived;` is taken **before** the `ConsToPrim` branch clears the flag, and the cubed-sphere call below takes the frozen variant on that one call |

**Order.** The state in the restart file was written *after* the straight run's own raise,
floors, ceiling and cache solve, so the restarted first call has to reproduce *that* state,
not the pre-raise one. That is exactly what redoing the algebra and leaving the cache alone
does: `w0` comes back out of the conserved variables the file restored, and `wtemp`/`wder`
are the numbers the straight run's own inversion produced.

No flag, no input parameter: the lifetime is the existing `Hydro::c2p_freeze_derived` /
`MHD::c2p_freeze_derived`, one call after a restart.

## The gate

```
./gate_cs.sh <athena binary> <input file> <tag> [meshblock nx2=nx3] [mpi ranks]
```

Same shape as `tests_mhd_rst/gate.sh`, and it reuses that directory's `cmpdumps.py`:
run **A** is 20 cycles uninterrupted, run **B** is 10 cycles + restart + 10, and the
per-cycle binary dumps are matched by the cycle number in their own header and compared
bit for bit.

Problem: a **radial shock tube on the cubed sphere** with the tabulated general EOS —
`cs_hyd_geneos.athinput` and `cs_mhd_geneos.athinput`, built-in pgen `shock_tube`
(so the gate binary needs no `-D PROBLEM=` and does not depend on any user problem
generator), 6 × 8 × 16 × 16 cells, the EOS block copied from `tests_mhd_rst`, so no
external table file either. Both states carry a small tangential velocity, so the metric
cross term is nonzero everywhere off the panel midlines; the MHD twin adds the
divergence-free azimuthal field the pgen lays down on this grid (`<problem>/bazi`), which
is tangential and therefore really needs the orthonormal rotation of `bcc`. The two states
straddle H2 dissociation and H ionisation (T ≈ 1.2e3 → 1e4 K), the part of the table where
the inversion is hardest.

`meshblock/nx2 = nx3 = 16` is one MeshBlock per panel (6), `= 8` is 2x2 per panel (24).

### Results (viper login node, Serial Kokkos, gcc 14)

| gate | `athena_ref` (tree at `0eb7e2c4`) | `athena_new` (after this change) |
| --- | --- | --- |
| cs **MHD**, 1 block/panel | **FAIL** at cycle 11, the first post-restart cycle: worst `ener` abs 8.53e-1, rel 3.24e-5; by cycle 20 abs 1.61e+3, rel 6.3e-3 | **PASS**, 21/21 bitwise |
| cs **MHD**, 2x2 blocks/panel | **FAIL** at cycle 11: `ener` abs 7.38e-1, rel 2.06e-5 | **PASS**, 21/21 |
| cs **hydro**, 1 block/panel | PASS (see below) | **PASS**, 21/21 |
| cs **hydro**, 2x2 blocks/panel | PASS (see below) | **PASS**, 21/21 |
| cs MHD, **ideal** EOS (control) | PASS | — |
| Cartesian MHD `tests_mhd_rst/gate.sh`, 1 and 4 blocks | PASS | **PASS** (no regression) |
| Cartesian hydro `tests_mhd_rst/gate_hyd.sh` | PASS | **PASS** (no regression) |
| cs MHD, 2x2 blocks/panel (24 MB) on **2 MPI ranks** | not run | **PASS**, 21/21 |
| cs hydro, 2x2 blocks/panel on **2 MPI ranks** | not run | **PASS**, 21/21 |

A straight 20-cycle run of both cs inputs is **bitwise identical** between `athena_ref` and
`athena_new`: the change moves no default answer.

**Where the MHD failure lived.** With the reference binary the conserved state *and* `bcc`
dumped immediately after `Driver::Initialize` (the restarted run's first dump, cycle 10)
are bitwise identical to the straight run's cycle 10, and the first post-restart timestep
`dt` is identical as well. The divergence therefore enters cycle 11 through the only state
the first `Fluxes` call reads that the dumps do not show: the `wtemp`/`wder` the
raise-velocity kernel had just re-solved. Freezing that call closes it.

## What is honestly still open

* **The hydro half of this hole never failed on this problem.** `GnomonicEquiangleRaiseVel`
  re-solves the cache exactly as the MHD one does, but on both cs hydro gates the
  reference binary already produced 21/21 bitwise cycles: with the warm starts this
  problem generates, the table root find returned its own converged root bit for bit. That
  is a property of the states, not a guarantee — the MHD twin, on the same grid and the
  same table, diverged at 3e-5 in the first cycle. The change removes the dependence for
  both.
* **Floored cells are reproduced through a round trip.** The frozen kernels rebuild
  `w0(IEN)` as `etot - ekin [- emag]`. Where a floor fired in the straight run, `etot` was
  itself written as `eint + ekin [+ emag]`, so the round trip can differ from the stored
  `eint` by an ULP in those cells (nothing else in the gate's state moves). No cell in
  these gates floors; a restart of a run that floors heavily in the same cycle may still
  show last-bit motion there.
* **FOFC** (`Hydro::FOFC` / `MHD::FOFC`) runs the same per-cell bodies over a trial state
  and was not touched; it is a startup fatal on this grid anyway.
* Everything `0eb7e2c4` listed as open stays open: the mode-3 RT Newton warm start with
  `problem/rt_impl_warm = 1`, and the pgen-side state of `deep_hot_jupiter_rt` (not
  audited — that file is being edited elsewhere).

## MPI

`gate_cs.sh <bin> <input> <tag> 8 2` runs the 24-MeshBlock case on 2 ranks; the binary is
`athena_new_mpi`, built with `-D Athena_ENABLE_MPI=ON` after
`module load gcc/14 openmpi/5.0`. Both cs gates pass there, 21/21 bitwise.

## Housekeeping

The gate scripts delete their binary dumps and restart files after the comparison (inode
quota); the `.cmp` tables and the run logs stay. The build directories
`build_csrst_{ref,new,mpi}` were removed after the runs; `athena_ref` and `athena_new` are
kept here so the table above can be reproduced without rebuilding.
