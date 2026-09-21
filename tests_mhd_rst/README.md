# tests_mhd_rst — the general-EOS **MHD** restart gate

Companion of the hydro work in `6ae9ccbe` ("box-restart-bitwise"), whose commit message
already said *"the same wder loss exists for general-EOS MHD"*.  This directory is the
gate for closing it.

## What was wrong

Under a tabulated (general) EOS the `e -> T` inversion is **not idempotent**: handed back
its own converged answer it takes another step and lands a few ULP away.  A straight run
never performs that inversion at the start of the cycle it is resuming — it carries
`wtemp` and `wder` (the derived `p`, `Gamma_1`, channels `IDPR`/`IDG1`) from the last
stage of the previous cycle — while `Driver::Initialize` runs one extra `ConsToPrim` after
a restart.  Restart files carried **MHD `wtemp`** already, but never **MHD `wder`**, and
the first `Fluxes` call of a restarted cycle reconstructs `wder`, not `w0(IEN)`.  So every
restart link boundary of a chained general-EOS MHD production (the deep hot Jupiter runs)
was a small discontinuity that then grew chaotically.

## What changed

| file | what |
| --- | --- |
| `src/eos/general_mhd_frozen.cpp` (new) | `GeneralMHD::ConsToPrimFrozen()`: the **algebraic half only** of `ConsToPrim` — d, v, bcc from the face field, `w0(IEN) = e - e_k - e_m` — in the same expressions and the same order as `SingleC2P_GeneralMHD`.  `wtemp`/`wder` are read by later kernels, never written here.  Its own translation unit for the reason given in `general_hyd_floors.cpp` (a second kernel in the same TU moves the default kernel's arithmetic by an ULP). |
| `src/eos/eos.hpp` | declaration on `GeneralMHD` |
| `src/CMakeLists.txt` | the new file |
| `src/mhd/mhd.hpp` | `MHD::c2p_freeze_derived`, set by the restart reader, cleared by the one call that honours it |
| `src/mhd/mhd_tasks.cpp` | `MHD::ConToPrim()` takes the frozen path once after a restart |
| `src/outputs/outputs.hpp`, `src/outputs/restart.cpp` | `outarray_wdpm` / `outarray_wdgm`: MHD `wder` written to the restart file, behind the hydro `wder` and in front of the mode-3 warm-start slabs |
| `src/pgen/pgen.cpp` | tail-size detection and the reader; sets `pmhd->c2p_freeze_derived` |

Face-centred `b0` and `bcc0` handling is untouched (30d21859 / 4cfdc329 stay as they are);
the frozen kernel rebuilds `bcc` from the restored face field with exactly the branch and
the weights `ConsToPrim` uses.

**Old restart files stay readable.**  Four tail lengths are accepted: both caches; the
hydro derived cache but not the MHD one (only a hydro+MHD ion-neutral run can see that);
`wtemp` only — which is what *every* general-EOS MHD restart file written so far is; and
neither.  The last three print a `### WARNING` saying the restart is not bitwise.
Verified: a file written by the reference binary loads under the new one with that
warning and runs.

## The gate

`./gate.sh <athena binary> <tag> [meshblock nx1] [mpi ranks]`

Problem: `inputs/tests/linear_wave_mhd_geneos_table.athinput` (built-in pgen
`linear_wave`, tabulated general EOS, hlld, 64 cells, physical cgs), copied here as
`lw_mhd_geneos.athinput` with `nlim = 20` and hst / rst / bin output blocks.  No
`-D PROBLEM=` build is needed, and no external EOS table file.

* run **A**: 20 cycles uninterrupted
* run **B**: 10 cycles, restart from the cycle-10 restart file, run to cycle 20
* `cmpdumps.py` matches the per-cycle binary dumps by the cycle number in their own header
  and compares the **data** bit for bit (the text header is not compared: it embeds the
  parameter dump, whose command line and `last_time` strings legitimately differ).

`gate_hyd.sh` is the same thing with `linear_wave_hydro_geneos_table.athinput`: the
no-regression check that the hydro half stays closed.

### Results (viper login node, Serial Kokkos, gcc 14)

| gate | reference binary (`athena_ref`, unmodified tree) | after the change (`athena_new`) |
| --- | --- | --- |
| MHD, 1 MeshBlock | **FAIL** at cycle 11 (the first post-restart cycle): `mom1` abs 1.12e-15, rel 7.5e-5; grows to abs 3.19e-15 by cycle 20 | **PASS**, 21/21 cycles bitwise |
| MHD, 4 MeshBlocks (`meshblock/nx1=16`) | **FAIL**, identical numbers | **PASS**, 21/21 |
| MHD, 4 MeshBlocks on 2 MPI ranks | not run | **PASS**, 21/21 (`build_mhdrst_mpi`, openmpi/5.0) |
| hydro (`gate_hyd.sh`) | PASS | PASS (no regression) |

A straight 20-cycle run is **bitwise identical** between `athena_ref` and `athena_new`:
the change moves no default answer.

## What is still NOT bitwise

* **Cubed sphere.**  `Coordinates::GnomonicEquiangleRaiseVelMHD`
  (`src/coordinates/coordinates.cpp:1012`, the `if (gen_)` block at
  `src/coordinates/coordinates.cpp:1109`) re-solves `wder`/`wtemp` from the corrected
  internal energy *after* `ConsToPrim`, so on a cubed-sphere grid the frozen cache is
  overwritten by a fresh inversion and the restart is still not bitwise.  The hydro
  counterpart (`GnomonicEquiangleRaiseVel`, called from `src/hydro/hydro_tasks.cpp:966`;
  MHD at `src/mhd/mhd_tasks.cpp:735`) has exactly the same hole and 6ae9ccbe did not close
  it either.  Closing it means threading a "freeze" flag into those two kernels.
  Spherical-polar and Cartesian grids are unaffected.
* **Mode-3 RT Newton warm start** with `problem/rt_impl_warm = 1` (the production
  setting): per-cell history partly outside the restart file — see the 6ae9ccbe message.
* **pgen-side state of `deep_hot_jupiter_rt`** was not audited here (that file is being
  edited by another agent).  The known class of defect is the one 6ae9ccbe fixed in
  `box_convection.cpp`: dump clocks re-armed at the restart time, and any cached
  background rebuilt on a cycle-count gate (`wb_cache_every`, `MHD::wb_cache_built`,
  `src/mhd/mhd.hpp:257`).

## Housekeeping

Binary dumps and restart files are deleted by the gate scripts after the comparison
(inode quota).  The build directories `build_mhdrst_{ref,new,mpi}` were removed after the
runs; `athena_ref` / `athena_new` / `athena_new_mpi` are kept here so the table above can
be reproduced without rebuilding.
