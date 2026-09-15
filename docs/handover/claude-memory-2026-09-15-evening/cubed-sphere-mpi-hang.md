---
name: cubed-sphere-mpi-hang
description: Cubed sphere on more than one MPI rank -- FIXED, 1c2e29d6; it hung in ClearRecv on an unmatched receive, same cube-vertex root cause as the seam EMF bug. 1/2/4 ranks now bitwise identical
metadata:
  type: project
---

**FIXED 2026-08-28 (sixteenth session), committed 1c2e29d6.** Found while regression-testing
[[cubed-sphere-seam-emf]]; it had been broken for as long as panel boundaries have existed.

## The symptom

`mpirun -np 2` on any cubed-sphere input hung **after the mesh was built and before cycle 0
completed** -- pure HYDRO as well as MHD, at HEAD 7b421818 as well as on the working tree,
so it was neither new nor in the EMF path. `gdb -p` on both ranks showed them parked in
`MeshBoundaryValues::ClearRecv` inside `MPI_Wait`. **Attaching gdb to the hung ranks is what
turned this from "MPI is broken somehow" into a one-line diagnosis** -- do that first.

## The cause: the SAME cube vertex, one more consequence

A reciprocity audit of the whole neighbour table (for each `(m,n)` with `gid >= 0`, check
`nghbr(nghbr(m,n).gid, dest).gid == gid(m)`) found the **only** non-reciprocal slots in the
entire table are the x2x3 edges at cube vertices -- 8 at 1 block/panel, 24 at 4:

| | reciprocal | NON-reciprocal |
|---|---|---|
| x2face / x3face | 96 | 0 |
| x2x3 edge | 72 | **24 (all cube vertices)** |
| everything else | 0 | 0 |

(The empty rows are not a gap: the radial direction has physical boundaries at both ends,
so a block has no x1-direction neighbour and the x1face / x1x2 / x3x1 / corner slots never
exist at all.)

A send is tagged with the **RECEIVER's** `(lid, dest)`, so a non-reciprocal slot is a posted
receive nobody satisfies -- and its own send collides with the legitimate sender for the
slot it targets. One unmatched receive hangs `ClearRecv`.

## The fix

The same one 2646f826 applied to the EMF: **skip the x2x3 edge exchange at a cube vertex,
send and receive together**, via `IsCubeVertexCorner` now promoted to `bvals.hpp`. Applied
in `PackAndSendCC/FC`, `RecvAndUnpackCC/FC`, their MPI send and `MPI_Test` loops, and
`InitRecv` / `ClearRecv` / `ClearSend`. Nothing is lost -- `FillPanelCornersCC/FC` already
overwrite exactly that ng x ng corner block from the two flanking face halos, **which is why
the serial answer is bit-identical**. `ClearFluxRecv`/`ClearFluxSend` needed no change: they
already guard on `MPI_REQUEST_NULL`, which is also what kept 2646f826's skipped flux
receives safe.

Verified: cs MHD and cs hydro **bitwise identical at np = 1, 2, 4**, and np=1 identical to
the serial binary; `linear_wave_mhd` identical to HEAD serially and at matched rank count
(vtk included).

## Reusable

**Run the reciprocity audit before anything else for a boundary bug.** `CS NEIGHBOUR
RECIPROCITY AUDIT` in `cs_test.cpp` (iprob=8, single rank, where every block is local)
reduced "the run hangs somewhere in MPI" to "these 24 slots and only these" in one
measurement, **needing no MPI build at all**.

## MPI build recipe (Viper login node, CPU)

```bash
source /etc/profile.d/modules.sh
module purge && module load gcc/14 openmpi/5.0 cmake/4.0
cmake -S /viper/u2/jinma/ATHENAK/athenak -B <build> -D PROBLEM=cs_test \
      -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx
mpirun -np 2 --oversubscribe <build>/src/athena -i ...
```

`openmpi/5.0` is **hierarchical** -- invisible until `gcc/14` is loaded -- and
`-D CMAKE_CXX_COMPILER=mpicxx` is required or cmake picks `/usr/bin/c++` and fails the MPI
probe. See [[viper-hip-build-recipe]]. Builds live in session a76ff032's scratchpad:
`bmpi/` (cs_test) and `bmpi_lw/` (built_in_pgens).
