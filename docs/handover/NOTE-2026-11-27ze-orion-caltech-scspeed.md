# NOTE 2026-11-27ze orion -> Caltech: Cartesian vs curved sc ray sweep, cost per build on H200

The USER asks which short-characteristics ray code is faster, the Cartesian sweep
(rays/sc_cart_sweep*.cpp) or the curved-mesh one (rays/sc_sph_sweep.cpp), and whether moving
Cartesian onto the curved path would slow it down. Same cell count, same 162-ray geodesic set
(Cartesian `vet_sc_set = gd`, sp `vet_gd_geo = 4`), same opacity per cell (uniform absorbing gas,
cell dtau 0.25, so no ray cut), fixed set (no rotation). Timing only: no physics check, no outputs.

Binary: your CUDA 12.9 + g++13 **suite** (built_in_pgens) binary of rel6 01deec542 (md5 8cc98fb3,
your NOTE 27v; src = release aa2646b92). No build needed.

Decks (this branch): `docs/handover/scspeed-1127/cart.athinput`, `sp.athinput` (64^3, 6 cycles,
dt_fixed, clean_verbose = 1, no output blocks). Optional per-kernel timer:
`scspeed-1127/kp_ktime.cpp` (a 90-line Kokkos tool, host-only, asks Kokkos for global fencing):
`g++ -O2 -std=c++17 -shared -fPIC -o libkp_ktime.so kp_ktime.cpp`.

Please, on one free H200 node (one rank per GPU, your usual GPU binding):

0. Smoke (1 GPU): each deck with `time/nlim=1`. If either FATALs, post the FATAL lines and stop.
1. **G1, 1 GPU, 128^3, one MeshBlock**, order A B A B:
   - A: `-i cart.athinput mesh/nx1=128 mesh/nx2=128 mesh/nx3=128 meshblock/nx1=128 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=32`
   - B: `-i sp.athinput` + the same 7 overrides
2. **G4, 4 GPUs, 128^3, 4 MeshBlocks split laterally** (x1 whole: Cartesian banded plane exchange
   every launch vs curved per-shell halo), order A B A B: as G1 with
   `meshblock/nx1=128 meshblock/nx2=64 meshblock/nx3=64`.
3. **G4x, 4 GPUs, 128^3, 2 MeshBlocks along x1** (`meshblock/nx1=64 meshblock/nx2=128 meshblock/nx3=64`),
   A then B (both sweeps then pass the x1 block face in order).
4. **G4L, 4 GPUs, 256^3** (`mesh/nx1..3=256 meshblock/nx1=256 meshblock/nx2=128 meshblock/nx3=128 rad_vet/kappa_p=64`),
   A then B. Skip if sp runs out of memory (post the error).
5. If you built the tool: rerun G1 A and B and G4 A and B once with
   `KOKKOS_TOOLS_LIBS=<path>/libkp_ktime.so KP_OUT=<rundir>/kp` (each rank writes kp.<rank>);
   post the top 25 lines of kp.0 of each.

Per run please post: rc, the per-cycle `elapsed`/wall lines (s/cycle over cycles 2..6), the
end lines `cpu time used` and `zone-cycles/cpu_second`, and these timer lines:
- Cartesian: `<rad_vet> vet_sc: formal solutions=... SC seconds=... per call ...`
- sp: `<rad_vet> vet_gd: N builds checked ...`, `<rad_vet> vet_gd: rank 0 fenced seconds: source ..., sweeps ... (of which per-shell halo+walls ..., of which MPI wait ...)`, `<rad_vet> vet_gd halo (rank 0): exchanges ...`
- any FATAL / CLEAN WARNING line, and the `vet_sc_set = gd: 162 rays` / vet_gd set log line.
Delete any rst/bin the runs write. Post a short NOTE back (caltech -> orion, "scspeed").
Do not touch vet-clean-release.

-- orion (SCSPEED, for the coordinator)
