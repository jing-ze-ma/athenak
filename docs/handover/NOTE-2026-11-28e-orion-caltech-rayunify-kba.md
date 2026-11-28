# NOTE 2026-11-28e orion -> Caltech: rayunify KBA pipelining (several MeshBlocks along r), timing arms

From RAYUNIFY (orion). Queue this after 28d. Idle nodes only.

**KB** = branch **`caltech-rayunify-37f021c12`** = vet-rayunify-n6 (e449d772c) + 2 commits:
- 8989c720c: several MeshBlocks along r are allowed on sp wedges (serialized level by level, exact).
- 37f021c12: KBA pipelining over direction groups closed under the wall maps. Key `rad_vet/vet_gd_kba_groups` (auto | N | 1; 1 = level by level).
- Orion CPU: the sweep is bitwise at any group count, rank count and level count.
- Not supported with several levels: the ray cut (logged off), the compact halo / band ring / IPC (dense halo), and the cubed sphere (FATAL).

Build like the earlier rounds: CUDA 12.9 + g++13, HOPPER90, star_vet. Run the tiny gate tiny_bsg_rel6, pass1 3.720e-11.

## Arms
bsg 27zy deck, cycles 61-120, outputs off. **First read the log line 'closed direction component(s)'** and post it verbatim.

| arm | GPUs | meshblock | keys | note |
|---|---|---|---|---|
| K0 | 4 | 256 128 128 (1 level) | as 27zy | baseline |
| K1 | 4 | 128 128 256 (2 levels x 2 lateral) | `mesh/layout_tile=128 rad_vet/vet_gd_kba_groups=1` | serialized |
| K2 | 4 | as K1 | `... vet_gd_kba_groups=auto` | pipelined |
| K3 | 8 | 128 128 128 (2 x 4) | groups 1 vs auto | |
| K4 (timing only, physics changes) | 4 | as K1 | `rad_vet/vet_gd_wall_interp=false rad_vet/vet_gd_wall_sym=false rad_vet/vet_gd_theta_seam_sym=false`, groups 1 vs auto | run if K2's log says 1 component, i.e. no pipelining with the default maps |

- The 'meshblock' column lists nx1 nx2 nx3. K1/K2 split r in 2 and keep the lateral 128 x 256 split in 2. If that layout is illegal, use 128 x 256 x 128 and say which.
- If the K3 8-GPU arm is not possible, run the scaling arm instead: 2 and 4 ranks along r, lateral one block (meshblock 128 / 64 x 256 x 256, layout_tile 256), groups 1 vs auto.

## Report per arm
1. s/cycle, and the vet_gd fenced line (sweeps, halo+walls, MPI wait).
2. The kba / x1 log lines.
3. Is K1 == K2 bitwise in the CLEAN lines? They must be.

-- RAYUNIFY (orion)
