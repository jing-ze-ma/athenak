# NOTE 2026-11-28f orion -> Caltech: 28e amended. KBA arms withdrawn, two arms kept for the cost of several MeshBlocks along r

From RAYUNIFY (orion). The USER has dropped direction-group KBA.

**Withdrawn:** arms K1, K2 and K4 of NOTE 28e. Do not run them. Do not build caltech-rayunify-37f021c12.

**Kept:** K0 and K3, run on branch **`caltech-rayunify-8989c720c`**. That branch allows several MeshBlocks along r on sp wedges; the sweep runs level by level and is exact. It measures what several levels along r cost.

**Build:** CUDA 12.9 + g++13, HOPPER90, star_vet. Tiny gate tiny_bsg_rel6, pass1 3.720e-11.

**Arms:** bsg 27zy deck, cycles 61-120, outputs off.

| arm | GPUs | meshblock (nx1 nx2 nx3) | layout_tile | levels along r |
|---|---|---|---|---|
| K0 | 4 | 256 128 128 | 128 | 1 |
| K3 | 8 | 128 128 128 | 128 | 2 (x 4 lateral) |

If K3 is not legal, run 4 GPUs at meshblock 128 128 256 and say which.

**Report per arm:**
- s/cycle;
- the vet_gd fenced line;
- the x1 / ray-cut log lines (with several levels the ray cut is off).

-- RAYUNIFY (orion)
