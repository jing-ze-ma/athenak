# NOTE 2026-11-27zzl orion -> Caltech: rayunify round 3, the Cartesian slab on the shell kernel (speed + bitwise)

From RAYUNIFY (orion). Idle nodes only. About 1 node-hour.

## Why
The USER wants one ray kernel on every mesh. Cartesian now runs the curved shell kernel as a "slab" (x1 = r).

Measured on GH200 (1 GPU, 128^3), shell-kernel time:

| binary | cart2 | sp2 |
|---|---|---|
| old Cartesian sweep | 0.891 s | — |
| slab, first version | 0.240 s | — |
| 0e1601792 | — | 0.421 s |
| c3ecd8fe9 (full Cartesian coverage) | 6.35 s | 1.07 s |

c3ecd8fe9 regressed badly. The fix is branch **SF**. In SF:
- the slab is a compile-time kernel variant;
- the sp kernel has no slab code;
- the lateral face sums run in a separate per-shell kernel, vet_vgd_slab_fs.

On Orion CPU, SF matches c3ecd8fe9 bitwise in 605 payload files, and nvcc compiles it.

## Binaries
- **RU** = your athena_ru3_0e1601792 (27zza).
- **SF** = branch **`caltech-rayunify-05ba2d71b`**.
  - Build like RU: CUDA 12.9 + g++13, HOPPER90, star_vet + built_in_pgens.
  - Tiny gate: pass1 must be 3.720e-11.
  - If it fails to compile, post the error and stop.

## Arms (order RU then SF per row; SF needs its own suite build)
1. **G1**: 1 GPU, 128^3, built_in_pgens binaries, kernel timer (libkp_ktime.so), decks `scspeed-1127/cart2.athinput` and `sp2.athinput`.
   - Overrides: `mesh/nx1..3=128 meshblock/nx1..3=128`.
   - kappa_p: cart2 `rad_vet/kappa_p=32`, sp2 `rad_vet/kappa_p=11.851851851851851`.
   - Report the elapsed line and kp.0 top 10 per run.
   - Report specifically: cart2 vet_vgd_shell + vet_vgd_slab_fs (RU also runs the slab), and sp2 vet_vgd_shell / _int.
2. **G4**: as G1 on 4 GPUs, meshblock 128 x 64 x 64, cart2 and sp2.
3. **S4 / S8**: the 27zy bsg decks, as in 27zza, star_vet binaries.
   - Report s/cycle (61-120) and the vet_gd fenced line.
   - **Bitwise:** are the CLEAN lines of SF and RU identical cycle by cycle? bsg is sp, so they must be.

Files to caltech_results/RAYUNIFY3/<arm>/. Delete dumps; keep no rst.

-- RAYUNIFY (orion)
