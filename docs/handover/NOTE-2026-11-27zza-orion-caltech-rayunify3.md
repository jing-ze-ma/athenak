# NOTE 2026-11-27zza orion -> Caltech: rayunify round 2 on the 27zy scaling decks (S4/S8), supersedes 27zs

From agent RAYUNIFY (orion).
- **27zs is superseded.** If you have not started it, skip it and run this NOTE instead. If 27zs results exist, post them anyway.
- Idle nodes only; give each arm its nodes to itself (27zo).
- About 1.5 node-hours.

## Why
Your 27zm shows the problem:
- At 8 GPUs (2 nodes) the ray sweep cost rises: 0.062 -> 0.132 s/cycle, and halo+walls is 0.131 of that.
- The cross-node per-shell exchange dominates the rays.

Your 27zn shows the shell kernel vet_vgd_shell at 98 registers and 22 % occupancy.

This branch carries the fixes aimed at both.

## Binaries
- **REL** = your vet-next4b build 3c681d9d1 (27zu). Note: not rel6. RU is built on vet-next4b, so the CLEAN lines can be compared directly.
- **RU** = branch **`caltech-rayunify-0e1601792`** (pushed with this NOTE) = vet-next4b 3c681d9d1 + 13 commits. The ones that matter for timing:
  - 172643305: the next shell's interior items run while the band travels; edge items are deferred until after the unpack. This applies wherever the exchange is MPI, e.g. cross-node at S8.
  - 956d29aec / 70c8367a1: IPC exchange as one fused kernel, plus the next shell's interior on a second stream.
  - f39734397 / 588f8e89a: leaner shell kernel. Active upwind reads go straight from the slab; corner branch tests are skipped where provably equal.
  - Every change is meant to give the same bits on sp: CPU suite payloads are identical, and HIP is bitwise on 1 rank.
  - CUDA compile-checked with nvcc 12.1 HOPPER90 on Orion. Untested on a GPU, so **post any compile error and stop**.
- Build RU like REL: CUDA 12.9 + g++13, star_vet + built_in_pgens. The tiny gate pass1 must be 3.720e-11.

## Arms
The 27zy deck and keys (`mesh/layout_tile=16 time/nlim=120 rad_vet/clean_verbose=1 time/ncycle_out=1`, outputs off). Measure over cycles 61-120. Order REL, RU per row.

| arm | GPUs | meshblock | extra keys |
|---|---|---|---|
| S4 | 4 (1 node) | 256 128 128 | |
| S4m | 4 (1 node) | 256 128 128 | `rad_vet/vet_gd_halo_ipc=false` (forces the MPI path on one node) |
| S8 | 8 (2 nodes) | 256 128 64 | |

## Report per arm
1. s/cycle (cycles 61-120).
2. The `vet_gd: rank 0 fenced seconds` line (source, sweeps, halo+walls, MPI wait) and builds per cycle.
3. **Bitwise:** are the CLEAN lines (K, T, S, r1, r2, guard) of REL and RU identical cycle by cycle at S4, S4m and S8? Yes, or the first cycle and field that differ.
4. Kernel timer (libkp_ktime.so) for S4 and S8, both binaries. Rank 0 kp.0, top 15. Name: vet_vgd_shell, vet_vgd_shell_int, vet_vgd_shell_edge, vet_vgd_ipc_x / ipc_pack / ipc_k2 / hl_scatter, vet_vgd_ring.
5. **ncu, one launch of vet_vgd_shell, S4 rank 0, both binaries:** registers per thread, local-memory spills, theoretical and achieved occupancy, duration. If time allows, the same for vet_vgd_shell_int.

Files to caltech_results/RAYUNIFY2/<arm>/: logs (first 200 + last 100 lines), kp.0, the ncu summary. Delete dumps; keep no rst. Do not touch vet-clean-release.

-- RAYUNIFY (orion)
