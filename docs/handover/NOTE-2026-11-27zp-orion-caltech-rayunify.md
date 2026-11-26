# NOTE 2026-11-27zp orion -> Caltech: rayunify round 1, ray-sweep speed + bitwise checks on H200 (1/4/8 GPUs)

From agent RAYUNIFY (orion). USER goal: one ray kernel on every mesh, faster; per-shell exchange latency is a first-class target. **Idle test nodes only; never stop or share a node with production.** About 1.5 node-hours.

## Binaries
- **REL** = your athena_rel6_01deec542 (aa2646b92 src; built_in_pgens and star_vet).
- **RU** = branch **`caltech-rayunify-97d437d68`** (pushed with this NOTE) = aa2646b92 + 4 commits:
  - b33edb8ba: Cartesian per-cell sums in 7 groups (one thread per cell and group), and a ray-parallel banded (several-MeshBlock) Cartesian sweep. Same bits by construction.
  - 3d1cade0a: sp IPC exchange, band-ring ops flushed before the hook, so the fused kernel `vet_vgd_ipc_x` fires (it never did with the ring on). CUDA only, same bits.
  - 7d201ea25: sp over MPI (no IPC). The next shell's interior items run behind this shell's pack; items whose stencil reaches the band are deferred to an edge list after the unpack. Exact, same bits.
  - 97d437d68: sp with IPC (CUDA). The next shell's interior runs on the second stream vgd_es2, concurrent with `vet_vgd_ipc_x`, joined by events. Same bits.
  - Orion CPU gate of b33edb8ba: full 61/61, and all 2992 suite payload files are byte-identical to aa2646b92. Gate of the same tree (6e2237287, differs from 97d437d68 only in one CUDA lambda qualifier) is running; RU is untested on CUDA: post the compile error if any and stop.
  - Build RU like REL: CUDA 12.9 + g++13, HOPPER90, never CUDA 12.2. Build both `built_in_pgens` (for the scspeed decks) and `star_vet` (for bsg). Run the tiny gate on star_vet: pass1 must be 3.720e-11. Post both md5 and the prov.

## Decks
`scspeed-1127/cart2.athinput` and `sp2.athinput` (SCSPEED's production-wedge decks: periodic theta/phi, GL 7 x 18 staggered = 126 rays). Non-uniform data for the bitwise arms: add
`HB = problem/vet_test=hot_blob problem/blob_rho=1.0 problem/blob_temp=2.0 problem/bg_rho=0.1 problem/bg_temp=1.0`, plus
- cart2: `problem/blob_x1=0.5 problem/blob_x2=0.5 problem/blob_x3=0.5 problem/blob_a1=0.25 problem/blob_a2=0.25 problem/blob_a3=0.25`
- sp2: `problem/blob_x1=2.8 problem/blob_x2=1.5707963 problem/blob_x3=0.5 problem/blob_a1=0.6 problem/blob_a2=0.25 problem/blob_a3=0.25`
Dumps for the bitwise arms: `output1/file_type=bin output1/variable=vet output1/bin_double=true output1/dcycle=6 time/nlim=6`.

## 1. Bitwise REL vs RU (do first; stop and post if any arm differs)
64^3 (`mesh/nx1..3=64`), HB on, dumps on. Compare the md5 of all .bin files (header excluded or included, the same for both).

| arm | deck | GPUs | meshblock | path exercised |
|---|---|---|---|---|
| B1c | cart2 | 1 | 64^3 | grouped sums |
| B4c | cart2 | 4 | 64 x 32 x 32 | banded ray-parallel sweep + grouped sums |
| B1s | sp2 | 1 | 64^3 | (unchanged) |
| B4s | sp2 | 4 | 64 x 32 x 32 | IPC: fused exchange + 2-stream overlap |
| B4sm | sp2 | 4 | 64 x 32 x 32, `rad_vet/vet_gd_halo_ipc=false` | MPI interior/edge split |

## 2. Speed REL vs RU (no dumps, uniform data as in 27zi, order A B A B)
As 27zi H1/H4/H4L, plus 8 GPUs (2 nodes). kappa_p: cart2 32 / sp2 11.851851851851851 at 128^3; 64 / 23.703703703703702 at 256^3.
- H1: 1 GPU, 128^3, cart2 and sp2.
- H4: 4 GPUs, 128^3, meshblock 128 x 64 x 64, cart2 and sp2; sp2 also with `vet_gd_halo_ipc=false` (H4m).
- H4L: 4 GPUs, 256^3, meshblock 256 x 128 x 128, cart2 and sp2.
- H8L: 8 GPUs (2 nodes), 256^3, meshblock 256 x 128 x 64, cart2 and sp2.
- Kernel timer (libkp_ktime.so) on H4 sp2, H4m sp2 and H4 cart2, both binaries: top 25 lines of kp.0.

Report per arm:
- s/cycle (cycles 2-6);
- the `vet_sc ... per call` line (Cartesian) and the `vet_gd: rank 0 fenced` + `vet_gd halo` lines (sp).

Report from kp:
- shell kernel time per build (vet_vgd_shell + vet_vgd_shell_int + vet_vgd_shell_edge);
- halo kernels (vet_vgd_ipc_x / ipc_pack / ipc_k2 / hl_scatter / ring);
- Cartesian vet_vet_sweep_ray, vet_vet_sweep_sum, vet_vet_mb_ray, vet_vet_mb_sum, vet_vet_mb_sweep.

Note: the kp tool fences every kernel, so it serialises the 2-stream overlap. Use s/cycle and the fenced-seconds line (no kp) for the overlap gain.

## 3. Production bsg restart (cost per cycle; bitwise)
Use the bsg rst of 27zc (cycle 29657) and its production keys, named. Run 20 cycles per arm, outputs off except the final hydro_w + user.hst, on the star_vet binaries, in the order REL, RU, REL, RU:
- **P4**: 4 GPUs (1 node), the production layout;
- **P8**: 8 GPUs (2 nodes), the same total mesh, the nearest legal layout (say which).

Report:
- s/cycle (mean over cycles 5-20) and the vet_gd fenced line;
- the per-shell overhead = (fenced sweeps - shell kernels - halo kernels) / shell steps, if you have kp for one extra arm per binary at P4;
- **bitwise**: the final hydro_w and user.hst of REL and RU must be byte-identical at P4 and at P8. Answer yes / no.

Files to caltech_results/RAYUNIFY/<arm>/ (logs: first 200 + last 100 lines, kp.0). Delete run dumps / rst after. Do not touch vet-clean-release.

-- RAYUNIFY (orion)
