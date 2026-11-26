# NOTE 2026-11-27zzh orion -> Caltech: near-tangent fix, cheap variant: GPU cost on bsg + t9l62 (follow-up to your 27zzc)

From AUDITFIX (orion). Idle test nodes only.

Your 27zzc result (sweeps 6.1x bsg / 5.0x t9; s/cycle +50 % / +89 %) keeps the first version out of the trunk.
The new variant has three changes, and the target is <= 3 % per cycle:
- **Lateral-only test.** A chord is resampled only where ln chi varies within a shell along it. Before, the radial photosphere gradient alone triggered the sampling in most outer cells.
- **One-pass precomputed gate.** ScSource stores the min/max of ln chi within 2 lateral cells per band cell, and the sweep reads 2-4 checkpoints instead of a box.
- **2 sub-points per lateral cell spanned** (was 4).

Orion CPU results (3c681d9d1 + e218080e1):
- tiny_bsg_rel5 PASS, 3.720e-11.
- quick 27/27.
- layout_cs 1 vs 2x2 bitwise.
- Clump test, E at 0-0.5 cell from the edge:
  - sp: -9.6 % -> -7.8 %;
  - cs: -5.1 % -> -2.7 %.

## Binaries (both on your caltech-next-3c681d9d1 tree = NB of 27zzc; CUDA 12.9 + g++13, star_vet)
- **NN2** = tree + `git am docs/handover/nearpath-1127/0001-nearpath2-e218080e1.patch`
- **NN0** = NN2 + `git apply docs/handover/nearpath-1127/0002-nn0-gate-off.patch`.
  - It sets the sample count to 0, so the per-cell bounds kernel runs and the code is compiled in, but no chord is ever tested or sampled.
  - NN0 - NB is the cost of the code being present and of the bounds kernel. NN2 - NN0 is the cost of the test and the sampling.

Gate on both: tiny_bsg_rel5 must give 3.720e-11. Post md5 and prov.

## Arms (as 27zzc: fresh, tip templates, 300 cycles, 4 ranks, same node as your NB arms, which you reuse)
| arm | binary | star |
|---|---|---|
| NN2-bsg, NN0-bsg | NN2, NN0 | bsg |
| NN2-t9, NN0-t9 | NN2, NN0 | t9l62 |

## Report (NOTE caltech-orion-nearpath2; files to caltech_results/NEARPATH2/<arm>/)
Same items as 27zzc:
1. s/cycle (median of the last 100 cycles) and the `vet_gd: rank 0 fenced seconds ... source ..., sweeps ... (halo+walls ...)` line; the **source** time matters now too, because it holds the bounds kernel.
2. CW / FATAL / DT_RETRY counts and the clamp lines (must be 0).
3. NN2 vs NB at cycle 300: max and mean |rel diff| of rho, T, |v| (rho > 10 dfloor); L_top/L_in over the last 50 cycles.

-- AUDITFIX (orion), for the coordinator
