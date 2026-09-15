---
name: rkl1-halo-every
description: RKL1 transverse halo exchanged every N substages (rad_tr_halo_every, commit 76b42ca8 rg-box-pcr, default 1): N=2 bitwise at 2 and 4 ranks, operator cost -51 % on the B star (259 -> 222 ms/cycle at 4 ranks); faces_only -9 % but refused with N>1; slabs/edge count irrelevant
metadata:
  type: project
---

Per-swap cost = round trip + launches + fences, independent of message count (slab test 1.3 %). Fix = fewer swaps: refill all nghost layers once per N substages, advance the recurrence on a shrinking ghost skin (Y_{j-1} via the 5-point cross, Y_{j-2} pointwise; N>=3 also exchanges Y_{j-2}); ghost coefficients built nghost deep when N>1. Gate arms m3pcr/prof/C4a-d, C2c vs F4/F2: ALL bitwise (hst, rad_implicit_ang lines, time/dt). ms/cycle at 4 ranks x 2 nodes: F4 259.3, faces_only 236.0, every=3 227.7, every=2 221.8 (best), floor without the operator 186.5. faces_only (skip_x2x3_diag in MeshBoundaryValues, 9 symmetric sites) is refused with N>1 (the skin reads the diagonal ghosts); rad_sts_split + N>1 is fatal. Binary build_hip_prod76b42ca8 (rebuild at HEAD 09-15). USE: hydro/rad_tr_halo_every=2 on the command line at a link break (restart format unchanged, bitwise-equivalent) - planned for the He chain links 2-3 (11710037/38 must be cancelled + resubmitted with the new BIN + override, sbatch snapshots the script) and for any later B-star chain. See [[rkl1-halo-gpu-mpi-race]], [[bstar-prod-cost-profile]], [[sbatch-snapshots-script]].
