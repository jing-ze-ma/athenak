---
name: caltech-h200-validation-0925
description: caltech-port on H200 (2026-09-25) -- GPU deterministic, restart bitwise, 2-GPU bitwise; dhj timing H200 ~1.9-2.2x faster than viper MI300A; the three porting bug classes found
metadata:
  type: project
---

Branch caltech-port (HEAD dfae3fef, merge on HOLD: [[caltech-port-hold-merge]]). Binary athena_p4_gpu (from 874e078f), WASP-121b 1x base arm, production keys.

**Validated 2026-09-25:**
- CPU gate bitwise vs the unported binary (c0926 vs p4).
- GPU run to run: bitwise.
- GPU restart (20 + 20 vs 40 cycles): rst data bitwise. Cosmetic only: the hst writes a duplicate row at the restart time and skips one output.
- 2 H200 vs 1 H200 (srun --mpi=pmix, HPC-X): rst bitwise. Only hst reductions differ, at 1e-16, except the near-cancelling 3-mom column.

**H200 ms/cycle** (cycles 52-196, xstep 8):

| grid | 1 GPU | 2 GPUs |
|---|---|---|
| nx1 76 | 22.4 | 13.3 |
| nx1 256 8-coef | 59.8 | 34.5 |

Viper, 2 MI300A, nx1 256: 75.3 ms with xstep 2, about 67 ms at xstep 8. So the H200 is about 1.9x faster. Going from 1 to 2 GPUs gives about 1.7x.

**Porting bug classes (APU hid them on viper):**
1. nvcc: generic/auto device lambdas.
2. Host code reading device Views: EOSTable now keeps host mirrors.
3. DualViews filled on the host but never synced to the device. mb_panel gave SILENTLY wrong ICs on 8/24 blocks.

Unchecked, off the dhj path: the geodesic-grid, cart_grid and spherical_surface iindcs/iwghts DualViews.

**How to apply:** for any new setup on the H200, compare GPU vs CPU hst/bin at t=0 and after a few cycles. Finishing without errors is not proof.
