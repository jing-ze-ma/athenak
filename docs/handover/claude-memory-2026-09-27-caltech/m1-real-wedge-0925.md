---
name: m1-real-wedge-0925
description: User 2026-09-25 (Caltech) -- build a REAL sp wedge (gravity, moving gas, He-envelope stratification) with implicit M1/VET, and verify and time it on GPU like the He box. Settles handover open decision #2. Test-size, not production.
metadata:
  type: project
---

User, 2026-09-25: "also try to have implicit vet working on a real wedge with rad hydro and do the same as you done on box convection".

**Before this**
- The only sp M1 setup was the static T-S4 test (`rad_m1_beam`, `m1_test = sph_atm`): no gravity, gas fixed.
- The 4 Msun presupernova He-star plan ([[he4-presn-global-plan]]) used grey two-stream RT, not M1.

**HANDED TO VIPER 09-25 PDT** (user: Caltech is GPU-bound): task doc docs/handover/TASK-2026-09-26-m1-real-wedge.md + pointer at the top of HANDOVER-2026-09-26.md, pushed c1dc8226. The Caltech agent was stopped before writing code (branch m1-wedge deleted).

**Work in progress (superseded)**
- Branch `m1-wedge`, worktree `/resnick/home/jingze/ATHENAK/wt_m1wedge`, work dir `/resnick/scratch/jingze/m1wedge/`.
- It reuses the nvcc fixes from branch `m1-port`: box_convection + rad_m1 GPU port, worktree `wt_m1port`.

**How to apply**
- Treat it as a development and performance test bed.
- Default design choices (star, radial range, wedge size) are to be confirmed by the user.
- Same gates as the box: CPU bitwise for existing problems, GPU vs CPU values, 1 vs 2 GPUs, restart, H200 lever timing, nsys profile.
