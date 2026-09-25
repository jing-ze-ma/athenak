---
name: m1-hesdirk2-cfl-recommendation
description: "For new implicit M1 (vet_sc/Eddington) runs with hesdirk2, decide the CFL: 0.6 (2x more accurate than be@0.3, 1.19x cheaper) or 0.9 (be@0.3 accuracy, 1.57x cheaper); be blows up at 0.9"
metadata:
  node_type: memory
  type: project
  originSessionId: 05818bc0-259d-4b2f-b317-f82d1e083443
  modified: 2026-09-24T17:32:25.024Z
---

Raise to the user when setting up ANY new simulation with implicit M1 / VET (He boxes, He slab, sp wedge once
hesdirk2 is supported there): choose `time/cfl_number` for time_scheme = hesdirk2 (the default on Cartesian implicit
M1 since m1-defaults2, 09-24; be stays default on sp / vet_col / cs, where hesdirk2 is refused).

Options (validation /viper/ptmp2/jinma/h2val_0924/README.md, 3-D He box vet_sc, GPU, per simulated second):
- cfl 0.6: 2x MORE accurate than be at cfl 0.3, 1.19x cheaper (0.232 vs 0.277 s/s).
- cfl 0.9: same accuracy as be at 0.3 (within 10 % in rho/eint/E), 1.57x cheaper (0.177 s/s).
- be at cfl 0.9 is NOT an option: He slab blows up (KE x150-3000); hesdirk2 bounded to cfl 1.2.
- Numbers predate m1-h2fast (hesdirk2 per step 1.85x -> 1.35x be), so the hesdirk2 margins are now larger.

**Why:** the user asked (09-24) to remember this so the CFL is decided consciously at the next setup, not left at 0.3.

**How to apply:** at setup, present the two options above and let the user pick; if dt is capped by something else
(hydro CFL, a source term, outputs), hesdirk2 buys accuracy only and costs 1.35x per step -- say so. Related:
[[rad-m1-design]].
