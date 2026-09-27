---
name: kokkos5-decision-0926
description: Kokkos 5.2.2 tested 2026-09-26 -- builds (C++20 + 5096-line DualView view_host/view_device rename), bitwise, but only 1-3 % faster on dhj and 0 % on the M1 He box (good node); decision (user) -- stay on 4.6.02; branch kokkos5 kept, not merged
metadata:
  type: project
---

**Measured on the good node hpc-sm-02-10, 3 interleaved reps each:**

| config | 4.6.02 | 5.2.2 | change |
|---|---|---|---|
| dhj nx1 256, 2 H200, median cycle | 16.63 ms | 16.25 ms | -2.3 % |
| dhj nx1 256, implicit-cycle p90 | | | -0.8 % |
| dhj nx1 76, 1 H200, median cycle | 12.23 ms | 11.54 ms | -5.6 % |
| dhj nx1 76, implicit-cycle p90 | | | -1.9 % |
| He box M1 full, 1 H200 | 68.3 ms | 69.7 ms | +2 % (noise) |

The earlier "-25 % box" was an artifact of the bad node hpc-sm-01-09.

**Changes Kokkos 5 needs** (branch `kokkos5`, commits 6b75b4fd and aec20106):
- C++20;
- DualView `h_view`/`d_view` renamed to `view_host()`/`view_device()` (5096 uses in 201 files);
- `HostMirror` renamed to `host_mirror_type`;
- `const ViewType&` in the multigrid templates;
- `<cfloat>` includes.
- GFX942_APU is still supported.

**Decision (user: "wait for the box result, then decide"):** stay on 4.6.02. Revisit only if a later Kokkos release gives a real gain or viper needs it. If revisited, do it as one isolated commit while no branches are open, with a HIP re-gate on viper.
