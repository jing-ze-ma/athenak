---
name: hip-dualview-sync-idiom
description: HIP build rejects DualView modify<Space>/sync<Space> templates; use modify_host/sync_device; and hipcc contraction breaks GPU bitwise identity when kernels are split
metadata:
  type: feedback
---

DualView sync must be written `modify_device(); sync_host();` / `modify_host(); sync_device();`,
never `template modify<DevExeSpace>()` / `template sync<HostMemSpace>()`: the templates fail a
Kokkos static_assert under hipcc (cdd7d2a5, again 27ca5b13). A CPU build accepts both, so code
written on orion breaks only when built on viper.

Also: splitting a kernel or moving arithmetic across inlined functions changes hipcc FMA
contraction; CPU stays bitwise, HIP differs at 1 ulp (77de5618, caad9247). Pure loop-invariant
hoists (1bd31298, 27ca5b13) stayed bitwise on HIP.

**Why:** the user asked (2026-09-11) that orion be told the convention so code written there
builds and stays bitwise on the GPU machine. Note pushed as
docs/handover/NOTE-2026-09-11-viper-gpu-code-conventions.md.

**How to apply:** use the space-free idiom in any new DualView code; when GPU bitwise identity
matters, verify on the GPU with an old-binary-twice control ([[validate-the-instrument]]).
