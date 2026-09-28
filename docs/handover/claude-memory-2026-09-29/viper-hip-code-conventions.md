---
name: viper-hip-code-conventions
description: "STANDING (viper note 09-11): GPU/HIP conventions for code written on orion: DualView modify_device()/sync_host() only (explicit-space templates fail under hipcc); hipcc changes round-off when arithmetic moves between kernels or inline functions, so CPU bit-identity does not imply GPU bit-identity; Debug build once for new View indexing"
metadata:
  type: feedback
---
From docs/handover/NOTE-2026-09-11-viper-gpu-code-conventions.md (fb69c554). Two things that pass
on orion CPU and break or change answers on viper HIP (ROCm 6.3.4, gfx942), each bitten twice:
1. DualView sync: NEVER `view.template modify<DevExeSpace>()` / `sync<HostMemSpace>()`; use
   `modify_device(); sync_host();` (or `modify_host(); sync_device();`). Hit in the hydro dt
   diagnostic (fixed cdd7d2a5) and the cs seam table (27ca5b13).
2. Splitting a kernel, or moving an expression into/out of an inline function, changes hipcc FMA
   contraction even with identical source: 1-ulp differences that grow chaotically in dhj runs.
   Hoisting loop-invariant work stays bitwise; kernel splits do not. A GPU last-bit difference
   after a refactor is expected; a run-to-run non-reproducible one is a race.
**Why:** the FOFC steps 3-5 factored GnomonicEquiangleRaiseVel(+MHD) bodies into inline headers and
split kernels; expect 1-ulp GPU drift there even though CPU is bit-identical.
**How to apply:** before calling GPU code done, build it on viper (or at least CUDA); state
"bitwise on CPU" explicitly and never claim GPU identity from it; Debug build once for new View
index schemes; copy lambda captures to locals. Related: [[viper-floors-legacy-gate]].
