---
name: rt-source-semi-implicit
description: "FIXED 048dff30: the two-stream RT source is now applied semi-implicitly (exponential relaxation at rate 4E/e); this cured the red giant dt collapse that rt_de_max's clamp only masked"
metadata:
  type: project
---

The two-stream source in `src/utils/two_stream_rt.hpp` used to be applied purely
explicitly, with `rt_de_max` clamping the overshoot. That clamp bounds damage but does
not stabilise anything: the clamped cell lags, its neighbour overshoots harder, and the
top of an atmosphere ends up oscillating cold-hot-cold.

**The fix (commit 048dff30).** Split `src = A - E(T)`, with `E = 4 sigma kappa_P rho T^4`
the cell's own emission. Both chain kernels now accumulate `E` into a new `rt_Em` array
as they sweep (the ck path is 17 % low because it carries `CK_DIFFUSIVITY = 1.66` where
the exact hemispheric factor is 2; immaterial for a rate estimate). The apply kernel then
takes `lambda = 4E/e` and `de = (src/lambda)*(1 - exp(-lambda*bdt))`.

**Why it is safe.** It is the identity when `lambda*bdt << 1`. A 200-cycle
`deep_hot_jupiter_rt` run is BITWISE IDENTICAL across the change -- that is the
regression test to repeat if this area is touched again (`inputs/mhd/
deep_hot_jupiter_rt_ideal_xe.athinput` with `nlim = 200`, diff the `.hst`).

**What it fixed.** The red giant cubed-sphere run lost its timestep from 32 s to 2e-3 s
between cycles 1400 and 1500, at t ~ 2.94e4 s. See [[red-giant-envelope-project]].

**Two traps that cost time here, worth remembering:**

1. The collapse looked like an inner-boundary problem because it appeared in the first
   run with `inner_bc = open`. It is NOT: a `wall` control collapsed at the SAME time to
   the same dt. Always run the control before believing a boundary is at fault.
2. The dt that collapsed was the CONDUCTION dt while hydro was still at 22 s. The cycle
   line does not say which module owns dt, so the same commit added a `### dt COLLAPSE`
   report to `Mesh::NewTimeStep` (fires on a >4x drop, only from the rank holding the
   global minimum, at most 20 lines). Use it first next time.
