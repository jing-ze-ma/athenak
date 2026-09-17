---
name: two-stream-rt-three-kernels-trap
description: "two_stream_rt.hpp has THREE mutually exclusive band-sweep kernels: rt_chain_grey (rt_grey), rt_chain_ck (rt_ck), rt_chain (else). A change to the sweep must be made in all three; red_giant uses rt_chain_grey. Cost a full round of the direct-source fix on 2026-09-09."
metadata:
  type: feedback
---

`picket_fence_two_stream_RT` in `src/utils/two_stream_rt.hpp` dispatches
`if (grey_on) -> "rt_chain_grey" (~1039, via launch_grey_chain)`,
`else if (ck_on) -> "rt_chain_ck" (~1164)`, `else -> "rt_chain" (~1389)`, and after the
band paths `return`s before the legacy picket-fence kernel (~1649).  They share structure
(down-sweep, up-sweep, e0 = -expm1(-x), alp/gm, bet, Fb_g += weight*(I_up - I_down),
Em_g) but are separate code.  I patched two of them and the one red_giant actually runs
(`rt_grey = true`) was untouched, so the new array read back as exact zeros -- and I then
spent a debugging round on "why is it 1e6 too small" when it was simply never written.

**Why:** grep for one distinctive line (`I_ir_up_c[cc] = (1.0-e0)*I_ir_up_c[cc]`) finds one
kernel; the grey one uses different variable names.  **How to apply:** before editing any
sweep, `grep -n 'par_for("rt_chain' two_stream_rt.hpp` and touch all three (and both
`#if RT_CACHE` variants); verify with `problem/rt_apply_debug` that the new quantity is
non-zero on the path the run takes.  Instrumentation that prints the raw accumulated
value is worth more than any amount of algebra.

Also recorded here: `<hydro>/wb_rmax` (per-cell radial cutoff of the well-balanced x1
reconstruction, src/hydro/hydro.hpp/.cpp/hydro_fluxes.cpp + the pgen gravity fallback)
was added and verified 2026-09-09: background free-falls at g to 5 digits like WB-off,
star bit-identical at t=0, provable no-op with the cutoff above the domain.
